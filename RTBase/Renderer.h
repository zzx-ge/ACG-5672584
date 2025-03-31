#pragma once

#include "Core.h"
#include "Sampling.h"
#include "Geometry.h"
#include "Imaging.h"
#include "Materials.h"
#include "Lights.h"
#include "Scene.h"
#include "GamesEngineeringBase.h"
#include <thread>
#include <functional>
#include <mutex>

#define MAX_DEPTH 8

class RayTracer
{
public:
	Scene* scene;
	GamesEngineeringBase::Window* canvas;
	Film* film;
	MTRandom *samplers;
	std::thread **threads;
	int numProcs;
	void init(Scene* _scene, GamesEngineeringBase::Window* _canvas)
	{
		scene = _scene;
		canvas = _canvas;
		film = new Film();
		film->init((unsigned int)scene->camera.width, (unsigned int)scene->camera.height, new BoxFilter());
		SYSTEM_INFO sysInfo;
		GetSystemInfo(&sysInfo);
		numProcs = sysInfo.dwNumberOfProcessors;
		threads = new std::thread*[numProcs];
		samplers = new MTRandom[numProcs];
		clear();
	}

	void clear()
	{
		film->clear();
	}
	Colour computeDirect(ShadingData shadingData, Sampler* sampler) {

		// 如果是纯镜面反射，直接返回0, 不参与直接光照计算
		if (shadingData.bsdf->isPureSpecular()) return Colour(0.0f);

		Colour result(0.0f);
		const float epsilon = 1e-3f;
		float allLightPdfSum = 0.f;
		
		// ========== 光源采样 ==========
		if (!scene->lights.empty()) {
			// 采样光源（包含PMF）
			float lightPMF;
			Light* light = scene->sampleLight(sampler, lightPMF); //check

			// 光源采样位置和方向
			Colour emitted;
			float lightPdf;
			Vec3 lightPos = light->sample(shadingData, sampler, emitted, lightPdf);

			if (lightPdf > 0 && !emitted.isZero()) {
				Vec3 wi = (lightPos - shadingData.x).normalize();

				// 可见性检测
				bool visible = scene->visible(shadingData.x, lightPos);
				//bool visible = true;
				if (visible) {
					// 几何项计算
					float distSq = (lightPos - shadingData.x).lengthSq();
					float cosTheta = Dot(wi, shadingData.sNormal);
					float cosLight = light->isArea() ? Dot(-wi, light->normal(shadingData, wi)) : 1.0f;
					float G = light->isArea() ? (cosTheta * cosLight) / distSq : cosTheta;

					// 面光源pdf转换为立体角pdf,后续其他光源pdf也需要转换
					if (light->isArea()) {
						float dist2 = (lightPos - shadingData.x).lengthSq();
						if (cosLight <= 1e-3f || dist2 <= 1e-6f) {
							//return Colour(0.f);
							lightPdf = 0.f;
						}
						lightPdf = lightPdf * dist2 / cosLight;
					}
					// Environment light
					else {
						lightPdf = light->PDF(shadingData, wi);
					}

					// 计算所有光源在wi方向的总pdf
					for (auto otherLight : scene->lights) {
						float pmf = scene->pmf(otherLight);
						float otherPdf = 0.f;
						if (otherLight->isArea()) {
							float otherCosLight = max(0.0f, Dot(-wi, otherLight->normal(shadingData, wi)));
							if (otherCosLight > 0) {
								float otherPdfArea = otherLight->PDF(shadingData, wi);
								otherPdf = otherPdfArea * distSq / otherCosLight;
							}
						}
						else {
							otherPdf = otherLight->PDF(shadingData, wi);
						}
						allLightPdfSum += pmf * otherPdf;
					}

					// 计算BSDF
					Colour bsdfVal = shadingData.bsdf->evaluate(shadingData, wi);
					float bsdfPdf = shadingData.bsdf->PDF(shadingData, wi);

					// MIS权重（考虑光源PMF）
					float misWeight = (lightPMF * lightPdf) / (allLightPdfSum + bsdfPdf + 1e-6f);
					//float misWeight = 1.f; // for debug
					result = result + bsdfVal * emitted * G * misWeight / (lightPMF * lightPdf + 1e-6f);
				}
			}
		}
		
		// ========== BSDF采样 ==========
		float bsdfPdf;
		Colour bsdfVal;
		Vec3 wi = shadingData.bsdf->sample(shadingData, sampler, bsdfVal, bsdfPdf);

		if (bsdfPdf > 0 && !bsdfVal.isZero()) {
			// 检测光线终点
			Ray visRay(shadingData.x + wi * epsilon, wi);
			IntersectionData lightIntersection = scene->traverse(visRay);

			// 获取辐射亮度
			Colour Le(0.0f);
			Vec3 lightpos;
			Vec3 n_light;
			bool envintersect = false;
			if (lightIntersection.t < FLT_MAX) {
				ShadingData lightShading = scene->calculateShadingData(lightIntersection, visRay);
				if (lightShading.bsdf->isLight()) {
					Le = lightShading.bsdf->emit();
					lightpos = lightShading.x;
					n_light = lightShading.sNormal;
				}
			}
			else{
				if (scene->background) {
					ShadingData lightShading;
					Le = scene->background->evaluate(lightShading, wi);
					envintersect = true;
				}
			}

			if (!Le.isZero()) {

				// MIS权重
				float misWeight = bsdfPdf / (bsdfPdf + allLightPdfSum + 1e-6f);
				float cosTheta = Dot(wi, shadingData.sNormal);
				result = result + bsdfVal * Le * cosTheta * misWeight / bsdfPdf;
			}
		}
		

		return result;
	}


	Colour estimateIndirectFromVPLs(const ShadingData& shadingData, Sampler* sampler) {
		Colour indirect(0.0f);

		for (const VPL& vpl : scene->vpls) {

			if (!scene->visible(shadingData.x + shadingData.sNormal * 1e-3f, vpl.position))
				continue;

			Vec3 wi = (vpl.position - shadingData.x).normalize();
			Vec3 wl = -wi;

			float dist2 = (vpl.position - shadingData.x).lengthSq();
			float cos1 = max(0.f, Dot(shadingData.sNormal, wi));
			float cos2 = max(0.f, Dot(vpl.normal, wl));
			float G = (cos1 * cos2) / (dist2 + 1e-6f);

			Colour f = shadingData.bsdf->evaluate(shadingData, wi);
			float pdf_bsdf = shadingData.bsdf->PDF(shadingData, wi);

			// Balance heuristic
			float weight = (pdf_bsdf) /
				(pdf_bsdf + vpl.pdf_light + 1e-6f);

			indirect = indirect + f * vpl.flux * G * weight;
		}

		return indirect / float(max(1, int(scene->vpls.size())));
	}

	/*
	std::vector<PathVertex> traceLight(Sampler* sampler) {
		std::vector<PathVertex> path;

		// 选择光源并采样初始光线
		float lightPdf;
		Light* light = scene->sampleLight(sampler, lightPdf);
		Vec3 lightPos = light->samplePositionFromLight(sampler, lightPdf);
		Vec3 lightDir = light->sampleDirectionFromLight(sampler, lightPdf);

		PathVertex vertex;
		ShadingData shadingdata;
		vertex.position = lightPos;
		vertex.normal = light->normal(shadingdata, lightPos);
		vertex.wi = lightDir;
		vertex.beta = light->totalIntegratedPower() * Dot(vertex.normal, lightDir) / lightPdf;
		vertex.pdf = lightPdf;
		vertex.bsdf = nullptr;
		path.push_back(vertex);

		// 传播光线并记录路径
		Ray ray(lightPos + lightDir * EPSILON, lightDir);
		for (int depth = 0; depth < 10; ++depth) {
			IntersectionData hit = scene->traverse(ray);
			if (hit.t == FLT_MAX) break;

			ShadingData shading = scene->calculateShadingData(hit, ray);
			vertex.position = shading.x;
			vertex.normal = shading.sNormal;
			vertex.bsdf = shading.bsdf;

			// 采样BSDF继续路径
			Vec3 wo = -ray.dir;
			Vec3 wi;
			Colour bsdfVal;
			float bsdfPdf;
			wi = shading.bsdf->sample(shading, sampler, bsdfVal, bsdfPdf);

			// 更新吞吐量和PDF
			vertex.beta = vertex.beta * bsdfVal * Dot(wi, shading.sNormal) / bsdfPdf;
			vertex.pdf = bsdfPdf;
			path.push_back(vertex);

			// 更新光线
			ray = Ray(shading.x + wi * EPSILON, wi);
		}
		return path;
	}
	*/
	/*
	void connectToCamera(const PathVertex& vertex) {
		Vec3 cameraPos = scene->camera.origin;
		Vec3 dirToCamera = (cameraPos - vertex.position).normalize();

		if (!scene->visible(vertex.position, cameraPos)) return;

		float u, v;
		if (!scene->camera.projectOntoCamera(vertex.position, u, v)) return;

		// 计算路径贡献
		Colour bsdfVal = vertex.bsdf ? vertex.bsdf->evaluate();
		float cosTheta = Dot(dirToCamera, vertex.normal);

	}
	*/

	// === Light Tracing integrated with MIS Path Tracing ===
	void traceLightPaths(Sampler* sampler, int numPaths, int maxDepth) {

		for (int pathIdx = 0; pathIdx < numPaths; ++pathIdx) {
			// Step 1: 光源采样起点
			float lightPdf;
			Light* light = scene->sampleLight(sampler, lightPdf);
			Vec3 lightPos = light->samplePositionFromLight(sampler, lightPdf);
			float dirPdf;
			Vec3 wi = light->sampleDirectionFromLight(sampler, dirPdf);

			// 将lightPdf转化为solid angle pdf, 距离先为1
			float cosThetaLight = max(0.f, Dot(wi, light->normal(ShadingData(), lightPos)));
			lightPdf = lightPdf / max(cosThetaLight, EPSILON);
			
			// 初始化路径贡献 (beta)
			float Pmf = scene->pmf(light);
			Colour beta = light->evaluate(ShadingData(), wi) * cosThetaLight / (dirPdf * Pmf);


			// 初始光线
			Ray ray(lightPos + wi * EPSILON, wi);

			// Step 2: 路径追踪
			for (int depth = 0; depth < maxDepth; ++depth) {
				IntersectionData hit = scene->traverse(ray);
				if (hit.t == FLT_MAX) break;

				ShadingData shading = scene->calculateShadingData(hit, ray);

				// 得到真正的lightPdf
				lightPdf *= (lightPos - shading.x).lengthSq();

				// Step 3: 连接到相机 (Pinhole相机)
				Vec3 dirToCam = (scene->camera.origin - shading.x).normalize();
				if (scene->visible(shading.x + shading.sNormal * EPSILON, scene->camera.origin)) {
					// 投影到相机
					float u, v;
					if (scene->camera.projectOntoCamera(shading.x, u, v)) {
						Colour bsdfVal = shading.bsdf->evaluate(shading, dirToCam);
						float cosTheta = max(0.0f, Dot(shading.sNormal, dirToCam));

						// MIS权重计算
						float bsdfPdf = shading.bsdf->PDF(shading, dirToCam);
						float allLightPdfSum = 0.0f;
						for (auto otherLight : scene->lights) {
							float pmf = scene->pmf(otherLight);
							float otherPdf = 0.f;
							if (otherLight->isArea()) {
								AreaLight* areaLight = static_cast<AreaLight*>(otherLight);
								Vec3 samplePoint = areaLight->triangle->sample(sampler, otherPdf);
								float distSq = (samplePoint - shading.x).lengthSq();
								float cosLight = max(0.0f, Dot(-dirToCam, areaLight->triangle->gNormal()));
								if (cosLight > EPSILON) {
									otherPdf = (1.0f / areaLight->triangle->area) * distSq / cosLight;
								}
							}
							else {
								otherPdf = otherLight->PDF(shading, dirToCam);
							}
							allLightPdfSum += pmf * otherPdf;
						}
						float misWeight = (dirPdf * Pmf) / (allLightPdfSum + bsdfPdf + 1e-6f);

						// 累加贡献到图像
						Colour contrib = beta * bsdfVal * cosTheta * misWeight;
						film->splat(u, v, beta * bsdfVal * cosTheta * misWeight);
					}
				}

				// Step 4: 继续传播路径
				Colour bsdfVal;
				float bsdfPdf;
				Vec3 newWi = shading.bsdf->sample(shading, sampler, bsdfVal, bsdfPdf);
				if (bsdfPdf < EPSILON || bsdfVal.isZero()) break;

				// 更新路径贡献
				float cosTheta = max(0.0f, Dot(newWi, shading.sNormal));
				beta = beta * bsdfVal * cosTheta / bsdfPdf;

				// Russian Roulette (optional)
				if (depth > 3) {
					float rr = min(0.9f, beta.maxComponent());
					if (sampler->next() > rr) break;
					beta = beta / rr;
				}

				// 更新光线
				ray = Ray(shading.x + newWi * EPSILON, newWi);
			}
		}
	}


	Colour pathTrace(Ray& r, Colour& pathThroughput, int depth, Sampler* sampler) {
		if (depth >= MAX_DEPTH) return Colour(0.0f);

		IntersectionData intersection = scene->traverse(r);
		ShadingData shadingData = scene->calculateShadingData(intersection, r);

		if (shadingData.t < FLT_MAX) {
			// 直接处理光源贡献，无需条件判断
			if (shadingData.bsdf->isLight()) {
				return pathThroughput * shadingData.bsdf->emit();
			}

			// ----------- 直接光照 -----------
			Colour direct = pathThroughput * computeDirect(shadingData, sampler);

			// ----------- 间接光照 -----------
			Colour indirectVPL = pathThroughput * estimateIndirectFromVPLs(shadingData, sampler);

			
			// ----------- 俄罗斯轮盘赌 -----------
			if (depth > 3) {
				float rrProb = min(0.95f, pathThroughput.Lum());
				if (sampler->next() > rrProb) return direct + indirectVPL;
				pathThroughput = pathThroughput / rrProb;
			}

			// ----------- 采样BSDF -----------
			float pdf;
			Colour bsdfVal;
			Vec3 wi = shadingData.bsdf->sample(shadingData, sampler, bsdfVal, pdf);

			if (pdf > 0.0f && !bsdfVal.isZero()) {
				float cosTheta = Dot(wi, shadingData.sNormal);
				Colour newThroughput = pathThroughput * (bsdfVal * cosTheta) / pdf;
				Ray newRay(shadingData.x + wi * EPSILON, wi);
				return direct + indirectVPL + pathTrace(newRay, newThroughput, depth + 1, sampler);
				//return indirectVPL + pathTrace(newRay, newThroughput, depth + 1, sampler);
			}
			return direct + indirectVPL;
			//return indirectVPL; // test radiuosity
		}

		// 背景处理
		ShadingData shadingdata;
		return scene->background ? pathThroughput * scene->background->evaluate(shadingdata, r.dir) : Colour(0.0f);
	}

	Colour direct(Ray& r, Sampler* sampler)
	{
		IntersectionData intersection = scene->traverse(r);
		ShadingData shadingData = scene->calculateShadingData(intersection, r);
		if (shadingData.t < FLT_MAX)
		{
			if (shadingData.bsdf->isLight())
			{
				return shadingData.bsdf->emit(shadingData, shadingData.wo);
			}
			return computeDirect(shadingData, sampler);
		}
		return Colour(0.0f, 0.0f, 0.0f);
	}
	Colour albedo(Ray& r)
	{
		IntersectionData intersection = scene->traverse(r);
		ShadingData shadingData = scene->calculateShadingData(intersection, r);
		if (shadingData.t < FLT_MAX)
		{
			if (shadingData.bsdf->isLight())
			{
				return shadingData.bsdf->emit(shadingData, shadingData.wo);
			}
			return shadingData.bsdf->evaluate(shadingData, Vec3(0, 1, 0));
		}
		return scene->background->evaluate(shadingData, r.dir);
	}
	Colour viewNormals(Ray& r)
	{
		IntersectionData intersection = scene->traverse(r);
		if (intersection.t < FLT_MAX)
		{
			ShadingData shadingData = scene->calculateShadingData(intersection, r);
			return Colour(fabsf(shadingData.sNormal.x), fabsf(shadingData.sNormal.y), fabsf(shadingData.sNormal.z));
		}
		return Colour(0.0f, 0.0f, 0.0f);
	}


	void render() {
		film->incrementSPP();

		// Get CPU core count for multi-threading
		int numThreads = std::thread::hardware_concurrency();

		std::vector<std::thread> workers;
		std::mutex filmMutex;  // Mutex for film updates

		
		// Lambda function for each thread to process part of the image
		auto renderTile = [&](int threadID, int startY, int endY) {
			MTRandom sampler(threadID); // Unique sampler per thread

			for (int y = startY; y < endY; y++) {
				for (int x = 0; x < film->width; x++) {
					Colour finalColor(0.f, 0.f, 0.f);
					Colour Albedo(0.f, 0.f, 0.f);
					Colour Viewnormal(0.f, 0.f, 0.f);
					// Path tracing: multiple samples per pixel (anti-aliasing)
					int numSamples = 16; // Adjust for quality
					for (int s = 0; s < numSamples; s++) {
						// Generate primary ray
						float px = x + sampler.next();
						float py = y + sampler.next();
						Ray ray = scene->camera.generateRay(px, py);

						// Path tracing
						Colour pathThroughput(1.f, 1.f, 1.f);
						finalColor = finalColor + pathTrace(ray, pathThroughput, 0, &sampler) + finalColor;
						Albedo = Albedo + albedo(ray);
						Viewnormal = Viewnormal + viewNormals(ray);
						Albedo = albedo(ray);
						Viewnormal = viewNormals(ray);
					}

					// Average color over samples
					finalColor = finalColor / float(numSamples);
					Albedo = Albedo / float(numSamples);
					Viewnormal = Viewnormal / float(numSamples);
					// Store final pixel color (Thread-Safe)
					{
						std::lock_guard<std::mutex> lock(filmMutex); // Protect shared memory
						film->splat(x, y, finalColor);
						film->AOV(x, y, Albedo, Vec3(Viewnormal.r, Viewnormal.g, Viewnormal.b));
					}

					// Tonemapping and gamma correction
					unsigned char r, g, b;
					film->tonemap(finalColor.r, finalColor.g, finalColor.b, r, g, b);
					{
						std::lock_guard<std::mutex> lock(filmMutex); // Ensure safe canvas updates
						canvas->draw(x, y, r, g, b);
					}
				}
			}
		};

		
		// Multi-threading setup with better load balancing
		int tileHeight = max(1, film->height / numThreads);
		for (int i = 0; i < numThreads; i++) {
			int startY = i * tileHeight;
			int endY = min((i + 1) * tileHeight, film->height);
			workers.emplace_back(renderTile, i, startY, endY);
		}

		// Join all threads
		for (auto& worker : workers) {
			worker.join();
		}
		
		/*

		// 光源路径追踪（同样是多线程实现）
		auto lightTraceWorker = [&](int pathsPerThread, int threadID) {
			MTRandom sampler(threadID + 1000); // 避免与相机采样重复
			traceLightPaths(&sampler, pathsPerThread, MAX_DEPTH);
		};

		workers.clear();
		int totalPaths = film->width * film->height;
		int pathsPerThread = totalPaths / numThreads;

		for (int i = 0; i < numThreads; ++i) {
			workers.emplace_back(lightTraceWorker, pathsPerThread, i);
		}

		// 等待光源路径追踪完成
		for (auto& worker : workers) worker.join();
		*/
	}

	int getSPP()
	{
		return film->SPP;
	}
	void saveHDR(std::string filename)
	{
		film->save(filename);
	}
	void savePNG(std::string filename)
	{
		stbi_write_png(filename.c_str(), canvas->getWidth(), canvas->getHeight(), 3, canvas->getBackBuffer(), canvas->getWidth() * 3);
	}
};