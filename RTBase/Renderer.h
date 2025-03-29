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
#include "shadermanager.h"

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

		// ========== 光源采样 ==========
		if (!scene->lights.empty()) {
			// 采样光源（包含PMF）
			float lightPMF;
			Light* light = scene->sampleLight(sampler, lightPMF);

			// 光源采样位置和方向
			Colour emitted;
			float lightPdf;
			Vec3 lightPos = light->sample(shadingData, sampler, emitted, lightPdf);

			if (lightPdf > 0 && !emitted.isZero()) {
				Vec3 wi = (lightPos - shadingData.x).normalize();

				// 可见性检测
				Ray shadowRay(shadingData.x + wi * epsilon, wi);
				float tMax = (lightPos - shadingData.x).length() - 2 * epsilon;
				bool visible = scene->traverse(shadowRay).t >= tMax;

				if (visible) {
					// 几何项计算
					float cosTheta = Dot(wi, shadingData.sNormal);
					float cosLight = light->isArea() ? Dot(-wi, light->normal(shadingData, wi)) : 1.0f;
					float G = light->isArea() ? (cosTheta * cosLight) / (lightPos - shadingData.x).lengthSq() : cosTheta;

					// 计算BSDF
					Colour bsdfVal = shadingData.bsdf->evaluate(shadingData, wi);
					float bsdfPdf = shadingData.bsdf->PDF(shadingData, wi);

					// MIS权重（考虑光源PMF）
					float misWeight = (lightPMF * lightPdf) / (lightPMF * lightPdf + bsdfPdf + 1e-6f);
					result = result + bsdfVal * emitted * G * misWeight / (lightPMF * lightPdf);
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
			if (lightIntersection.t < FLT_MAX) {
				ShadingData lightShading = scene->calculateShadingData(lightIntersection, visRay);
				if (lightShading.bsdf->isLight()) {
					Le = lightShading.bsdf->emit(lightShading, -wi);
				}
			}
			else if (scene->background) {
				ShadingData lightShading;
				Le = scene->background->evaluate(lightShading, wi);
			}

			if (!Le.isZero()) {
				// 计算光源PDF（包含PMF）
				float lightPdfSum = 0.0f;
				for (auto light : scene->lights) {
					float pmf = scene->pmf(light);
					lightPdfSum += pmf * light->PDF(shadingData, wi);
				}

				if (scene->background) {
					float bgpmf = scene->pmf(scene->background);
					lightPdfSum += bgpmf * scene->background->PDF(shadingData, wi);
				}

				// MIS权重
				float misWeight = bsdfPdf / (bsdfPdf + lightPdfSum + 1e-6f);
				float cosTheta = Dot(wi, shadingData.sNormal);
				result = result + bsdfVal * Le * cosTheta * misWeight / bsdfPdf;
			}
		}

		return result;
	}


	Colour estimateIndirectFromVPLs(const ShadingData& shadingData, Sampler* sampler) {
		Colour indirect(0.0f);

		for (const VPL& vpl : scene->vpls) {
			if (!scene->visible(shadingData.x, vpl.position))
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
			float weight = (pdf_bsdf * pdf_bsdf) /
				(pdf_bsdf * pdf_bsdf + vpl.pdf_light * vpl.pdf_light + 1e-6f);

			indirect = indirect + f * vpl.flux * G * weight;
		}

		return indirect / float(max(1, int(scene->vpls.size())));
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
				return direct + indirectVPL +
					pathTrace(newRay, newThroughput, depth + 1, sampler);
			}
			return direct + indirectVPL;
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

		ShaderManager shadermanager;

		std::vector<std::thread> workers;
		std::mutex filmMutex;  // Mutex for film updates

		// Lambda function for each thread to process part of the image
		auto renderTile = [&](int threadID, int startY, int endY) {
			MTRandom sampler(threadID); // Unique sampler per thread

			for (int y = startY; y < endY; y++) {
				for (int x = 0; x < film->width; x++) {
					Colour finalColor(0.f, 0.f, 0.f);

					// Path tracing: multiple samples per pixel (anti-aliasing)
					int numSamples = 16; // Adjust for quality
					for (int s = 0; s < numSamples; s++) {
						// Generate primary ray
						float px = x + sampler.next();
						float py = y + sampler.next();
						Ray ray = scene->camera.generateRay(px, py);

						// Path tracing
						Colour pathThroughput(1.f, 1.f, 1.f);
						finalColor = pathTrace(ray, pathThroughput, 0, &sampler) + finalColor;
					}

					// Average color over samples
					finalColor = finalColor / float(numSamples);

					// Store final pixel color (Thread-Safe)
					{
						std::lock_guard<std::mutex> lock(filmMutex); // Protect shared memory
						film->splat(x, y, finalColor);
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