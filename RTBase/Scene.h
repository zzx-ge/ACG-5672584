#pragma once

#include "Core.h"
#include "Sampling.h"
#include "Geometry.h"
#include "Imaging.h"
#include "Materials.h"
#include "Lights.h"

class Camera
{
public:
	Matrix projectionMatrix;
	Matrix inverseProjectionMatrix;
	Matrix camera;
	Matrix cameraToView;
	float width = 0;
	float height = 0;
	Vec3 origin;
	Vec3 viewDirection;
	float Afilm;

	void init(Matrix ProjectionMatrix, int screenwidth, int screenheight)
	{
		projectionMatrix = ProjectionMatrix;
		inverseProjectionMatrix = ProjectionMatrix.invert();
		width = (float)screenwidth;
		height = (float)screenheight;
		float Wlens = (2.0f / ProjectionMatrix.a[1][1]);
		float aspect = ProjectionMatrix.a[0][0] / ProjectionMatrix.a[1][1];
		float Hlens = Wlens * aspect;
		Afilm = Wlens * Hlens;
	}
	void updateView(Matrix V)
	{
		camera = V;
		cameraToView = V.invert();
		origin = camera.mulPoint(Vec3(0, 0, 0));
		viewDirection = inverseProjectionMatrix.mulPointAndPerspectiveDivide(Vec3(0, 0, 1));
		viewDirection = camera.mulVec(viewDirection);
		viewDirection = viewDirection.normalize();
	}
	// Add code here
	Ray generateRay(float x, float y)
	{
		// Step 1: Convert to Normalized Device Coordinates (NDC)
		float ndcX = (2.0f * x / width) - 1.0f;
		float ndcY = 1.0f - (2.0f * y / height); // Flip Y-axis

		// Step 2: Convert to Camera Space (Clip Space -> Camera Space)
		Vec3 clipCoords(ndcX, ndcY, -1.0f, 1.0f);
		Vec3 cameraCoords = inverseProjectionMatrix.mulPoint(clipCoords);

		// Step 3: Perspective Divide to Get 3D Point
		Vec3 pointCamera = Vec3(cameraCoords.x / cameraCoords.w,
			cameraCoords.y / cameraCoords.w,
			cameraCoords.z / cameraCoords.w);

		// Step 4: Convert to World Space
		Vec3 pointWorld = camera.mulPoint(Vec3(pointCamera.x, pointCamera.y, pointCamera.z, 1.0f));

		// Step 5: Compute Ray Direction
		Vec3 dir = (pointWorld - origin).normalize();

		return Ray(origin, dir);
	}
	bool projectOntoCamera(const Vec3& p, float& x, float& y)
	{
		Vec3 pview = cameraToView.mulPoint(p);
		Vec3 pproj = projectionMatrix.mulPointAndPerspectiveDivide(pview);
		x = (pproj.x + 1.0f) * 0.5f;
		y = (pproj.y + 1.0f) * 0.5f;
		if (x < 0 || x > 1.0f || y < 0 || y > 1.0f)
		{
			return false;
		}
		x = x * width;
		y = 1.0f - y;
		y = y * height;
		return true;
	}
};

class Scene
{
public:
	std::vector<Triangle> triangles;
	std::vector<int> triangleIndices;
	std::vector<BSDF*> materials;
	std::vector<Light*> lights;
	Light* background = NULL;
	BVHNode* bvh = NULL;
	Camera camera;
	AABB bounds;
	std::vector<VPL> vpls;

	void build()
	{
		// Add BVH building code here
		if (bvh == nullptr) {
			bvh = new BVHNode();
		}
		//bvh->build(triangles);
		for (int i = 0; i < triangles.size(); i++) {
			triangleIndices.push_back(i);
		}
		bvh->build(triangleIndices, triangles);

		// Do not touch the code below this line!
		// Build light list
		for (int i = 0; i < triangles.size(); i++)
		{
			if (materials[triangles[i].materialIndex]->isLight())
			{
				AreaLight* light = new AreaLight();
				light->triangle = &triangles[i];
				light->emission = materials[triangles[i].materialIndex]->emission;
				lights.push_back(light);
			}
		}

		for (int i = 0; i < materials.size(); ++i) {
			if (materials[i]->isLight()) {
				std::cout << "Material " << i << " is light. Emission: "
					<< materials[i]->emission.r << ", "
					<< materials[i]->emission.g << ", "
					<< materials[i]->emission.b << std::endl;
			}
		}

		int numVPLs = 1000;
		MTRandom initSampler(42);
		std::vector<VPL>vpls = generateVPLs(numVPLs, &initSampler);
	}

	IntersectionData traverse(const Ray& ray)
	{
		IntersectionData intersection;
		intersection.t = FLT_MAX;
		if (bvh) {
			bvh->traverse(ray, triangles, intersection);
		}

		else {
			for (int i = 0; i < triangles.size(); i++)
			{
				float t;
				float u;
				float v;
				if (triangles[i].rayIntersect(ray, t, u, v))
				{
					if (t < intersection.t)
					{
						intersection.t = t;
						intersection.ID = i;
						intersection.alpha = u;
						intersection.beta = v;
						intersection.gamma = 1.0f - (u + v);
					}
				}
			}
		}
		return intersection;
	}

	Light* sampleLight(Sampler* sampler, float& pmf) {
		// 合并所有光源（包含背景光）
		std::vector<Light*> allLights = lights;
		if (background && std::find(allLights.begin(), allLights.end(), background) == allLights.end()) {
			allLights.push_back(background);
		}

		if (allLights.empty()) {
			pmf = 0.0f;
			return nullptr;
		}

		// 根据光源功率计算采样权重
		std::vector<float> weights;
		float totalPower = 0.0f;
		for (auto light : allLights) {
			float power = light->totalIntegratedPower();
			weights.push_back(power);
			totalPower += power;
		}

		// 处理全零功率（如只有背景光）
		if (totalPower <= 0.0f) {
			pmf = 1.0f / allLights.size();
			int idx = sampler->next() * allLights.size();
			return allLights[idx];
		}

		// 使用功率分布采样
		float xi = sampler->next() * totalPower;
		float accum = 0.0f;
		for (int i = 0; i < allLights.size(); ++i) {
			accum += weights[i];
			if (xi <= accum) {
				pmf = weights[i] / totalPower;
				return allLights[i];
			}
		}

		pmf = weights.back() / totalPower;
		return allLights.back();
	}

	// Do not modify any code below this line
	void init(std::vector<Triangle> meshTriangles, std::vector<BSDF*> meshMaterials, Light* _background)
	{
		for (int i = 0; i < meshTriangles.size(); i++)
		{
			triangles.push_back(meshTriangles[i]);
			bounds.extend(meshTriangles[i].vertices[0].p);
			bounds.extend(meshTriangles[i].vertices[1].p);
			bounds.extend(meshTriangles[i].vertices[2].p);
		}
		for (int i = 0; i < meshMaterials.size(); i++)
		{
			materials.push_back(meshMaterials[i]);
		}
		background = _background;
		if (background->totalIntegratedPower() > 0)
		{
			lights.push_back(background);
		}
	}

	bool visible(const Vec3& p1, const Vec3& p2)
	{
		Ray ray;
		Vec3 dir = p2 - p1;
		float maxT = dir.length() - (2.0f * EPSILON);
		dir = dir.normalize();
		ray.init(p1 + (dir * EPSILON), dir);
		return bvh->traverseVisible(ray, triangles, maxT);
	}

	Colour emit(Triangle* light, ShadingData shadingData, Vec3 wi)
	{
		return materials[light->materialIndex]->emit(shadingData, wi);
	}

	ShadingData calculateShadingData(IntersectionData intersection, Ray& ray)
	{
		ShadingData shadingData = {};
		if (intersection.t < FLT_MAX)
		{
			shadingData.x = ray.at(intersection.t);
			shadingData.gNormal = triangles[intersection.ID].gNormal();
			triangles[intersection.ID].interpolateAttributes(intersection.alpha, intersection.beta, intersection.gamma, shadingData.sNormal, shadingData.tu, shadingData.tv);
			shadingData.bsdf = materials[triangles[intersection.ID].materialIndex];
			shadingData.wo = -ray.dir;
			if (shadingData.bsdf->isTwoSided())
			{
				if (Dot(shadingData.wo, shadingData.sNormal) < 0)
				{
					shadingData.sNormal = -shadingData.sNormal;
				}
				if (Dot(shadingData.wo, shadingData.gNormal) < 0)
				{
					shadingData.gNormal = -shadingData.gNormal;
				}
			}
			shadingData.frame.fromVector(shadingData.sNormal);
			shadingData.t = intersection.t;
		} else
		{
			shadingData.wo = -ray.dir;
			shadingData.t = intersection.t;
		}
		return shadingData;
	}

	std::vector<VPL> generateVPLs(int numPaths, Sampler* sampler) {
		std::vector<VPL> vpls;

		for (int i = 0; i < numPaths; ++i) {
			// 采样光源
			float pmf;
			Light* light = this->sampleLight(sampler, pmf);
			if (light->isArea()) {
				Colour Le;
				Vec3 origin, normal;
				float lightPdf;
				float dirPdf;

				// 从光源上采样一个点
				origin = light->samplePositionFromLight(sampler, lightPdf);
				ShadingData sd;
				normal = light->normal(sd, origin);

				// 构造 dummy shadingData（用于方向采样）
				ShadingData dummy;
				dummy.x = origin;
				dummy.sNormal = normal; // 任意单位向上法线
				dummy.wo = -normal;     // 任意方向（可选为 -normal）
				dummy.frame.fromVector(dummy.sNormal);

				// 从光源上采样出射方向
				Vec3 wi = light->sample(dummy, sampler, Le, dirPdf);

				// 发射光线
				Ray ray(origin + normal * EPSILON, wi);
				IntersectionData intersection = this->traverse(ray);
				ShadingData shadingData = this->calculateShadingData(intersection, ray);

				if (intersection.t < FLT_MAX && !shadingData.bsdf->isLight() && !shadingData.bsdf->isPureSpecular()) {
					Colour fr = shadingData.bsdf->evaluate(shadingData, -wi);
					float dist2 = (shadingData.x - origin).lengthSq();
					float cosLight = std::max(Dot(normal, wi), 0.f);
					if (cosLight > 1e-6f && dist2 > 1e-6f) {
						dirPdf = dirPdf * dist2 / cosLight;
						lightPdf = lightPdf * dist2 / cosLight;
					}
					else {
						dirPdf = 0.f;
						lightPdf = 0.f;
					}
					float cosTheta = std::max(0.f, Dot(shadingData.sNormal, -wi));
					float G = (cosLight * cosTheta) / dist2;

					VPL vpl;
					vpl.position = shadingData.x;
					vpl.normal = shadingData.sNormal;
					vpl.bsdf = shadingData.bsdf;
					vpl.wi = -wi;
					vpl.pdf_light = lightPdf * pmf * dirPdf;
					vpl.flux = Le * fr * G / (vpl.pdf_light + 1e-6f);

					vpls.push_back(vpl);
				}
			}	
		}

		return vpls;
	}

	float pmf(Light* queryLight) {
		// 合并所有光源（包含背景光）
		std::vector<Light*> allLights = lights;
		if (background && std::find(allLights.begin(), allLights.end(), background) == allLights.end()) {
			allLights.push_back(background);
		}

		// 计算总功率
		float totalPower = 0.0f;
		for (auto light : allLights) {
			totalPower += light->totalIntegratedPower();
		}

		// 处理全零功率（回退到均匀分布）
		if (totalPower <= 0.0f) {
			return 1.0f / allLights.size();
		}

		// 返回查询光源的PMF
		return queryLight->totalIntegratedPower() / totalPower;
	}
};