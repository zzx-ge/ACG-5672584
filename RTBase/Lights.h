#pragma once

#include "Core.h"
#include "Geometry.h"
#include "Materials.h"
#include "Sampling.h"

#pragma warning( disable : 4244)

struct VPL {
	Vec3 position;
	Vec3 normal;
	BSDF* bsdf;      // Material at the VPL
	Colour flux;     // Le * fr * cos / lightPdf
	Vec3 wi;         // Direction from VPL to the surface (shading point)
	float pdf_light; // Light sampling PDF
};

class SceneBounds
{
public:
	Vec3 sceneCentre;
	float sceneRadius;
};

class Light
{
public:
	virtual Vec3 sample(const ShadingData& shadingData, Sampler* sampler, Colour& emittedColour, float& pdf) = 0;
	virtual Colour evaluate(const ShadingData& shadingData, const Vec3& wi) = 0;
	virtual float PDF(const ShadingData& shadingData, const Vec3& wi) = 0;
	virtual bool isArea() = 0;
	virtual Vec3 normal(const ShadingData& shadingData, const Vec3& wi) = 0;
	virtual float totalIntegratedPower() = 0;
	virtual Vec3 samplePositionFromLight(Sampler* sampler, float& pdf) = 0;
	virtual Vec3 sampleDirectionFromLight(Sampler* sampler, float& pdf) = 0;
};

class AreaLight : public Light
{
public:
	Triangle* triangle = NULL;
	Colour emission;
	Vec3 lastSamplePoint;
	Vec3 sample(const ShadingData& shadingData, Sampler* sampler, Colour& emittedColour, float& pdf)
	{
		emittedColour = emission;
		lastSamplePoint = triangle->sample(sampler, pdf);
		return lastSamplePoint;
	}
	Colour evaluate(const ShadingData& shadingData, const Vec3& wi)
	{
		if (Dot(wi, triangle->gNormal()) < 0)
		{
			return emission;
		}
		return Colour(0.0f, 0.0f, 0.0f);
	}
	float PDF(const ShadingData& shadingData, const Vec3& wi) {

		return 1.0f / triangle->area;

	}
	bool isArea()
	{
		return true;
	}
	Vec3 normal(const ShadingData& shadingData, const Vec3& wi)
	{
		return triangle->gNormal();
	}
	float totalIntegratedPower()
	{
		return (triangle->area * emission.Lum());
	}
	Vec3 samplePositionFromLight(Sampler* sampler, float& pdf)
	{
		return triangle->sample(sampler, pdf);
	}
	Vec3 sampleDirectionFromLight(Sampler* sampler, float& pdf)
	{
		// Add code to sample a direction from the light
		Vec3 wi = Vec3(0, 0, 1);
		pdf = 1.0f;
		Frame frame;
		frame.fromVector(triangle->gNormal());
		return frame.toWorld(wi);
	}
};

class BackgroundColour : public Light
{
public:
	Colour emission;
	BackgroundColour(Colour _emission)
	{
		emission = _emission;
	}
	Vec3 sample(const ShadingData& shadingData, Sampler* sampler, Colour& reflectedColour, float& pdf)
	{
		Vec3 wi = SamplingDistributions::uniformSampleSphere(sampler->next(), sampler->next());
		pdf = SamplingDistributions::uniformSpherePDF(wi);
		reflectedColour = emission;
		return wi;
	}
	Colour evaluate(const ShadingData& shadingData, const Vec3& wi)
	{
		return emission;
	}
	float PDF(const ShadingData& shadingData, const Vec3& wi) override
	{
		return SamplingDistributions::uniformSpherePDF(wi);
	}
	bool isArea() override
	{
		return false;
	}
	Vec3 normal(const ShadingData& shadingData, const Vec3& wi)
	{
		return -wi;
	}
	float totalIntegratedPower()
	{
		return emission.Lum() * 4.0f * M_PI;
	}
	Vec3 samplePositionFromLight(Sampler* sampler, float& pdf)
	{
		Vec3 p = SamplingDistributions::uniformSampleSphere(sampler->next(), sampler->next());
		p = p * use<SceneBounds>().sceneRadius;
		p = p + use<SceneBounds>().sceneCentre;
		pdf = 4 * M_PI * use<SceneBounds>().sceneRadius * use<SceneBounds>().sceneRadius;
		return p;
	}
	Vec3 sampleDirectionFromLight(Sampler* sampler, float& pdf)
	{
		Vec3 wi = SamplingDistributions::uniformSampleSphere(sampler->next(), sampler->next());
		pdf = SamplingDistributions::uniformSpherePDF(wi);
		return wi;
	}
};

class EnvironmentMap : public Light
{
public:
	Texture* env;
	Distribution2D* distribution;
	EnvironmentMap(Texture* _env)
	{
		env = _env;
		std::vector<std::vector<float>> luminanceMap(env->height, std::vector<float>(env->width));
		for (int y = 0; y < env->height; ++y) {
			for (int x = 0; x < env->width; ++x) {
				Colour c = env->texels[y * env->width + x];
				float lum = c.Lum();

				// 经纬度重要性采样需要加 sin(θ) 权重
				float theta = M_PI * (float(y) + 0.5f) / float(env->height);
				luminanceMap[y][x] = lum * sinf(theta);
			}
		}

		distribution = new Distribution2D(luminanceMap);
	}
	Vec3 sample(const ShadingData& shadingData, Sampler* sampler, Colour& reflectedColour, float& pdf) override {
		int x, y;
		float u = sampler->next();
		float v = sampler->next();
		distribution->sample(u, v, x, y, pdf);

		// 转为经纬角度
		float uCoord = (x + 0.5f) / float(env->width);
		float vCoord = (y + 0.5f) / float(env->height);

		float phi = uCoord * 2.0f * M_PI;
		float theta = vCoord * M_PI;

		//achor
		pdf = pdf / (2.0f * M_PI * M_PI * sinf(theta));

		// 方向向量
		Vec3 wi = Vec3(
			sin(theta) * cos(phi),
			cos(theta),
			sin(theta) * sin(phi)
		);

		reflectedColour = evaluate(shadingData, wi);
		return wi;
	}
	Colour evaluate(const ShadingData& shadingData, const Vec3& wi) override
	{
		float u = atan2f(wi.z, wi.x);
		u = (u < 0.0f) ? u + (2.0f * M_PI) : u;
		u = u / (2.0f * M_PI);
		float v = acosf(wi.y) / M_PI;
		return env->sample(u, v);
	}
	float PDF(const ShadingData& shadingData, const Vec3& wi) override {
		float theta = acosf(clamp(wi.y, -1.0f, 1.0f));
		float phi = atan2f(wi.z, wi.x);
		if (phi < 0.0f) phi += 2.0f * M_PI;

		float u = phi / (2.0f * M_PI);
		float v = theta / M_PI;

		int x = clamp(int(u * env->width), 0, env->width - 1);
		int y = clamp(int(v * env->height), 0, env->height - 1);

		float sinTheta = sinf(theta);
		if (sinTheta == 0.0f) return 0.0f;

		float marginalPdf, conditionalPdf;
		marginalPdf = distribution->marginal.func[y] / distribution->marginal.integral;
		conditionalPdf = distribution->conditional[y].func[x] / distribution->conditional[y].integral;

		return (marginalPdf * conditionalPdf) / (2.0f * M_PI * M_PI * sinTheta);
	}
	bool isArea()
	{
		return false;
	}
	Vec3 normal(const ShadingData& shadingData, const Vec3& wi) override
	{
		return -wi;
	}
	float totalIntegratedPower() override
	{
		float total = 0;
		for (int i = 0; i < env->height; i++)
		{
			float st = sinf(((float)i / (float)env->height) * M_PI);
			for (int n = 0; n < env->width; n++)
			{
				total += (env->texels[(i * env->width) + n].Lum() * st);
			}
		}
		total = total / (float)(env->width * env->height);
		return total * 4.0f * M_PI;
	}
	Vec3 samplePositionFromLight(Sampler* sampler, float& pdf)
	{
		// Samples a point on the bounding sphere of the scene. Feel free to improve this.
		Vec3 p = SamplingDistributions::uniformSampleSphere(sampler->next(), sampler->next());
		p = p * use<SceneBounds>().sceneRadius;
		p = p + use<SceneBounds>().sceneCentre;
		pdf = 1.0f / (4 * M_PI * SQ(use<SceneBounds>().sceneRadius));
		return p;
	}
	Vec3 sampleDirectionFromLight(Sampler* sampler, float& pdf) override {
		int x, y;
		float u = sampler->next();
		float v = sampler->next();
		distribution->sample(u, v, x, y, pdf);

		float uCoord = (x + 0.5f) / float(env->width);
		float vCoord = (y + 0.5f) / float(env->height);
		float phi = uCoord * 2.0f * M_PI;
		float theta = vCoord * M_PI;

		Vec3 wi = Vec3(
			sin(theta) * cos(phi),
			cos(theta),
			sin(theta) * sin(phi)
		);

		// Convert pdf to solid angle domain
		pdf = pdf / (2.0f * M_PI * M_PI * sinf(theta));
		return wi;
	}

};