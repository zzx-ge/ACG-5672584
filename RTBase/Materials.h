#pragma once

#include "Core.h"
#include "Imaging.h"
#include "Sampling.h"
#include "cmath"

#pragma warning( disable : 4244)

struct RefractionInfo {
	float etaI, etaT, eta;      // 折射率
	float cosThetaO;            // 出射方向与法线夹角的余弦值
	bool entering;              // 是否从外向内
	Vec3 normal;                // 修正后的法线方向
};

float Clamp(float x, float a, float b) {
	return std::max(a, std::min(x, b));
}

Vec3 Reflect(const Vec3& v, const Vec3& n) {
	return v - n * 2.0f * Dot(v, n);
}

bool Refract(const Vec3& wo, Vec3& wi, const RefractionInfo& info) {
	float sin2ThetaI = std::max(0.0f, 1.0f - info.cosThetaO * info.cosThetaO);
	float sin2ThetaT = info.eta * info.eta * sin2ThetaI;

	if (sin2ThetaT >= 1.0f)
		return false; // 全反射

	float cosThetaT = std::sqrt(1.0f - sin2ThetaT);

	// 折射方向计算
	wi = -wo * info.eta + info.normal * (info.eta * info.cosThetaO - cosThetaT);
	return true;
}



class BSDF;

class ShadingData
{
public:
	Vec3 x;
	Vec3 wo;
	Vec3 sNormal;
	Vec3 gNormal;
	float tu;
	float tv;
	Frame frame;
	BSDF* bsdf;
	float t;
	ShadingData() {}
	ShadingData(Vec3 _x, Vec3 n)
	{
		x = _x;
		gNormal = n;
		sNormal = n;
		bsdf = NULL;
	}
};

class ShadingHelper
{
public:
	// GGX Normal Distribution Function
	static float Dggx(Vec3 h, float alpha)
	{
		// Add code here
		float a2 = alpha * alpha;
		float NdotH = std::max(h.z, 0.f);
		float denom = (NdotH * NdotH) * (a2 - 1.0f) + 1.0f;
		return a2 / (M_PI * denom * denom);
		//return 1.0f;
	}

	// GGX Lambda term (Used in Smith G)
	static float lambdaGGX(Vec3 wi, float alpha)
	{
		// Add code here
		float absCosTheta = std::max(wi.z, 1e-5f);
		float tanTheta = sqrtf(1.0f - absCosTheta * absCosTheta) / absCosTheta;
		float a = 1.0f / (alpha * tanTheta);
		if (a < 1.6f)
		{
			return (1.0f - 1.259f * a + 0.396f * a * a) / (3.535f * a + 2.181f * a * a);
		}
		else
		{
			return 0.0f;
		}
		//return 1.0f;
	}

	// Smith 's G term
	static float Gggx(Vec3 wi, Vec3 wo, float alpha)
	{
		// Add code here
		return 1.0f / (1.0f + lambdaGGX(wi, alpha) + lambdaGGX(wo, alpha));
		//return 1.0f;
	}

	// Fresnel term for dielectric materials using schlick approximation
	static float fresnelDielectric(float cosTheta, float iorInt, float iorExt)
	{
		// If exiting the material
		// If cosTheta is negative, the ray is inside the material
		bool exiting = cosTheta > 0.0f;
		if (exiting) {
			std::swap(iorInt, iorExt);
			cosTheta = fabsf(cosTheta);
		}
		// Calculate refracted direction using snell's law
		float eta = iorInt / iorExt;
		float sin2ThetaT = eta * eta * (1.0f - cosTheta * cosTheta);
		// Total internal reflection
		if (sin2ThetaT > 1.0f) return 1.0f;
		float cosThetaT = sqrtf(1.0f - sin2ThetaT);
		// Calculate Fresnel (|| and T)
		float rParallel = (iorExt * cosTheta - iorInt * cosThetaT) / (iorExt * cosTheta + iorInt * cosThetaT);
		float rPerpendicular = (iorInt * cosTheta - iorExt * cosThetaT) / (iorInt * cosTheta + iorExt * cosThetaT);
		// return average^2
		float reflectance = 0.5f * (rParallel * rParallel + rPerpendicular * rPerpendicular);
		return reflectance;
		//return 1.0f;
	}
	static Colour fresnelCondutor(float cosTheta, Colour ior, Colour k)
	{
		// Clamp to valid cosine range
		cosTheta = clamp(cosTheta, 0.0f, 1.0f);
		float cos2 = cosTheta * cosTheta;
		Colour cos2C(cos2);
		Colour sin2C = Colour(1.0f) - cos2C;

		// Compute η² and k²
		Colour ior2 = ior * ior;
		Colour k2 = k * k;

		// Compute intermediate terms
		Colour t0 = ior2 - k2 - sin2C;
		Colour a2b2 = Colour::sqrt(t0 * t0 + ior2 * k2 * 4.0f);

		Colour t1 = a2b2 + cos2C;
		Colour Rs = (a2b2 - ior * cosTheta * 2.0f + cos2C) /
			(t1 + ior * cosTheta * 2.0f);

		Colour t2 = cos2C * a2b2 + sin2C * sin2C;
		Colour Rp = Rs * (t2 - ior * cosTheta * 2.0f + cos2C * cos2C) /
			(t2 + ior * cosTheta * 2.0f + cos2C * cos2C);

		// Return the average reflectance
		return (Rs + Rp) * 0.5;
	}
	static Colour fresnelSchlickConductor(float cosTheta, Colour ior, Colour k)
	{
		cosTheta = clamp(cosTheta, 0.0f, 1.0f);

		// Compute F0 = ((ior - 1)^2 + k^2) / ((ior + 1)^2 + k^2)
		Colour one(1.0f);
		Colour iorMinusOne = ior - one;
		Colour iorPlusOne = ior + one;

		Colour iorMinusOne2 = iorMinusOne * iorMinusOne;
		Colour iorPlusOne2 = iorPlusOne * iorPlusOne;

		Colour numerator = iorMinusOne2 + k * k;
		Colour denominator = iorPlusOne2 + k * k;

		Colour F0 = numerator / denominator;

		// Schlick's approximation
		float oneMinusCos = 1.0f - cosTheta;
		float factor = std::pow(oneMinusCos, 5.0f);
		return F0 + (Colour(1.0f) - F0) * factor;
	}
};

class BSDF
{
public:
	Colour emission;
	virtual Vec3 sample(const ShadingData& shadingData, Sampler* sampler, Colour& reflectedColour, float& pdf) = 0;
	virtual Colour evaluate(const ShadingData& shadingData, const Vec3& wi) = 0;
	virtual float PDF(const ShadingData& shadingData, const Vec3& wi) = 0;
	virtual bool isPureSpecular() = 0;
	virtual bool isTwoSided() = 0;
	bool isLight()
	{
		return emission.Lum() > 0 ? true : false;
	}
	void addLight(Colour _emission)
	{
		emission = _emission;
	}
	Colour emit(const ShadingData& shadingData, const Vec3& wi)
	{
		return emission;
	}
	Colour emit()
	{
		return emission;
	}
	virtual float mask(const ShadingData& shadingData) = 0;
};


class DiffuseBSDF : public BSDF
{
public:
	Texture* albedo;
	DiffuseBSDF() = default;
	DiffuseBSDF(Texture* _albedo)
	{
		albedo = _albedo;
	}
	Vec3 sample(const ShadingData& shadingData, Sampler* sampler, Colour& reflectedColour, float& pdf) override
	{// Cosine-weighted sampling over hemisphere
		Vec3 localWi = SamplingDistributions::cosineSampleHemisphere(sampler->next(), sampler->next());
		pdf = localWi.z / M_PI;

		Vec3 wi = shadingData.frame.toWorld(localWi);

		// Evaluate BRDF × cosθ
		Colour albedoVal = albedo->sample(shadingData.tu, shadingData.tv);
		reflectedColour = (albedoVal / M_PI);
		return wi;
	}
	Colour evaluate(const ShadingData& shadingData, const Vec3& wi) override
	{
		// Add correct evaluation code here
		Vec3 localWi = shadingData.frame.toLocal(wi);
		if (localWi.z <= 0.0f) return Colour(0.0f);
		return albedo->sample(shadingData.tu, shadingData.tv) / M_PI;
	}
	float PDF(const ShadingData& shadingData, const Vec3& wi) override
	{
		// Add correct PDF code here
		Vec3 localWi = shadingData.frame.toLocal(wi);
		if (localWi.z <= 0.0f) return 0.0f;
		return SamplingDistributions::cosineHemispherePDF(localWi);
	}
	bool isPureSpecular() override
	{
		return false;
	}
	bool isTwoSided() override
	{
		return true;
	}
	float mask(const ShadingData& shadingData) override
	{
		return albedo->sampleAlpha(shadingData.tu, shadingData.tv);
	}
};

class MirrorBSDF : public BSDF
{
public:
	Texture* albedo;
	MirrorBSDF() = default;
	MirrorBSDF(Texture* _albedo)
	{
		albedo = _albedo;
	}
	Vec3 sample(const ShadingData& shadingData, Sampler* sampler, Colour& reflectedColour, float& pdf) override
	{
		Vec3 woLocal = shadingData.frame.toLocal(shadingData.wo);

		Vec3 wiLocal(-woLocal.x, -woLocal.y, woLocal.z);

		Vec3 wi = shadingData.frame.toWorld(wiLocal);

		Colour albedoVal = albedo->sample(shadingData.tu, shadingData.tv);

		float cosTheta = std::max(0.0f, Dot(wi, shadingData.sNormal));
		if (cosTheta > 0.0f)
		{
			reflectedColour = albedoVal;
		}
		else
		{
			reflectedColour = Colour(0.0f);
		}

		pdf = 1.0f;

		return wi;
	}
	Colour evaluate(const ShadingData& shadingData, const Vec3& wi) override
	{
		// Add correct evaluation code here
		return Colour(0.f);
	}
	float PDF(const ShadingData& shadingData, const Vec3& wi) override
	{
		// Add correct PDF code here
		return 0.f;
	}
	bool isPureSpecular() override
	{
		return true;
	}
	bool isTwoSided() override
	{
		return true;
	}
	float mask(const ShadingData& shadingData) override
	{
		return albedo->sampleAlpha(shadingData.tu, shadingData.tv);
	}
};


class ConductorBSDF : public BSDF
{
public:
	Texture* albedo;
	Colour eta;
	Colour k;
	float alpha;
	ConductorBSDF() = default;
	ConductorBSDF(Texture* _albedo, Colour _eta, Colour _k, float roughness)
	{
		albedo = _albedo;
		eta = _eta;
		k = _k;
		alpha = 1.62142f * sqrtf(roughness);
	}
	Vec3 sample(const ShadingData& shadingData, Sampler* sampler, Colour& reflectedColour, float& pdf) override
	{
		// Replace this with Conductor sampling code
		Vec3 wi = SamplingDistributions::cosineSampleHemisphere(sampler->next(), sampler->next());
		pdf = wi.z / M_PI;
		reflectedColour = albedo->sample(shadingData.tu, shadingData.tv) / M_PI;
		wi = shadingData.frame.toWorld(wi);
		return wi;
	}
	Colour evaluate(const ShadingData& shadingData, const Vec3& wi) override
	{
		// Replace this with Conductor evaluation code
		return albedo->sample(shadingData.tu, shadingData.tv) / M_PI;
	}
	float PDF(const ShadingData& shadingData, const Vec3& wi) override
	{
		// Replace this with Conductor PDF
		Vec3 wiLocal = shadingData.frame.toLocal(wi);
		return SamplingDistributions::cosineHemispherePDF(wiLocal);
	}
	bool isPureSpecular() override
	{
		return false;
	}
	bool isTwoSided() override
	{
		return true;
	}
	float mask(const ShadingData& shadingData) override
	{
		return albedo->sampleAlpha(shadingData.tu, shadingData.tv);
	}
};

class GlassBSDF : public BSDF
{
public:
	Texture* albedo;
	float intIOR;
	float extIOR;
	GlassBSDF() = default;
	GlassBSDF(Texture* _albedo, float _intIOR, float _extIOR)
	{
		albedo = _albedo;
		intIOR = _intIOR;
		extIOR = _extIOR;
	}
	Vec3 sample(const ShadingData& shadingData, Sampler* sampler, Colour& reflectedColour, float& pdf) override
	{
		Vec3 wo = shadingData.wo;
		Vec3 n = shadingData.sNormal;

		// 提取折射信息
		RefractionInfo info = prepareRefraction(wo, n, extIOR, intIOR);
		float F = ShadingHelper::fresnelDielectric(info.cosThetaO, info.etaI, info.etaT);

		Colour baseColour = albedo->sample(shadingData.tu, shadingData.tv);

		if (sampler->next() < F) {
			// 反射分支
			Vec3 wi = Reflect(-wo, info.normal);
			reflectedColour = baseColour * F;
			pdf = F;
			return wi;
		}
		else {
			// 折射分支
			Vec3 wi;
			if (!Refract(wo, wi, info)) {
				// 发生全反射 fallback
				wi = Reflect(-wo, info.normal);
				reflectedColour = baseColour;
				pdf = 1.0f;
				return wi;
			}

			reflectedColour = baseColour * (1.0f - F);
			pdf = 1.0f - F;
			return wi;
		}
	}

	Colour evaluate(const ShadingData& shadingData, const Vec3& wi) override
	{
		// Add correct evaluation code here
		return 0.f;
	}
	float PDF(const ShadingData& shadingData, const Vec3& wi) override
	{
		// Add correct PDF code here
		return 0.f;
	}
	bool isPureSpecular() override
	{
		return true;
	}
	bool isTwoSided() override
	{
		return false;
	}
	float mask(const ShadingData& shadingData) override
	{
		return albedo->sampleAlpha(shadingData.tu, shadingData.tv);
	}

	RefractionInfo prepareRefraction(const Vec3& wo, const Vec3& n, float extIOR, float intIOR) {
		RefractionInfo info;
		info.entering = Dot(wo, n) > 0.0f;
		info.normal = info.entering ? n : -n;
		info.etaI = info.entering ? extIOR : intIOR;
		info.etaT = info.entering ? intIOR : extIOR;
		info.eta = info.etaI / info.etaT;
		info.cosThetaO = Clamp(Dot(wo, info.normal), -1.0f, 1.0f);
		return info;
	}
};

class DielectricBSDF : public BSDF
{
public:
	Texture* albedo;
	float intIOR;
	float extIOR;
	float alpha;
	DielectricBSDF() = default;
	DielectricBSDF(Texture* _albedo, float _intIOR, float _extIOR, float roughness)
	{
		albedo = _albedo;
		intIOR = _intIOR;
		extIOR = _extIOR;
		alpha = 1.62142f * sqrtf(roughness);
	}
	Vec3 sample(const ShadingData& shadingData, Sampler* sampler, Colour& reflectedColour, float& pdf) override
	{
		// Replace this with Dielectric sampling code
		Vec3 wi = SamplingDistributions::cosineSampleHemisphere(sampler->next(), sampler->next());
		pdf = wi.z / M_PI;
		reflectedColour = albedo->sample(shadingData.tu, shadingData.tv) / M_PI;
		wi = shadingData.frame.toWorld(wi);
		return wi;
	}
	Colour evaluate(const ShadingData& shadingData, const Vec3& wi) override
	{
		// Replace this with Dielectric evaluation code
		return albedo->sample(shadingData.tu, shadingData.tv) / M_PI;
	}
	float PDF(const ShadingData& shadingData, const Vec3& wi) override
	{
		// Replace this with Dielectric PDF
		Vec3 wiLocal = shadingData.frame.toLocal(wi);
		return SamplingDistributions::cosineHemispherePDF(wiLocal);
	}
	bool isPureSpecular() override
	{
		return false;
	}
	bool isTwoSided() override
	{
		return false;
	}
	float mask(const ShadingData& shadingData) override
	{
		return albedo->sampleAlpha(shadingData.tu, shadingData.tv);
	}
};

class OrenNayarBSDF : public BSDF
{
public:
	Texture* albedo;
	float sigma;
	OrenNayarBSDF() = default;
	OrenNayarBSDF(Texture* _albedo, float _sigma)
	{
		albedo = _albedo;
		sigma = _sigma;
	}
	Vec3 sample(const ShadingData& shadingData, Sampler* sampler, Colour& reflectedColour, float& pdf) override
	{
		// Replace this with OrenNayar sampling code
		Vec3 wi = SamplingDistributions::cosineSampleHemisphere(sampler->next(), sampler->next());
		pdf = wi.z / M_PI;
		reflectedColour = albedo->sample(shadingData.tu, shadingData.tv) / M_PI;
		wi = shadingData.frame.toWorld(wi);
		return wi;
	}
	Colour evaluate(const ShadingData& shadingData, const Vec3& wi) override
	{
		// Replace this with OrenNayar evaluation code
		return albedo->sample(shadingData.tu, shadingData.tv) / M_PI;
	}
	float PDF(const ShadingData& shadingData, const Vec3& wi) override
	{
		// Replace this with OrenNayar PDF
		Vec3 wiLocal = shadingData.frame.toLocal(wi);
		return SamplingDistributions::cosineHemispherePDF(wiLocal);
	}
	bool isPureSpecular() override
	{
		return false;
	}
	bool isTwoSided() override
	{
		return true;
	}
	float mask(const ShadingData& shadingData) override
	{
		return albedo->sampleAlpha(shadingData.tu, shadingData.tv);
	}
};

class PlasticBSDF : public BSDF
{
public:
	Texture* albedo;
	float intIOR;
	float extIOR;
	float alpha;
	PlasticBSDF() = default;
	PlasticBSDF(Texture* _albedo, float _intIOR, float _extIOR, float roughness)
	{
		albedo = _albedo;
		intIOR = _intIOR;
		extIOR = _extIOR;
		alpha = 1.62142f * sqrtf(roughness);
	}
	float alphaToPhongExponent()
	{
		return (2.0f / SQ(std::max(alpha, 0.001f))) - 2.0f;
	}
	Vec3 sample(const ShadingData& shadingData, Sampler* sampler, Colour& reflectedColour, float& pdf) override
	{
		// Replace this with Plastic sampling code
		Vec3 wi = SamplingDistributions::cosineSampleHemisphere(sampler->next(), sampler->next());
		pdf = wi.z / M_PI;
		reflectedColour = albedo->sample(shadingData.tu, shadingData.tv) / M_PI;
		wi = shadingData.frame.toWorld(wi);
		return wi;
	}
	Colour evaluate(const ShadingData& shadingData, const Vec3& wi) override
	{
		// Replace this with Plastic evaluation code
		return albedo->sample(shadingData.tu, shadingData.tv) / M_PI;
	}
	float PDF(const ShadingData& shadingData, const Vec3& wi) override
	{
		// Replace this with Plastic PDF
		Vec3 wiLocal = shadingData.frame.toLocal(wi);
		return SamplingDistributions::cosineHemispherePDF(wiLocal);
	}
	bool isPureSpecular() override
	{
		return false;
	}
	bool isTwoSided() override
	{
		return true;
	}
	float mask(const ShadingData& shadingData) override
	{
		return albedo->sampleAlpha(shadingData.tu, shadingData.tv);
	}
};

class LayeredBSDF : public BSDF
{
public:
	BSDF* base;
	Colour sigmaa;
	float thickness;
	float intIOR;
	float extIOR;
	LayeredBSDF() = default;
	LayeredBSDF(BSDF* _base, Colour _sigmaa, float _thickness, float _intIOR, float _extIOR)
	{
		base = _base;
		sigmaa = _sigmaa;
		thickness = _thickness;
		intIOR = _intIOR;
		extIOR = _extIOR;
	}
	Vec3 sample(const ShadingData& shadingData, Sampler* sampler, Colour& reflectedColour, float& pdf) override
	{
		// Add code to include layered sampling
		return base->sample(shadingData, sampler, reflectedColour, pdf);
	}
	Colour evaluate(const ShadingData& shadingData, const Vec3& wi) override
	{
		// Add code for evaluation of layer
		return base->evaluate(shadingData, wi);
	}
	float PDF(const ShadingData& shadingData, const Vec3& wi) override
	{
		// Add code to include PDF for sampling layered BSDF
		return base->PDF(shadingData, wi);
	}
	bool isPureSpecular() override
	{
		return base->isPureSpecular();
	}
	bool isTwoSided() override
	{
		return true;
	}
	float mask(const ShadingData& shadingData) override
	{
		return base->mask(shadingData);
	}
};