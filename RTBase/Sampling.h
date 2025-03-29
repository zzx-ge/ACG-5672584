#pragma once

#include "Core.h"
#include <random>
#include <algorithm>

class Sampler
{
public:
	virtual float next() = 0;
};

class MTRandom : public Sampler
{
public:
	std::mt19937 generator;
	std::uniform_real_distribution<float> dist;
	explicit MTRandom(unsigned int seed = 10) : dist(0.0f, 1.0f)
	{
		generator.seed(seed);
	}
	float next() override
	{
		return dist(generator);
	}
};

// Note all of these distributions assume z-up coordinate system
class SamplingDistributions
{
public:
	static Vec3 uniformSampleHemisphere(float r1, float r2)
	{
		// Add code here
		return Vec3(0, 0, 1);
	}
	static float uniformHemispherePDF(const Vec3 wi)
	{
		// Add code here
		return 1.0f;
	}

	static Vec3 cosineSampleHemisphere(float r1, float r2)
	{
		float z = sqrt(1.0f - r1);  // Corrected formula
		float r = sqrt(r1);         // Compute circle radius at z
		float phi = 2.0f * M_PI * r2;

		float x = r * cos(phi);
		float y = r * sin(phi);

		return Vec3(x, y, z);       // Return uniformly sampled point on hemisphere
	}
	static float cosineHemispherePDF(const Vec3 wi)
	{
		// Add code here
		return wi.z / M_PI;
	}
	static Vec3 uniformSampleSphere(float r1, float r2)
	{
		float z = 1.0f - 2.0f * r1;       // Map r1 to [-1,1] for uniform z
		float r = sqrt(1.0f - z * z);     // Compute circle radius at z
		float phi = 2.0f * M_PI * r2;     // Uniform azimuthal angle

		float x = r * cos(phi);
		float y = r * sin(phi);

		return Vec3(x, y, z);             // Return uniformly sampled point on sphere
	}
	static float uniformSpherePDF(const Vec3& wi)
	{
		// Add code here
		return 1.0f/(4.0f * M_PI);
	}
};

class Distribution1D {
public:
	std::vector<float> func;      // original function
	std::vector<float> cdf;       // cdf
	float integral = 0.0f;
	Distribution1D() {}

	Distribution1D(const std::vector<float>& f) {
		int n = f.size();
		func = f;
		cdf.resize(n + 1);
		cdf[0] = 0.0f;

		// buil CDF
		for (int i = 0; i < n; i++) {
			cdf[i + 1] = cdf[i] + func[i];
		}

		integral = cdf[n];
		for (int i = 0; i <= n; i++) {
			cdf[i] /= integral; // Normalize
		}
	}

	// sample from cdf
	int sample(float u, float& pdf) const {
		auto it = std::upper_bound(cdf.begin(), cdf.end(), u);
		int idx = std::max(0, int(it - cdf.begin() - 1));

		float du = (u - cdf[idx]) / (cdf[idx + 1] - cdf[idx]);
		pdf = func[idx] / integral;
		return idx;
	}
};

class Distribution2D {
public:
	std::vector<Distribution1D> conditional; // line sample
	Distribution1D marginal;                // line setection

	Distribution2D(const std::vector<std::vector<float>>& image) {
		int height = image.size();
		int width = image[0].size();
		std::vector<float> marginalWeights;

		for (int y = 0; y < height; ++y) {
			Distribution1D row(image[y]);
			conditional.push_back(row);
			marginalWeights.push_back(row.integral);
		}

		marginal = Distribution1D(marginalWeights);
	}

	// sample from 2D distribution
	void sample(float u1, float u2, int& x, int& y, float& pdf) {
		float pdfY, pdfX;

		y = marginal.sample(u1, pdfY);
		x = conditional[y].sample(u2, pdfX);

		pdf = pdfX * pdfY;
	}
};

