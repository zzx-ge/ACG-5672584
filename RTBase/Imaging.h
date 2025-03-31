#pragma once

#include "Core.h"
#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"
#define STB_IMAGE_WRITE_IMPLEMENTATION
#define __STDC_LIB_EXT1__
#include "stb_image_write.h"


// Stop warnings about buffer overruns if size is zero. Size should never be zero and if it is the code handles it.
#pragma warning( disable : 6386)

template <typename T>
T clamp(T value, T minVal, T maxVal) {
	return (value < minVal) ? minVal : (value > maxVal) ? maxVal : value;
}


class Texture
{
public:
	Colour* texels;
	float* alpha;
	int width;
	int height;
	int channels;
	void loadDefault()
	{
		width = 1;
		height = 1;
		channels = 3;
		texels = new Colour[1];
		texels[0] = Colour(1.0f, 1.0f, 1.0f);
	}
	void load(std::string filename)
	{
		alpha = NULL;
		if (filename.find(".hdr") != std::string::npos)
		{
			float* textureData = stbi_loadf(filename.c_str(), &width, &height, &channels, 0);
			if (width == 0 || height == 0)
			{
				loadDefault();
				return;
			}
			texels = new Colour[width * height];
			for (int i = 0; i < (width * height); i++)
			{
				texels[i] = Colour(textureData[i * channels], textureData[(i * channels) + 1], textureData[(i * channels) + 2]);
			}
			stbi_image_free(textureData);
			return;
		}
		unsigned char* textureData = stbi_load(filename.c_str(), &width, &height, &channels, 0);
		if (width == 0 || height == 0)
		{
			loadDefault();
			return;
		}
		texels = new Colour[width * height];
		for (int i = 0; i < (width * height); i++)
		{
			texels[i] = Colour(textureData[i * channels] / 255.0f, textureData[(i * channels) + 1] / 255.0f, textureData[(i * channels) + 2] / 255.0f);
		}
		if (channels == 4)
		{
			alpha = new float[width * height];
			for (int i = 0; i < (width * height); i++)
			{
				alpha[i] = textureData[(i * channels) + 3] / 255.0f;
			}
		}
		stbi_image_free(textureData);
	}
	Colour sample(const float tu, const float tv) const
	{
		Colour tex;
		float u = std::max(0.0f, fabsf(tu)) * width;
		float v = std::max(0.0f, fabsf(tv)) * height;
		int x = (int)floorf(u);
		int y = (int)floorf(v);
		float frac_u = u - x;
		float frac_v = v - y;
		float w0 = (1.0f - frac_u) * (1.0f - frac_v);
		float w1 = frac_u * (1.0f - frac_v);
		float w2 = (1.0f - frac_u) * frac_v;
		float w3 = frac_u * frac_v;
		x = x % width;
		y = y % height;
		Colour s[4];
		s[0] = texels[y * width + x];
		s[1] = texels[y * width + ((x + 1) % width)];
		s[2] = texels[((y + 1) % height) * width + x];
		s[3] = texels[((y + 1) % height) * width + ((x + 1) % width)];
		tex = (s[0] * w0) + (s[1] * w1) + (s[2] * w2) + (s[3] * w3);
		return tex;
	}
	float sampleAlpha(const float tu, const float tv) const
	{
		if (alpha == NULL)
		{
			return 1.0f;
		}
		float tex;
		float u = std::max(0.0f, fabsf(tu)) * width;
		float v = std::max(0.0f, fabsf(tv)) * height;
		int x = (int)floorf(u);
		int y = (int)floorf(v);
		float frac_u = u - x;
		float frac_v = v - y;
		float w0 = (1.0f - frac_u) * (1.0f - frac_v);
		float w1 = frac_u * (1.0f - frac_v);
		float w2 = (1.0f - frac_u) * frac_v;
		float w3 = frac_u * frac_v;
		x = x % width;
		y = y % height;
		float s[4];
		s[0] = alpha[y * width + x];
		s[1] = alpha[y * width + ((x + 1) % width)];
		s[2] = alpha[((y + 1) % height) * width + x];
		s[3] = alpha[((y + 1) % height) * width + ((x + 1) % width)];
		tex = (s[0] * w0) + (s[1] * w1) + (s[2] * w2) + (s[3] * w3);
		return tex;
	}
	~Texture()
	{
		delete[] texels;
		if (alpha != NULL)
		{
			delete alpha;
		}
	}
};

class ImageFilter
{
public:
	virtual float filter(const float x, const float y) const = 0;
	virtual int size() const = 0;
};

class BoxFilter : public ImageFilter
{
public:
	float filter(float x, float y) const
	{
		if (fabsf(x) <= 0.5f && fabs(y) <= 0.5f)
		{
			return 1.0f;
		}
		return 0;
	}
	int size() const
	{
		return 0;
	}
};

class Gaussian_filter : public ImageFilter
{
private:
	float sigma;
	int filterSize;
public:
	Gaussian_filter(float sigma) : sigma(sigma), filterSize(static_cast<int>(std::ceil(3 * sigma))) {}

	// Gaussian function implementation
	float filter(const float x, const float y) const override {
		float factor = 1.0f / (2.0f * M_PI * sigma * sigma);
		float exponent = -(x * x + y * y) / (2.0f * sigma * sigma);
		return factor * std::exp(exponent);
	}
	// Return the size of the filter (usually 3 * sigma for reasonable accuracy)
	int size() const override {
		return filterSize;
	}
};

class MitchellNetravaliFilter : public ImageFilter {
private:
	float B, C;  // Mitchell-Netravali parameters
	int filterSize;  // Precomputed filter size

	// Mitchell-Netravali function
	float mitchellNetravali1D(float x) const {
		x = std::abs(x);
		if (x < 1) {
			return ((1.0f / 6.0f) * ((12 - 9 * B - 6 * C) * x * x * x +
				(-18 + 12 * B + 6 * C) * x * x +
				(6 - 2 * B)));
		}
		else if (x < 2) {
			return ((1.0f / 6.0f) * ((-B - 6 * C) * x * x * x +
				(6 * B + 30 * C) * x * x +
				(-12 * B - 48 * C) * x +
				(8 * B + 24 * C)));
		}
		return 0.0f;
	}

public:
	// Constructor with default Mitchell-Netravali parameters (B=1/3, C=1/3)
	explicit MitchellNetravaliFilter(float B = 1.0f / 3.0f, float C = 1.0f / 3.0f)
		: B(B), C(C), filterSize(2) {
	} // Filter is typically 2 pixels wide

// Compute the 2D filter response as a separable function
	float filter(const float x, const float y) const override {
		return mitchellNetravali1D(x) * mitchellNetravali1D(y);
	}

	// Return the filter size (usually 2 for Mitchell-Netravali)
	int size() const override {
		return filterSize;
	}
};

class Film
{
public:
	Colour* film;
	unsigned int width;
	unsigned int height;
	int SPP;
	ImageFilter* filter;
	void splat(const float x, const float y, const Colour& L) {
		float filterWeights[25]; // Storage to cache weightsunsigned int indices[25]; // Store indices to minimize computations unsigned int used = 0;
		unsigned int indices[25];
		unsigned int used = 0;
		float total = 0;
		int size = filter->size();
		for (int i = -size; i <= size; i++) {
			for (int j = -size; j <= size; j++) {
				int px = (int)x + j;
				int py = (int)y + i;
				if (px >= 0 && px < width && py >= 0 && py < height) {
					indices[used] = (py * width) + px;
					filterWeights[used] = filter->filter(j, i); total += filterWeights[used];
					used++;
				}
			}
		}
		for (int i = 0; i < used; i++) {
			film[indices[i]] = film[indices[i]] + (L * filterWeights[i] / total);
		}
	}
	void tonemap(float r, float g, float b, unsigned char& outR, unsigned char& outG, unsigned char& outB, float exposure = 1.f)
	{
		// Apply exposure correction
		r *= exposure;
		g *= exposure;
		b *= exposure;

		const float epsilon = 1e-6f;

		// Luminance-based Reinhard Tonemapping
		float L = 0.2126f * r + 0.7152f * g + 0.0722f * b + epsilon; // Compute luminance
		float scale = L / (1.0f + L); // Reinhard Tonemapping: L_out = L / (1 + L)

		// Reinhard Tonemapping: L_out = L / (1 + L)
		if (L > 0.0f)
		{
			r *= scale / L;
			g *= scale / L;
			b *= scale / L;
		}

		// Gamma correction (linear to sRGB)
		r = pow(clamp(r, 0.0f, 1.0f), 1.0f / 2.2f);
		g = pow(clamp(g, 0.0f, 1.0f), 1.0f / 2.2f);
		b = pow(clamp(b, 0.0f, 1.0f), 1.0f / 2.2f);

		// Convert to 8-bit (0-255)
		outR = static_cast<unsigned char>(clamp(r * 255.0f, 0.0f, 255.0f));
		outG = static_cast<unsigned char>(clamp(g * 255.0f, 0.0f, 255.0f));
		outB = static_cast<unsigned char>(clamp(b * 255.0f, 0.0f, 255.0f));
	}
	// Do not change any code below this line
	void init(int _width, int _height, ImageFilter* _filter)
	{
		width = _width;
		height = _height;
		film = new Colour[width * height];
		clear();
		filter = _filter;
	}
	void clear()
	{
		memset(film, 0, width * height * sizeof(Colour));
		SPP = 0;
	}
	void incrementSPP()
	{
		SPP++;
	}
	void save(std::string filename)
	{
		Colour* hdrpixels = new Colour[width * height];
		for (unsigned int i = 0; i < (width * height); i++)
		{
			hdrpixels[i] = film[i] / (float)SPP;
		}
		stbi_write_hdr(filename.c_str(), width, height, 3, (float*)hdrpixels);
		delete[] hdrpixels;
	}
};