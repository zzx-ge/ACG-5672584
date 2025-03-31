#pragma once

#include "Core.h"
#include "Sampling.h"

struct IntersectionData
{
	unsigned int ID;
	float t;
	float alpha;
	float beta;
	float gamma;
};

struct PathVertex {
	Vec3 position;
	Vec3 normal;
	Vec3 wi;
	Colour beta;
	float pdf;
	const BSDF* bsdf;
};

class Ray
{
public:
	Vec3 o;
	Vec3 dir;
	Vec3 invDir;
	Ray()
	{
	}
	Ray(Vec3 _o, Vec3 _d)
	{
		init(_o, _d);
	}
	void init(Vec3 _o, Vec3 _d)
	{
		o = _o;
		dir = _d;
		invDir = Vec3(1.0f / dir.x, 1.0f / dir.y, 1.0f / dir.z);
	}
	Vec3 at(const float t) const
	{
		return (o + (dir * t));
	}
};

class Plane
{
public:
	Vec3 n;
	float d;
	void init(Vec3& _n, float _d)
	{
		n = _n;
		d = _d;
	}
	// Add code here
	bool rayIntersect(Ray& r, float& t)
	{
		return false;
	}
};

#define EPSILON 0.001f

class Triangle
{
public:
	Vertex vertices[3];
	Vec3 e1; // Edge 1
	Vec3 e2; // Edge 2
	Vec3 n; // Geometric Normal
	float area; // Triangle area
	float d; // For ray triangle if needed
	unsigned int materialIndex;
	Vec3 samplepoint = Vec3(-1.f, -1.f, -1.f);
	void init(Vertex v0, Vertex v1, Vertex v2, unsigned int _materialIndex)
	{
		materialIndex = _materialIndex;
		vertices[0] = v0;
		vertices[1] = v1;
		vertices[2] = v2;
		e1 = vertices[1].p - vertices[0].p;
		e2 = vertices[2].p - vertices[0].p;
		n = e1.cross(e2).normalize();
		area = e1.cross(e2).length() * 0.5f;
		d = Dot(n, vertices[0].p);
	}
	Vec3 centre() const
	{
		return (vertices[0].p + vertices[1].p + vertices[2].p) / 3.0f;
	}
	// Add code here
	bool rayIntersect(const Ray& r, float& t, float& u, float& v) const
	{
		const float TOLERANCE = 1e-6f; // Tolerance to avoid floating point precision issues

		// Step 1: Compute determinant
		Vec3 P = r.dir.cross(e2);
		float det = e1.dot(P);

		// If determinant is near zero, the ray is parallel to the triangle
		if (fabs(det) < TOLERANCE)
			return false;

		float invDet = 1.0f / det;

		// Step 3: Compute barycentric coordinate u
		Vec3 T = r.o - vertices[0].p;
		u = T.dot(P) * invDet;
		if (u < 0.0f || u > 1.0f)
			return false;

		// Step 4: Compute barycentric coordinate v
		Vec3 Q = T.cross(e1);
		v = r.dir.dot(Q) * invDet;
		if (v < 0.0f || (u + v) > 1.0f)
			return false;

		// Step 5: Compute intersection distance t
		t = e2.dot(Q) * invDet;

		// If t < 0, the intersection is behind the ray
		return t >= 0.0f;
	}
	// Just for the intersection check
	bool rayIntersect(const Ray& r) const{
		const float TOLERANCE = 1e-6f; // Tolerance to avoid floating point precision issues

		// Step 1: Compute determinant
		Vec3 P = r.dir.cross(e2);
		float det = e1.dot(P);

		// If determinant is near zero, the ray is parallel to the triangle
		if (fabs(det) < TOLERANCE)
			return false;

		float invDet = 1.0f / det;

		// Step 3: Compute barycentric coordinate u
		Vec3 T = r.o - vertices[0].p;
		float u = T.dot(P) * invDet;
		if (u < 0.0f || u > 1.0f)
			return false;

		// Step 4: Compute barycentric coordinate v
		Vec3 Q = T.cross(e1);
		float v = r.dir.dot(Q) * invDet;
		if (v < 0.0f || (u + v) > 1.0f)
			return false;

		// Step 5: Compute intersection distance t
		float t = e2.dot(Q) * invDet;

		// If t < 0, the intersection is behind the ray
		return t >= 0.0f;
	}

	bool rayIntersect(const Ray& r, IntersectionData& hit) const{
		const float TOLERANCE = 1e-6f; // Tolerance to avoid floating point precision issues

		// Step 1: Compute determinant
		Vec3 P = r.dir.cross(e2);
		float det = e1.dot(P);

		// If determinant is near zero, the ray is parallel to the triangle
		if (fabs(det) < TOLERANCE)
			return false;

		float invDet = 1.0f / det;

		// Step 3: Compute barycentric coordinate u
		Vec3 T = r.o - vertices[0].p;
		float u = T.dot(P) * invDet;
		if (u < 0.0f || u > 1.0f)
			return false;
		hit.alpha = u;

		// Step 4: Compute barycentric coordinate v
		Vec3 Q = T.cross(e1);
		float v = r.dir.dot(Q) * invDet;
		if (v < 0.0f || (u + v) > 1.0f)
			return false;
		hit.beta = v;
		hit.gamma = 1.0f - (u + v);

		// Step 5: Compute intersection distance t
		float t = e2.dot(Q) * invDet;
		hit.t = t;

		// If t < 0, the intersection is behind the ray
		return t >= 0.0f;
	}

	void interpolateAttributes(const float alpha, const float beta, const float gamma, Vec3& interpolatedNormal, float& interpolatedU, float& interpolatedV) const
	{
		interpolatedNormal = vertices[0].normal * alpha + vertices[1].normal * beta + vertices[2].normal * gamma;
		interpolatedNormal = interpolatedNormal.normalize();
		interpolatedU = vertices[0].u * alpha + vertices[1].u * beta + vertices[2].u * gamma;
		interpolatedV = vertices[0].v * alpha + vertices[1].v * beta + vertices[2].v * gamma;
	}
	// Add code here
	Vec3 sample(Sampler* sampler, float& pdf)
	{
		pdf = 1 / area;
		float r1 = sampler->next();
		float r2 = sampler->next();
		float sr1 = sqrtf(r1);
		float alpha = 1.0f - sr1;
		float beta = r2 * sr1;
		samplepoint = vertices[0].p * (1 - alpha - beta) + vertices[1].p * alpha + vertices[2].p * beta;
		return samplepoint;
		//return Vec3(vertices[0].p + vertices[1].p * alpha + (vertices[2].p - vertices[0].p) * beta);
	}

	Vec3 samplePoint() const
	{
		if (samplepoint.x != -1) return samplepoint;
		else return this->centre();
	}

	Vec3 gNormal()
	{
		return (n * (Dot(vertices[0].normal, n) > 0 ? 1.0f : -1.0f));
	}

	Vec3 centroid() const {
		return (vertices[0].p + vertices[1].p + vertices[2].p) / 3;
	}

	bool inside(const Vec3 point, float epsilon = 1e-4f) {
		// 1. 投影判断：点是否在三角形平面上（用法线投影）
		float planeDist = Dot(n, point - vertices[0].p);
		if (std::abs(planeDist) > epsilon)
			return false;

		// 2. 使用重心坐标判断是否在三角形内部
		Vec3 v0 = vertices[1].p - vertices[0].p;
		Vec3 v1 = vertices[2].p - vertices[0].p;
		Vec3 v2 = point - vertices[0].p;

		float d00 = Dot(v0, v0);
		float d01 = Dot(v0, v1);
		float d11 = Dot(v1, v1);
		float d20 = Dot(v2, v0);
		float d21 = Dot(v2, v1);

		float denom = d00 * d11 - d01 * d01;
		if (std::abs(denom) < 1e-8f)
			return false; // 退化三角形

		float v = (d11 * d20 - d01 * d21) / denom;
		float w = (d00 * d21 - d01 * d20) / denom;
		float u = 1.0f - v - w;

		return (u >= -epsilon && v >= -epsilon && w >= -epsilon);
	}
};

class AABB
{
public:
	Vec3 max;
	Vec3 min;
	AABB()
	{
		reset();
	}
	void reset()
	{
		max = Vec3(-FLT_MAX, -FLT_MAX, -FLT_MAX);
		min = Vec3(FLT_MAX, FLT_MAX, FLT_MAX);
	}
	void extend(const Vec3 p)
	{
		max = Max(max, p);
		min = Min(min, p);
	}
	void extend(const Triangle& tri) {
		extend(tri.vertices[0].p);
		extend(tri.vertices[1].p);
		extend(tri.vertices[2].p);
	}
	// Add code here
	bool rayAABB(const Ray& r, float& t)
	{
		float tMin = 0.0f;
		float tMax = FLT_MAX;

		for (int i = 0; i < 3; i++)
		{
			float invD = 1.0f / r.dir[i]; // Inverse direction
			float t0 = (min[i] - r.o[i]) * invD;
			float t1 = (max[i] - r.o[i]) * invD;

			if (invD < 0.0f)
				std::swap(t0, t1); // Ensure t0 is always min

			tMin = std::max(tMin, t0);
			tMax = std::min(tMax, t1);

			if (tMax < tMin)
				return false;
		}

		t = tMin; // Return the entry point
		return true;
	}
	// Add code here
	bool rayAABB(const Ray& r)
	{
		float invDirX = 1.0f / r.dir.x;
		float invDirY = 1.0f / r.dir.y;
		float invDirZ = 1.0f / r.dir.z;

		float txmin = (min.x - r.o.x) * invDirX;
		float txmax = (max.x - r.o.x) * invDirX;
		if (txmin > txmax) std::swap(txmin, txmax);

		float tymin = (min.y - r.o.y) * invDirY;
		float tymax = (max.y - r.o.y) * invDirY;
		if (tymin > tymax) std::swap(tymin, tymax);

		if ((txmin > tymax) || (tymin > txmax)) return false;

		float tzmin = (min.z - r.o.z) * invDirZ;
		float tzmax = (max.z - r.o.z) * invDirZ;
		if (tzmin > tzmax) std::swap(tzmin, tzmax);

		if ((txmin > tzmax) || (tzmin > txmax)) return false;

		return true;
	}
	// Add code here
	float area()
	{
		Vec3 size = max - min;
		return ((size.x * size.y) + (size.y * size.z) + (size.x * size.z)) * 2.0f;
	}
	static AABB Union(const AABB& a, const AABB& b) {
		AABB result;
		result.min = Vec3(
			std::min(a.min.x, b.min.x),
			std::min(a.min.y, b.min.y),
			std::min(a.min.z, b.min.z)
		);
		result.max = Vec3(
			std::max(a.max.x, b.max.x),
			std::max(a.max.y, b.max.y),
			std::max(a.max.z, b.max.z)
		);
		return result;
	}
};

class Sphere
{
public:
	Vec3 centre;
	float radius;
	void init(Vec3& _centre, float _radius)
	{
		centre = _centre;
		radius = _radius;
	}
	// Add code here
	bool rayIntersect(Ray& r, float& t)
	{
		return false;
	}
};


#define MAXNODE_TRIANGLES 8
#define TRAVERSE_COST 1.0f
#define TRIANGLE_COST 2.0f
#define BUILD_BINS 32

class BVHNode {
public:
	AABB bounds;
	BVHNode* l;
	BVHNode* r;
	std::vector<int> triangleIndices;  // Store indices instead of full triangle objects

	BVHNode() : l(nullptr), r(nullptr) {
		bounds.reset();
	}

	void build(std::vector<int>& inputIndices, std::vector<Triangle>& triangles, int depth = 0) {
		bounds.reset();
		for (int index : inputIndices) bounds.extend(triangles[index]);
		triangleIndices = inputIndices;

		if (inputIndices.size() <= 2 || depth > 32) return;

		const int numBuckets = 16;
		struct BucketInfo {
			int count = 0;
			AABB bounds;
		};

		int bestAxis = -1, bestSplit = -1;
		float bestCost = std::numeric_limits<float>::max();

		for (int axis = 0; axis < 3; ++axis) {
			// Compute bounds of triangle centroids
			AABB centroidBounds;
			for (int idx : inputIndices)
				centroidBounds.extend(triangles[idx].centroid());

			if (centroidBounds.max[axis] - centroidBounds.min[axis] < 1e-6f) continue;

			// Initialize buckets
			BucketInfo buckets[numBuckets];
			float scale = numBuckets / (centroidBounds.max[axis] - centroidBounds.min[axis]);

			for (int idx : inputIndices) {
				int b = std::min(numBuckets - 1, int((triangles[idx].centroid()[axis] - centroidBounds.min[axis]) * scale));
				buckets[b].count++;
				buckets[b].bounds.extend(triangles[idx]);
			}

			// Try all splits
			for (int i = 1; i < numBuckets; ++i) {
				AABB b0, b1;
				int count0 = 0, count1 = 0;
				for (int j = 0; j < i; ++j) {
					b0 = AABB::Union(b0, buckets[j].bounds);
					count0 += buckets[j].count;
				}
				for (int j = i; j < numBuckets; ++j) {
					b1 = AABB::Union(b1, buckets[j].bounds);
					count1 += buckets[j].count;
				}

				float cost = 0.125f + (count0 * b0.area() + count1 * b1.area()) / bounds.area();
				if (cost < bestCost) {
					bestCost = cost;
					bestAxis = axis;
					bestSplit = i;
				}
			}
		}

		// If no good split found, return as leaf
		if (bestAxis == -1) return;

		// Partition triangles
		AABB centroidBounds;
		for (int idx : inputIndices)
			centroidBounds.extend(triangles[idx].centroid());
		float scale = numBuckets / (centroidBounds.max[bestAxis] - centroidBounds.min[bestAxis]);

		auto mid = std::partition(inputIndices.begin(), inputIndices.end(), [&](int idx) {
			int b = std::min(numBuckets - 1, int((triangles[idx].centroid()[bestAxis] - centroidBounds.min[bestAxis]) * scale));
			return b < bestSplit;
			});

		std::vector<int> leftIndices(inputIndices.begin(), mid);
		std::vector<int> rightIndices(mid, inputIndices.end());

		if (leftIndices.empty() || rightIndices.empty()) return;

		l = new BVHNode();
		r = new BVHNode();
		l->build(leftIndices, triangles, depth + 1);
		r->build(rightIndices, triangles, depth + 1);
	}


	void traverse(const Ray& ray, const std::vector<Triangle>& triangles, IntersectionData& intersection) {
		// Check if ray intersects the bounding box
		if (!bounds.rayAABB(ray)) return;
		// If leaf node and intersect, check all triangles
		if (!l && !r) {
			for (int index : triangleIndices) {
				IntersectionData hit;
				if (triangles[index].rayIntersect(ray, hit) && hit.t < intersection.t) {
					intersection = hit;
					intersection.ID = index;
				}
			}
			return;
		}
		if (l) l->traverse(ray, triangles, intersection);
		if (r) r->traverse(ray, triangles, intersection);
	}

	bool traverseVisible(const Ray& ray, const std::vector<Triangle>& triangles, float maxT) {
		if (!bounds.rayAABB(ray)) return true;
		if (!l && !r) {
			for (int index : triangleIndices) {
				IntersectionData hit;
				if (triangles[index].rayIntersect(ray, hit) && hit.t < maxT) return false;
			}
			return true;
		}
		if (l && !l->traverseVisible(ray, triangles, maxT)) return false;
		if (r && !r->traverseVisible(ray, triangles, maxT)) return false;
		return true;
	}
};
