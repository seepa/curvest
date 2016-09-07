/*
 * Copyright (C) 2016, Patrick Seemann
 *
 * This software may be modified and distributed under the terms
 * of the BSD 3-Clause license. See the LICENSE.txt file for details.
 */

#include <iostream>
#include <sstream>
#include <vector>
#include <stdio.h>
#include <stdlib.h>
#include <set>
#include <cassert>
#include <exception>

#include "utils.h"
#include "defines.h"

#include "mve/mesh.h"
#include "mve/mesh_info.h"
#include "math/matrix.h"
#include "math/matrix_tools.h"

NAMESPACE_CURVEST_BEGIN
NAMESPACE_UTILS_BEGIN

bool raySphereIntersection(math::Vec3f point, math::Vec3f direction,
    math::Vec3f center, float radius, std::vector<math::Vec3f>& intersections)
{
    // make sure direction is normalized
    assert(EPSILON_EQUAL(std::abs(direction.norm()), 1.0, 1e-5) &&
           "[raySphereIntersection] Ray direction should be normalized!");

    if (EPSILON_EQUAL(std::abs((point - center).norm()), 0.0, 1e-13)) {
        math::Vec3f p1 = point + radius * direction;
        intersections.push_back(p1);
        return true;
    }

    double b, c, discriminant;
    b = 2.0 * direction.dot(point - center);
    c = ((point - center).dot(point - center)) - radius * radius;
    discriminant = b * b - 4.0 * c;
    if (discriminant < 0.0) {
        // no itersection
        return false;
    } else if (EPSILON_EQUAL(std::abs(discriminant), 0.0, 1e-13)) {
        // one intersection
        std::cout << 1.0 - discriminant << std::endl;
        std::cout << "Warning: [raySphereIntersection] Only one "
                  << "intersection point found!" << std::endl;
        std::exit(1);
        float t = (-1.0 * b) / 2.0;
        math::Vec3f p1 = point + t * direction;
        intersections.push_back(p1);
        return true;
    } else {
        // two intersections
        float t1, t2;
        math::Vec3f p1, p2;
        t1 = (-1.0 * b + std::sqrt(discriminant)) / 2.0;
        t2 = (-1.0 * b - std::sqrt(discriminant)) / 2.0;
        p1 = point + t1 * direction;
        p2 = point + t2 * direction;
        intersections.push_back(p1);
        intersections.push_back(p2);
        return true;
    }
}

math::Vec3f raySphereIntersectionInside(
    math::Vec3f point, math::Vec3f direction, math::Vec3f center, float radius)
{
    std::vector<math::Vec3f> intersections;
    bool res =
        raySphereIntersection(point, direction, center, radius, intersections);
    assert(res &&
           "Ray assumed to be inside sphere, but no intersection was found!");

    // return the intersection point lays in the given direction
    if (intersections.size() == 1) {
        return intersections[0];
    } else {
        if ((intersections[0] - point).norm() >= 0.0f)
            return intersections[0];
        else
            return intersections[1];
    }
}

bool fitCubicPolynomial(std::vector<math::Vec2f> const& data,
    math::Vec4f& coeffsOut, float& errorOut)
{
    // perform LEAST SQUARES fitting

    // setup matrices
    std::vector<double> A, ATrans, b;
    for (std::size_t i = 0; i < data.size(); ++i) {
        double x, y;
        x = data[i][0];
        y = data[i][1];

        A.push_back(1);
        A.push_back(x);
        A.push_back(x * x);
        A.push_back(x * x * x);

        b.push_back(y);
    }

    // transpose A
    for (std::size_t i = 0; i < A.size(); i += 4) ATrans.push_back(A[i]);
    for (std::size_t i = 1; i < A.size(); i += 4) ATrans.push_back(A[i]);
    for (std::size_t i = 2; i < A.size(); i += 4) ATrans.push_back(A[i]);
    for (std::size_t i = 3; i < A.size(); i += 4) ATrans.push_back(A[i]);

    // ATrans * b
    std::vector<double> M1(4, 0.0);
    math::matrix_multiply(&ATrans[0], 4, ATrans.size() / 4, &b[0], 1, &M1[0]);
    math::Vec4f ATrans_b;
    for (std::size_t i = 0; i < 4; ++i) ATrans_b[i] = M1[i];

    // ATrans * A
    std::vector<double> M2(16, 0.0);
    math::matrix_multiply(&ATrans[0], 4, ATrans.size() / 4, &A[0], 4, &M2[0]);

    math::Matrix4d M3;
    for (std::size_t i = 0; i < 4; ++i) {
        M3(i, 0) = M2[i * 4];
        M3(i, 1) = M2[i * 4 + 1];
        M3(i, 2) = M2[i * 4 + 2];
        M3(i, 3) = M2[i * 4 + 3];
    }

    // (ATrans * A) ^ -1
    math::Matrix4d ATrans_A_Inv = math::matrix_inverse(M3);

    math::Vec4f C;
    C = ATrans_A_Inv.mult(ATrans_b);

    // calculate fit error
    errorOut = 0.0f;
    float min, max;
    min = std::numeric_limits<float>::max();
    max = -std::numeric_limits<float>::max();
    for (std::size_t i = 0; i < data.size(); ++i) {
        float x = data[i][0];
        float y = data[i][1];
        float fittedY = C[0] + C[1] * x + C[2] * x * x + C[3] * x * x * x;
        float dist = std::abs(y - fittedY);
        errorOut += dist;

        if (y < min) min = y;
        if (y > max) max = y;
    }
    errorOut /= data.size();
    errorOut /= std::abs(max - min);
    coeffsOut = C;
    return true;
}

float percentile(std::vector<float> const& vectorIn, float percent)
{
    // copy the vectorIn because nth_element will rearrange it
    std::vector<float> vector(vectorIn.size());
    for (std::size_t i = 0; i < vector.size(); ++i) vector[i] = vectorIn[i];
    auto nth = vector.begin() + (percent * vector.size()) / 100;
    std::nth_element(vector.begin(), nth, vector.end());
    return *nth;
}

void smoothMeshProperty(mve::TriangleMesh::ConstPtr mesh,
    std::vector<float>& data, std::size_t iters, float lambda)
{
    mve::TriangleMesh::VertexList const& mVerts = mesh->get_vertices();
    mve::MeshInfo const vinfo(mesh);
    assert(data.size() == mVerts.size());

    for (std::size_t iter = 0; iter < iters; ++iter) {
        // backup
        std::vector<float> backup(mVerts.size(), 0.0f);
        for (std::size_t i = 0; i < mVerts.size(); ++i) backup[i] = data[i];

        // smooth
        for (std::size_t i = 0; i < mVerts.size(); ++i) {
            mve::MeshInfo::VertexInfo const& vinf = vinfo[i];
            float sumValues = 0.0f;
            if (vinf.verts.size() <= 2) continue;

            std::size_t validNeighbors = 0;
            for (std::size_t k = 0; k < vinf.verts.size(); ++k) {
                if (backup[vinf.verts[k]] != -1.0f) {
                    sumValues += backup[vinf.verts[k]] - backup[i];
                    validNeighbors++;
                }
            }
            data[i] = backup[i] + lambda * sumValues / (float)validNeighbors;
        }
    }
}

void percentileClamp(
    std::vector<float>& data, float percentLow, float percentHigh)
{
    float minPercentile = utils::percentile(data, percentLow);
    float maxPercentile = utils::percentile(data, percentHigh);
    std::cout << "low percentile: " << minPercentile
              << " high percentile: " << maxPercentile << std::endl;

    for (std::size_t i = 0; i < data.size(); ++i) {
        if (data[i] < minPercentile) data[i] = minPercentile;
        if (data[i] > maxPercentile) data[i] = maxPercentile;
    }
}

std::string getPrettyTimeString(std::size_t time)
{
    std::size_t ms, m, s;
    ms = time % 1000;
    m = (time / 1000) / 60;
    s = (time / 1000) % 60;

    std::stringstream ss;
    if (m != 0) ss << m << "m ";
    if (s != 0) ss << s << "s ";
    if (m == 0) ss << ms << "ms";

    return ss.str();
}

void getInfoCubicPoly(math::Vec4f const& coeffs, float& x1, float& x2,
    float& inflectionPoint, bool outputDebugInfo = false)
{
    float p = 2.0f * coeffs[2] / (6.0f * coeffs[3]);
    float q = coeffs[1] / (3.0f * coeffs[3]);
    float descriminant = p * p - q;

    inflectionPoint = (-2.0f * coeffs[2]) / (6.0f * coeffs[3]);

    // check inflection point
    float p1, p2;
    p1 = 6.0f * coeffs[3] * (inflectionPoint + 0.1f) + 2.0f * coeffs[2];
    p2 = 6.0f * coeffs[3] * (inflectionPoint - 0.1f) + 2.0f * coeffs[2];

    if (p1 * p2 < 0.0f) {
        if (outputDebugInfo)
            std::cout << "Inflection: " << inflectionPoint << std::endl;
    } else {
        if (outputDebugInfo)
            std::cout << "No Inflection point" << std::endl;
    }

    if (EPSILON_EQUAL(descriminant, 0.0f, 1e-13)) {
        std::cout << "Warning: descriminant is zero!" << std::endl;
    }
    // extrema points
    x1 = -p + std::sqrt(descriminant);
    x2 = -p - std::sqrt(descriminant);

    // make sure x1 <= x2
    if (x1 > x2) {
        float tmp = x1;
        x1 = x2;
        x2 = tmp;
    }
}

NAMESPACE_UTILS_END
NAMESPACE_CURVEST_END
