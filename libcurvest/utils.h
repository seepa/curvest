/*
 * Copyright (C) 2016, Patrick Seemann
 *
 * This software may be modified and distributed under the terms
 * of the BSD 3-Clause license. See the LICENSE.txt file for details.
 */

#pragma once

#include <vector>
#include <iostream>
#include <set>

#include "mve/mesh.h"
#include "mve/mesh_info.h"
#include "math/vector.h"

#include "defines.h"

NAMESPACE_CURVEST_BEGIN
NAMESPACE_UTILS_BEGIN

struct Result {
    unsigned int vertexId;
    std::string name;
    std::string fittedFunction;

    // pair.first is the set option string
    // pair.second is the corresponding unset option string
    std::vector<std::pair<std::string, std::string> > setOptions;
    std::vector<math::Vec2f> data;
    bool useCustomLabels;
    std::string xLabel;
    std::string yLabel;
};

// Define a 'less than' struct for sorting Result structs
struct Result_less_than_key {
    inline bool operator()(math::Vec2f const& r1, math::Vec2f const& r2)
    {
        return (r1[0] < r2[0]);
    }
};

// Intersects a ray with a sphere
// http://www.csee.umbc.edu/~olano/435f02/ray-sphere.html
//
// Returns false if there is no intersection, otherwise true.
// If the return value is true, "intersections" holds either one or
// two intersection points.
bool raySphereIntersection(math::Vec3f point, math::Vec3f direction,
    math::Vec3f center, float radius, std::vector<math::Vec3f>& intersections);

// Helper: Intersects a ray with a sphere, but the given point
// must lay inside the sphere.
//
// Returns the intersection point, which lays in the given direction
math::Vec3f raySphereIntersectionInside(
    math::Vec3f point, math::Vec3f direction, math::Vec3f center, float radius);

// Fit cubic polynomial to the given data
bool fitCubicPolynomial(std::vector<math::Vec2f> const& data,
    math::Vec4f& coeffsOut, float& errorOut);

// Get a percentile from the given data.
// The input vector will be copied for calculating the percentile.
//
// percent: The percentage of the percentile, e.g. 1.0f for the
//          1-percentile or 90.0f for the 90-percentile
float percentile(std::vector<float> const& vectorIn, float percent);

// Clamps the given data vector to the given min and max percentile
void percentileClamp(
    std::vector<float>& data, float percentLow, float percentHigh);

// Smoothes the values in data (i.e. vertex values) based on the
// connectivity of the given mesh.
void smoothMeshProperty(mve::TriangleMesh::ConstPtr mesh,
    std::vector<float>& data, std::size_t iters, float lambda);

// Get a pretty string from millisecons time
std::string getPrettyTimeString(std::size_t time);

// Calculates extreme and inflection points of the given cubic polynomial.
void getInfoCubicPoly(math::Vec4f const& coeffs, float& x1, float& x2,
    float& inflectionPoint, bool outputDebugInfo);

NAMESPACE_UTILS_END
NAMESPACE_CURVEST_END
