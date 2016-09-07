/*
 * Copyright (C) 2016, Patrick Seemann
 *
 * This software may be modified and distributed under the terms
 * of the BSD 3-Clause license. See the LICENSE.txt file for details.
 */

#pragma once

#include <iostream>
#include <set>
#include <vector>

#define NAMESPACE_CURVEST_BEGIN \
    namespace curvest           \
    {
#define NAMESPACE_CURVEST_END }

#define NAMESPACE_UTILS_BEGIN \
    namespace utils           \
    {
#define NAMESPACE_UTILS_END }

NAMESPACE_CURVEST_BEGIN

// Default parameters
#define DEFAULT_SUBDIV_ITERS 2
#define DEFAULT_RADIUS_FACTOR_INITIAL 1.0f
#define DEFAULT_RADIUS_FACTOR_INC 1.3f
#define DEFAULT_RADIUS_FACTOR_MAX 10.0f
#define DEFAULT_EDGE_SMOOTH 0.2f
#define DEFAULT_PLANAR_THRESHOLD 0.2f
#define DEFAULT_INITIAL_RADIUS_SMOOTH_ITERS 5

// Dimensions of gnuplot graphs (for debugging)
#define GRAPH_HEIGHT 300
#define GRAPH_WIDTH 350

// Smooth final radius before computing final mean curvature?
#define SMOOTH_FINAL_RADIUS 1
#define FINAL_RADIUS_SMOOTHING_ITERS 2

// Clamp final mean curvature using percentiles?
#define CLAMP_FINAL_CURVATURE 1
#define PERCENTILE_LOW 1.0f    // %
#define PERCENTILE_HIGH 99.0f  // %

// Number of subdivision iterations of the sphere when computing
// the adaptive intersection with the mesh
#define INTERSECTION_SUBDIV_ITERS 6

// Helper macros
#define EPSILON_EQUAL(a, b, e) ((std::abs(a - b) < (e)))
#define IS_NAN(x) ((x) != (x))

NAMESPACE_CURVEST_END
