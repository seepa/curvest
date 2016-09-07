/*
 * Copyright (C) 2016, Patrick Seemann
 *
 * This software may be modified and distributed under the terms
 * of the BSD 3-Clause license. See the LICENSE.txt file for details.
 */

#pragma once

#include <set>

#include "mve/mesh.h"
#include "mve/mesh_info.h"

#include "defines.h"
#include "curvest.h"

NAMESPACE_CURVEST_BEGIN
/*
 * Performs a circular region growing to create a surface patch. Additional
 * triangles will be added to the patch so that the patch border is closer
 * to the radius and thus more circular.
 */
class RegionGrower
{
    struct Edge {
        Edge(std::size_t idx0, std::size_t idx1) : a(idx0), b(idx1) {}
        std::size_t a, b;
    };

   public:

    /*
     * Creates a region grower centered at the vertex with id 'vertId'.
     */
    RegionGrower(std::size_t const vertId, mve::TriangleMesh::ConstPtr mesh,
        mve::MeshInfo const* vertexInfoList, Estimator::Config const& conf);

    /*
     * Returns a surface patch of size 'radius'.
     */
    void getSurfacePatch(float const radius, mve::TriangleMesh::Ptr patchMesh);

    /*
     * Returns true iff the the edge defined by idx{0,1} is a boundary edge
     * of the mesh.
     */
    bool isBorderEdge(std::size_t idx0, std::size_t idx1,
        mve::TriangleMesh::ConstPtr patchMesh, mve::MeshInfo* vInfoList,
        const float& radius);

    /*
     * Adds a triangle to the mesh. The given edge is used as one side. The
     * third vertex is placed on the circle.
     */
    bool performFlap(
        const Edge& e, mve::TriangleMesh::Ptr patchMesh, float const& radius);

   private:
    std::size_t const centerVertId;
    mve::TriangleMesh::ConstPtr mesh;
    mve::MeshInfo const* vertInfoList;
    Estimator::Config const& conf;
};

NAMESPACE_CURVEST_END
