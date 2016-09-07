/*
 * Copyright (C) 2016, Patrick Seemann
 *
 * This software may be modified and distributed under the terms
 * of the BSD 3-Clause license. See the LICENSE.txt file for details.
 */

#pragma once

#include <limits.h>
#include <vector>
#include <cassert>

#include "defines.h"
#include "mve/mesh.h"
#include "mve/mesh_info.h"

NAMESPACE_CURVEST_BEGIN
NAMESPACE_UTILS_BEGIN

// define NotAnIndex
#define NAI std::numeric_limits<std::size_t>::max()

inline void getOppositeVertices(mve::TriangleMesh::ConstPtr mesh,
    mve::MeshInfo const& mInfo, std::size_t const& v0Id,
    std::size_t const& v1Id, std::size_t& f1Id, std::size_t& f2Id,
    std::size_t& oId1, std::size_t& oId2)
{
    mve::TriangleMesh::FaceList const& faces = mesh->get_faces();
    std::vector<std::size_t> adjFaces;
    mInfo.get_faces_for_edge(v0Id, v1Id, &adjFaces);

    assert(
        adjFaces.size() == 2 && "Edge must have exactly two adjacent faces!");

    // the face which contains the edge v0Id -> v1Id is face1
    f1Id = adjFaces[0];
    f2Id = adjFaces[1];

    // check if this decision was wrong, that is, the edge v0Id -> v1Id is
    // actually in the other face
    for (std::size_t k = 0; k < 3; ++k) {
        std::size_t nextId = (k + 1) % 3;
        if (faces[adjFaces[1] * 3 + k] == v0Id &&
            faces[adjFaces[1] * 3 + nextId] == v1Id) {
            // switch the faces
            f1Id = adjFaces[1];
            f2Id = adjFaces[0];
            break;
        }
    }

    // get two vertices opposite of the edge
    oId1 = NAI;
    oId2 = NAI;
    for (std::size_t k = 0; k < 3; ++k) {
        std::size_t const& id = faces[f1Id * 3 + k];
        if (id != v0Id && id != v1Id) {
            oId1 = id;
            break;
        }
    }
    for (std::size_t k = 0; k < 3; ++k) {
        std::size_t const& id = faces[f2Id * 3 + k];
        if (id != v0Id && id != v1Id) {
            oId2 = id;
            break;
        }
    }
    assert((oId1 != NAI && oId2 != NAI) && "Could not find opposite vertices.");
}



// Loop subdivision for spheres only !!
//
// Each old face will be subdivided into four new faces: Three new
// vertices will be added per face, that is, one per edge.
//
// The assumption here is that we start with an icosahedron/sphere, so there
// is no mesh boundary or other special cases that have to be taken care of.
//
// The only case we care about is an edge (v,vn) that has exactly two
// neighboring faces:
//
//            o1
//           /   \
//          / f1  \
//         v------vn
//          \ f2  /
//           \   /
//            o2
//
// A new vertex will be created on the edge (v, vn).
//
inline void loopSubdivide(mve::TriangleMesh::Ptr& mesh)
{
    mve::TriangleMesh::VertexList& verts = mesh->get_vertices();
    mve::TriangleMesh::FaceList& faces = mesh->get_faces();
    std::size_t const numVerts = verts.size();
    std::size_t const numFaces = faces.size() / 3;
    mve::MeshInfo mInfo(mesh);
    std::vector<std::size_t> newVertIds(numFaces * 3, NAI);

    for (std::size_t i = 0; i < numFaces; ++i) {
        std::size_t const fId = i * 3;

        for (std::size_t j = 0; j < 3; ++j) {
            std::size_t const nextId = (j + 1) % 3;

            // current edge (v, vn)
            std::size_t const vId = faces[fId + j];
            std::size_t const vnId = faces[fId + nextId];

            // get ids of opposite vertices o1 and o2
            std::size_t f1Id, f2Id, o1Id, o2Id;
            getOppositeVertices(mesh, mInfo, vId, vnId, f1Id, f2Id, o1Id, o2Id);

            // skip if we have seen this edge already
            if (newVertIds[f1Id * 3 + j] != NAI) continue;

            math::Vec3f const& v = verts[vId];
            math::Vec3f const& vn = verts[vnId];
            math::Vec3f const& vo1 = verts[o1Id];
            math::Vec3f const& vo2 = verts[o2Id];

            // calculate and add new edge vertex
            float a = 3.0f / 8.0f;
            float b = 1.0f / 8.0f;
            math::Vec3f newVert = a * (v + vn) + b * (vo1 + vo2);
            verts.push_back(newVert);

            // save new vertex id
            std::size_t offset = NAI;
            for (std::size_t k = 0; k < 3; ++k) {
                if (faces[f2Id * 3 + k] == vnId) {
                    offset = k;
                    break;
                }
            }
            assert(offset != NAI);
            assert(newVertIds[f1Id * 3 + j] == NAI &&
                   newVertIds[f2Id * 3 + offset] == NAI);

            newVertIds[f1Id * 3 + j] = verts.size() - 1;
            newVertIds[f2Id * 3 + offset] = verts.size() - 1;
        }
    }

    // calculate new positions of old vertices
    // the new position of a vertex depends on all its surrounding vertices
    for (std::size_t i = 0; i < numVerts; ++i) {
        // loop over adjacent vertices
        mve::MeshInfo::VertexInfo const& vInfo = mInfo[i];
        math::Vec3f w(0.0f);
        std::size_t const N = vInfo.verts.size();
        for (std::size_t j = 0; j < N; ++j) {
            w += verts[vInfo.verts[j]];
        }
        assert(N >= 3);
        float beta;
        if (N == 3)
            beta = 3.0f / 16.0f;
        else if (N >= 3)
            beta = 3.0f / (N * 8.0f);

        verts[i] = (1.0f - N * beta) * verts[i] + beta * w;
    }

    // now subdivide each face into four new faces
    for (std::size_t i = 0; i < numFaces; ++i) {
        std::size_t const fId = i * 3;

        // original vertices
        std::size_t vIdx0 = faces[fId];
        std::size_t vIdx1 = faces[fId + 1];
        std::size_t vIdx2 = faces[fId + 2];

        // added vertices on each of the three edges
        std::size_t wIdx0 = newVertIds[fId];
        std::size_t wIdx1 = newVertIds[fId + 1];
        std::size_t wIdx2 = newVertIds[fId + 2];

        assert(wIdx0 != NAI && wIdx1 != NAI && wIdx2 != NAI);

        faces[fId] = vIdx0;
        faces[fId + 1] = wIdx0;
        faces[fId + 2] = wIdx2;

        faces.push_back(wIdx0);
        faces.push_back(vIdx1);
        faces.push_back(wIdx1);

        faces.push_back(wIdx1);
        faces.push_back(vIdx2);
        faces.push_back(wIdx2);

        faces.push_back(wIdx0);
        faces.push_back(wIdx1);
        faces.push_back(wIdx2);
    }
}

inline mve::TriangleMesh::Ptr getSphere(std::size_t const& iters)
{
    // create an icosahedron
    mve::TriangleMesh::Ptr sphere = mve::TriangleMesh::create();
    mve::TriangleMesh::VertexList& verts = sphere->get_vertices();
    verts.push_back(math::Vec3f(0.0f, -0.525731f, 0.850651f));
    verts.push_back(math::Vec3f(0.0f, 0.525731f, 0.850651f));
    verts.push_back(math::Vec3f(0.0f, -0.525731f, -0.850651f));
    verts.push_back(math::Vec3f(0.0f, 0.525731f, -0.850651f));
    verts.push_back(math::Vec3f(0.850651f, 0.0f, 0.525731f));
    verts.push_back(math::Vec3f(0.850651f, 0.0f, -0.525731f));
    verts.push_back(math::Vec3f(-0.850651f, 0.0f, 0.525731f));
    verts.push_back(math::Vec3f(-0.850651f, 0.0f, -0.525731f));
    verts.push_back(math::Vec3f(0.525731f, 0.850651f, 0.0f));
    verts.push_back(math::Vec3f(0.525731f, -0.850651f, 0.0f));
    verts.push_back(math::Vec3f(-0.525731f, 0.850651f, 0.0f));
    verts.push_back(math::Vec3f(-0.525731f, -0.850651f, 0.0f));

    mve::TriangleMesh::FaceList& faces = sphere->get_faces();
    faces = {0, 4, 1, 0, 9, 4, 9, 5, 4, 4, 5, 8, 4, 8, 1, 8, 10, 1, 8, 3, 10, 5,
        3, 8, 5, 2, 3, 2, 7, 3, 7, 10, 3, 7, 6, 10, 7, 11, 6, 11, 0, 6, 0, 1, 6,
        6, 1, 10, 9, 0, 11, 9, 11, 2, 9, 2, 5, 7, 2, 11};

    // subdivide icosahedron using Loop subdivision
    for (std::size_t i = 0; i < iters; ++i) loopSubdivide(sphere);

    // scale sphere to radius 1
    for (std::size_t i = 0; i < verts.size(); ++i)
        verts[i] = (verts[i] / verts[i].norm());

    return sphere;
}

NAMESPACE_UTILS_END
NAMESPACE_CURVEST_END
