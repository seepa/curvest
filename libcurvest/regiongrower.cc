/*
 * Copyright (C) 2016, Patrick Seemann
 *
 * This software may be modified and distributed under the terms
 * of the BSD 3-Clause license. See the LICENSE.txt file for details.
 */

#include "regiongrower.h"
#include <iostream>
#include <map>

#include <set>

#include "mve/mesh.h"
#include "mve/mesh_info.h"
#include "mve/mesh_io_ply.h"

#include "utils.h"

NAMESPACE_CURVEST_BEGIN

// When a vertex is this distance away from the radius, still use it to
// triangulate the patch. If the vertex is farther away, the triangle containing
// this vertex is cut and retriangulated to fit inside the patch radius
#define RADIUS_EPSILON_CHECK 1e-8

// Min length an edge needs to have in order to perform a flap on that edge
#define MIN_EDGE_LENGTH 1e-8

// Min distance from a patch border edge to the radius
#define MIN_EDGE_TO_RADIUS_DIST 1e-5

typedef mve::MeshBase::VertexID VertexId;

RegionGrower::RegionGrower(std::size_t const vertId,
    mve::TriangleMesh::ConstPtr mesh, const mve::MeshInfo* vertexInfoList,
    Estimator::Config const& conf)
    : centerVertId(vertId), mesh(mesh), vertInfoList(vertexInfoList), conf(conf)
{
}

void RegionGrower::getSurfacePatch(
    float const radius, mve::TriangleMesh::Ptr patchMesh)
{
    mve::MeshBase::VertexList const& mVerts = this->mesh->get_vertices();
    mve::TriangleMesh::FaceList const& mFaces = this->mesh->get_faces();

    std::set<std::size_t> vertexIdsToVisit;
    vertexIdsToVisit.insert(this->centerVertId);

    std::set<std::size_t> visitedVertIds;
    mve::TriangleMesh::FaceList& newFaces = patchMesh->get_faces();
    mve::TriangleMesh::VertexList& newVertices = patchMesh->get_vertices();

    // maps old vertex ids to new vertex ids
    std::map<VertexId, VertexId> oldToNewIds;

    // keeps track of new verts inserted on an previous edge
    std::map<std::pair<VertexId, VertexId>, VertexId> newVertIds;

    // if we are on the mesh border, we can't perform the patch-refinement
    // at the end
    bool seenBorderVertex = false;
    while (!vertexIdsToVisit.empty()) {
        std::size_t currVertexId = *(vertexIdsToVisit.begin());
        vertexIdsToVisit.erase(vertexIdsToVisit.begin());
        visitedVertIds.insert(currVertexId);

        mve::MeshInfo::VertexInfo const& vertInfo =
            (*this->vertInfoList)[currVertexId];
        mve::MeshInfo::AdjacentFaces const& adjFaces = vertInfo.faces;

        for (std::size_t i = 0; i < adjFaces.size(); ++i) {
            // figure out ids of both neighbor verts
            std::size_t foff = adjFaces[i] * 3;
            std::size_t k = 0;
            for (; k < 3; ++k)
                if (mFaces[foff + k] == currVertexId) break;
            std::size_t neighborId, nextNeighborId;
            neighborId = mFaces[foff + (k + 1) % 3];
            nextNeighborId = mFaces[foff + (k + 2) % 3];

            // check if neighor or nextneighbor is on a border
            {
                mve::MeshInfo::VertexInfo const& neighborInfo =
                    this->vertInfoList->at(neighborId);
                mve::MeshInfo::VertexInfo const& nextNeighborInfo =
                    this->vertInfoList->at(neighborId);

                if (neighborInfo.vclass ==
                        mve::MeshInfo::VertexClass::VERTEX_CLASS_BORDER ||
                    nextNeighborInfo.vclass ==
                        mve::MeshInfo::VertexClass::VERTEX_CLASS_BORDER) {
                    seenBorderVertex = true;
                }
            }

            math::Vec3f v0, v1, v2, n;
            v0 = mVerts[currVertexId];
            v1 = mVerts[neighborId];
            v2 = mVerts[nextNeighborId];
            n = (v1 - v0).cross(v2 - v0);
            if (n.norm() == 0) continue;

            if (visitedVertIds.find(neighborId) != visitedVertIds.end())
                continue;
            if (visitedVertIds.find(nextNeighborId) != visitedVertIds.end())
                continue;

            bool neighborInside, nextNeighborInside;
            neighborInside = nextNeighborInside = false;

            double dist1 =
                (mVerts[this->centerVertId] - mVerts[neighborId]).norm();
            double check1 = dist1 - (double)radius;
            bool useNeighbor = false;

            // check if vertex far enough inside the radius
            if (check1 <= -RADIUS_EPSILON_CHECK)
                neighborInside = true;
            else {
                // check if it was within the epsilon range
                if (std::abs(check1) >= 0.0 &&
                    std::abs(check1) <= RADIUS_EPSILON_CHECK)
                    useNeighbor = true;
            }

            double dist2 =
                (mVerts[this->centerVertId] - mVerts[nextNeighborId]).norm();
            bool useNextNeighbor = false;
            double check2 = dist2 - (double)radius;
            if (check2 <= -RADIUS_EPSILON_CHECK)
                nextNeighborInside = true;
            else {
                if (std::abs(check2) >= 0.0 &&
                    std::abs(check2) <= RADIUS_EPSILON_CHECK)
                    useNextNeighbor = true;
            }

            // case 1: triangle fully inside
            if (neighborInside && nextNeighborInside) {
                std::size_t idx0, idx1, idx2;

                auto res0 = oldToNewIds.insert(
                    std::make_pair(currVertexId, newVertices.size()));

                if (res0.second) {
                    idx0 = newVertices.size();
                    newVertices.push_back(mVerts[currVertexId]);
                } else
                    idx0 = res0.first->second;

                auto res1 = oldToNewIds.insert(
                    std::make_pair(neighborId, newVertices.size()));
                if (res1.second) {
                    idx1 = newVertices.size();
                    newVertices.push_back(mVerts[neighborId]);
                } else
                    idx1 = res1.first->second;

                auto res2 = oldToNewIds.insert(
                    std::make_pair(nextNeighborId, newVertices.size()));
                if (res2.second) {
                    idx2 = newVertices.size();
                    newVertices.push_back(mVerts[nextNeighborId]);
                } else
                    idx2 = res2.first->second;

                newFaces.push_back(idx0);
                newFaces.push_back(idx1);
                newFaces.push_back(idx2);

                vertexIdsToVisit.insert(neighborId);
                vertexIdsToVisit.insert(nextNeighborId);

                continue;
            }

            // case2: triangle has two vertices outside
            if (!(neighborInside || nextNeighborInside)) {
                // create two new faces
                math::Vec3f dirNeighbor, dirNextNeighbor;
                dirNeighbor = mVerts[neighborId] - mVerts[currVertexId];
                dirNextNeighbor = mVerts[nextNeighborId] - mVerts[currVertexId];

                dirNeighbor /= dirNeighbor.norm();
                dirNextNeighbor /= dirNextNeighbor.norm();

                // either calculate two new vertices or reuse the outside
                // ones, if they were really close to the radius
                math::Vec3f ip1, ip2;
                if (useNeighbor)
                    ip1 = mVerts[neighborId];
                else
                    ip1 =
                        utils::raySphereIntersectionInside(mVerts[currVertexId],
                            dirNeighbor, mVerts[this->centerVertId], radius);

                if (useNextNeighbor)
                    ip2 = mVerts[nextNeighborId];
                else
                    ip2 = utils::raySphereIntersectionInside(
                        mVerts[currVertexId], dirNextNeighbor,
                        mVerts[this->centerVertId], radius);

                std::size_t idx0, idx1, idx2;
                auto res0 = oldToNewIds.insert(
                    std::make_pair(currVertexId, newVertices.size()));
                if (res0.second) {
                    idx0 = newVertices.size();
                    newVertices.push_back(mVerts[currVertexId]);
                } else
                    idx0 = res0.first->second;

                auto edgeNeighbor = std::make_pair(currVertexId, neighborId);
                auto res1 = newVertIds.insert(
                    std::make_pair(edgeNeighbor, newVertices.size()));
                if (res1.second) {
                    idx1 = newVertices.size();
                    newVertices.push_back(ip1);
                } else
                    idx1 = res1.first->second;

                auto edgeNextNeighbor =
                    std::make_pair(currVertexId, nextNeighborId);
                auto res2 = newVertIds.insert(
                    std::make_pair(edgeNextNeighbor, newVertices.size()));
                if (res2.second) {
                    idx2 = newVertices.size();
                    newVertices.push_back(ip2);
                } else
                    idx2 = res2.first->second;

                newFaces.push_back(idx0);
                newFaces.push_back(idx1);
                newFaces.push_back(idx2);
                continue;
            }

            // case 3: triangle has one vertex outside
            if (neighborInside || nextNeighborInside) {
                std::size_t outside, inside1, inside2;
                if (neighborInside) {
                    outside = nextNeighborId;
                    inside1 = currVertexId;
                    inside2 = neighborId;
                    vertexIdsToVisit.insert(neighborId);
                } else {
                    outside = neighborId;
                    inside1 = nextNeighborId;
                    inside2 = currVertexId;
                    vertexIdsToVisit.insert(nextNeighborId);
                }

                math::Vec3f dirInside1, dirInside2;
                dirInside1 = mVerts[outside] - mVerts[inside1];
                dirInside2 = mVerts[outside] - mVerts[inside2];

                dirInside1 /= dirInside1.norm();
                dirInside2 /= dirInside2.norm();

                // either calculate two new vertices or reuse the outside
                // vertex, if it was really close to the radius
                math::Vec3f ip1, ip2;
                if ((neighborId == outside && useNeighbor) ||
                    nextNeighborId == outside && useNextNeighbor) {
                    ip1 = mVerts[outside];
                } else {
                    ip1 = utils::raySphereIntersectionInside(mVerts[inside1],
                        dirInside1, mVerts[this->centerVertId], radius);
                    ip2 = utils::raySphereIntersectionInside(mVerts[inside2],
                        dirInside2, mVerts[this->centerVertId], radius);
                }

                // just add one triangle
                if (useNeighbor || useNextNeighbor) {
                    std::size_t idx0;  // outside
                    std::size_t idx1;  // inside1
                    std::size_t idx2;  // inside2

                    auto res0 = oldToNewIds.insert(
                        std::make_pair(outside, newVertices.size()));

                    if (res0.second) {
                        idx0 = newVertices.size();
                        newVertices.push_back(ip1);
                    } else
                        idx0 = res0.first->second;

                    auto res1 = oldToNewIds.insert(
                        std::make_pair(inside1, newVertices.size()));
                    if (res1.second) {
                        idx1 = newVertices.size();
                        newVertices.push_back(mVerts[inside1]);
                    } else
                        idx1 = res1.first->second;

                    auto res2 = oldToNewIds.insert(
                        std::make_pair(inside2, newVertices.size()));
                    if (res2.second) {
                        idx2 = newVertices.size();
                        newVertices.push_back(mVerts[inside2]);
                    } else
                        idx2 = res2.first->second;

                    newFaces.push_back(idx0);
                    newFaces.push_back(idx1);
                    newFaces.push_back(idx2);
                    continue;
                }

                // else: add two triangles
                std::size_t idx0, idx1, idx2, idx3;
                auto edgeInside1 = std::make_pair(inside1, outside);
                auto res0 = newVertIds.insert(
                    std::make_pair(edgeInside1, newVertices.size()));
                if (res0.second) {
                    idx3 = newVertices.size();
                    newVertices.push_back(ip1);
                } else
                    idx3 = res0.first->second;

                auto edgeInside2 = std::make_pair(inside2, outside);
                auto res1 = newVertIds.insert(
                    std::make_pair(edgeInside2, newVertices.size()));
                if (res1.second) {
                    idx2 = newVertices.size();
                    newVertices.push_back(ip2);
                } else
                    idx2 = res1.first->second;

                auto res2 = oldToNewIds.insert(
                    std::make_pair(inside1, newVertices.size()));
                if (res2.second) {
                    idx0 = newVertices.size();
                    newVertices.push_back(mVerts[inside1]);
                } else
                    idx0 = res2.first->second;

                auto res3 = oldToNewIds.insert(
                    std::make_pair(inside2, newVertices.size()));
                if (res3.second) {
                    idx1 = newVertices.size();
                    newVertices.push_back(mVerts[inside2]);
                } else
                    idx1 = res3.first->second;

                newFaces.push_back(idx1);
                newFaces.push_back(idx2);
                newFaces.push_back(idx3);

                newFaces.push_back(idx1);
                newFaces.push_back(idx3);
                newFaces.push_back(idx0);
            }
        }
    }

    // refine the patch if it doesn't cross the mesh border
    if (!seenBorderVertex) {
        // get all border edges
        mve::MeshInfo oldVInfoList(patchMesh);
        oldVInfoList.initialize(patchMesh);
        mve::TriangleMesh::FaceList const& patchFaces = patchMesh->get_faces();
        std::vector<Edge> borderEdges;
        for (std::size_t i = 0; i < patchFaces.size(); i = i + 3) {
            std::size_t const& vidx0 = patchFaces[i];
            std::size_t const& vidx1 = patchFaces[i + 1];
            std::size_t const& vidx2 = patchFaces[i + 2];

            if (isBorderEdge(vidx0, vidx1, patchMesh, &oldVInfoList, radius)) {
                borderEdges.push_back(Edge(vidx0, vidx1));
            }
            if (isBorderEdge(vidx1, vidx2, patchMesh, &oldVInfoList, radius)) {
                borderEdges.push_back(Edge(vidx1, vidx2));
            }
            if (isBorderEdge(vidx2, vidx0, patchMesh, &oldVInfoList, radius)) {
                borderEdges.push_back(Edge(vidx2, vidx0));
            }
        }

        // Add triangles to the border until all borderEdges are close enough
        //to the radius
        while (!borderEdges.empty()) {
            // get next edge
            Edge edge = borderEdges[borderEdges.size() - 1];
            borderEdges.pop_back();

            // save next vertex id, if a triangle will be created
            std::size_t newVertId = patchMesh->get_vertices().size();
            if (performFlap(edge, patchMesh, radius)) {
                // a new triangle has been added
                borderEdges.push_back(Edge(edge.a, newVertId));
                borderEdges.push_back(Edge(newVertId, edge.b));
            }
        }
    }
}

bool RegionGrower::isBorderEdge(std::size_t idx0, std::size_t idx1,
    mve::TriangleMesh::ConstPtr patchMesh, mve::MeshInfo* vInfoList,
    float const& radius)
{
    std::vector<std::size_t> adjFaceIds;
    (*vInfoList).get_faces_for_edge(idx0, idx1, &adjFaceIds);

    if (adjFaceIds.size() > 1) return false;

    mve::TriangleMesh::VertexList const& pVerts = patchMesh->get_vertices();
    mve::TriangleMesh::VertexList const& mVerts = this->mesh->get_vertices();

    float d0, d1;
    d0 = (pVerts[idx0] - mVerts[this->centerVertId]).norm();
    d1 = (pVerts[idx1] - mVerts[this->centerVertId]).norm();
    return EPSILON_EQUAL(d0, radius, 1e-7) && EPSILON_EQUAL(d1, radius, 1e-7);
}

bool RegionGrower::performFlap(
    Edge const& e, mve::TriangleMesh::Ptr patchMesh, float const& radius)
{
    mve::TriangleMesh::VertexList& pVerts = patchMesh->get_vertices();
    mve::TriangleMesh::FaceList& pFaces = patchMesh->get_faces();
    mve::TriangleMesh::VertexList const& mVerts = this->mesh->get_vertices();

    math::Vec3f v0, v1, eMiddle;
    v0 = pVerts[e.a];
    v1 = pVerts[e.b];
    math::Vec3f v = (v1 - v0);
    float eLength = v.norm();

    if (IS_NAN(eLength) || eLength < MIN_EDGE_LENGTH) return false;

    eMiddle = v0 + 0.5f * eLength * (v / eLength);

    math::Vec3f dirVector = (eMiddle - mVerts[this->centerVertId]);
    float dirVectorLength = dirVector.norm();
    float distToCenter = dirVectorLength;
    float diffToRadius = std::abs(distToCenter - radius);

    // check if triangle needs to be added
    if (diffToRadius > MIN_EDGE_TO_RADIUS_DIST) {
        // find point on radius
        math::Vec3f pointOnRadius, dir;
        dir = dirVector / dirVectorLength;
        pointOnRadius = mVerts[this->centerVertId] + radius * dir;

        unsigned int newVertId = pVerts.size();
        pVerts.push_back(pointOnRadius);

        pFaces.push_back(e.a);
        pFaces.push_back(newVertId);
        pFaces.push_back(e.b);

        return true;
    }
    return false;
}

NAMESPACE_CURVEST_END
