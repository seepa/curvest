/*
 * Copyright (C) 2016, Patrick Seemann
 *
 * This software may be modified and distributed under the terms
 * of the BSD 3-Clause license. See the LICENSE.txt file for details.
 */

#include <set>
#include <iostream>
#include <sstream>
#include <cmath>
#include <algorithm>
#include <limits>
#include <assert.h>

#include "defines.h"
#include "utils.h"
#include "curvest.h"
#include "regiongrower.h"
#include "sphere.h"
#include "KDTree.h"
#include "gnuplot_helpers.h"

#include "mve/mesh_io_ply.h"
#include "mve/mesh.h"
#include "mve/mesh_info.h"
#include "math/vector.h"
#include "util/timer.h"

NAMESPACE_CURVEST_BEGIN

void Estimator::compute()
{
    std::cout << "Computing multi-scale curvature field:" << std::endl;

    this->initialize();

    util::WallTimer timer;
    std::size_t tComputeInitialRadius, tRefineRadii, tSmoothFinalRadii,
        tComputeFinalMC, total;

    std::vector<std::size_t> skippedVerts;
    if (!this->config.useGivenConfidences) {
        // compute initial radius per vertex
        std::cout << "Computing initial radii..." << std::endl;
        this->computeInitialRadius();
        tComputeInitialRadius = timer.get_elapsed();
        timer.reset();

        // compute mean curvature at multiple radii and decide on final radius
        // per vertex
        std::cout << "Perform radius refinement..." << std::endl;
        this->refineRadii(skippedVerts);
        tRefineRadii = timer.get_elapsed();
        timer.reset();

// smooth final radii
#if SMOOTH_FINAL_RADIUS
        if (!this->config.debugVertex) {
            std::cout << "Smoothing final radii..." << std::endl;
            this->smoothFinalRadii();
            tSmoothFinalRadii = timer.get_elapsed();
            timer.reset();
        }
#endif
    } else {
        std::cout << "Using given confidences as final radii." << std::endl;
    }

    // compute mean curvature based on final radii
    std::cout << "Computing mean curvatures for final radii..." << std::endl;
    this->computeFinalMeanCurvature();
    tComputeFinalMC = timer.get_elapsed();

#if CLAMP_FINAL_CURVATURE
        // clamp mean curvature using percentiles
        if (!this->config.debugVertex)
            this->clampVertexValues();
#endif

    total = tComputeInitialRadius + tComputeInitialRadius + tRefineRadii +
            tSmoothFinalRadii + tComputeFinalMC;

    std::cout << "Statistics: " << std::endl;
    std::cout << "   Time initial radius: "
              << utils::getPrettyTimeString(tComputeInitialRadius) << std::endl;
    std::cout << "   Time refine radii: " << utils::getPrettyTimeString(
                                                 tRefineRadii) << std::endl;
    std::cout << "   Time smooth final radii: "
              << utils::getPrettyTimeString(tSmoothFinalRadii) << std::endl;
    std::cout << "   Time compute final curvature: "
              << utils::getPrettyTimeString(tComputeFinalMC) << std::endl;
    std::cout << "Total: " << utils::getPrettyTimeString(total) << std::endl;

    if (!skippedVerts.empty()) {
        std::cout << "Warning: skipped vertices: ";
        for (std::size_t i = 0; i < skippedVerts.size(); ++i)
            std::cout << skippedVerts[i] << " ";
        std::cout << std::endl;
    }
}

void Estimator::computeFixedRadius(float const radius)
{
    std::cout << "Computing single-scale curvature field with radius " << radius
              << "..." << std::endl;
    this->initialize();
    util::WallTimer timer;
    mve::TriangleMesh::ConfidenceList& vertexConfidences =
        this->mesh->get_vertex_confidences();
    for (std::size_t i = 0; i < vertexConfidences.size(); ++i)
        vertexConfidences[i] = radius;
    this->computeFinalMeanCurvature();
#if CLAMP_FINAL_CURVATURE
    // clamp mean curvature using percentiles
    if (!this->config.debugVertex)
        this->clampVertexValues();
#endif
    std::cout << "Statistics: " << std::endl;
    std::cout << "    Time curvature computation: "
              << utils::getPrettyTimeString(timer.get_elapsed()) << std::endl;
}

void Estimator::initialize()
{
    // resize vertex values and confidences
    mve::TriangleMesh::VertexList const& mVerts = this->mesh->get_vertices();
    this->mesh->get_vertex_confidences().resize(mVerts.size());
    this->mesh->get_vertex_values().resize(mVerts.size());

    // create a sphere
    this->initSphere();
}

void Estimator::computeInitialRadius()
{
    mve::TriangleMesh::VertexList const& mVerts = this->mesh->get_vertices();
    mve::TriangleMesh::ConfidenceList& vertexConfidences =
        this->mesh->get_vertex_confidences();

#pragma omp parallel for schedule(dynamic)
    for (std::size_t i = 0; i < mVerts.size(); ++i) {
        mve::MeshInfo::VertexInfo const& vinf = this->meshInfo[i];
        if (vinf.verts.size() <= 2) {
            vertexConfidences[i] = -1.0f;
            continue;
        }
        vertexConfidences[i] = 0.0f;
        for (std::size_t k = 0; k < vinf.verts.size(); ++k)
            vertexConfidences[i] += (mVerts[i] - mVerts[vinf.verts[k]]).norm();
        vertexConfidences[i] /= vinf.verts.size();
    }

    utils::smoothMeshProperty(
        mesh, vertexConfidences, this->config.radSmoothIters, 1.0f);
}

void Estimator::refineRadii(std::vector<std::size_t>& skippedVerts)
{
    mve::TriangleMesh::VertexList const& mVerts = this->mesh->get_vertices();
    mve::TriangleMesh::ConfidenceList& vertexConfidences =
        mesh->get_vertex_confidences();
    mve::TriangleMesh::ValueList& vertexValues = mesh->get_vertex_values();

    // DEBUGGING
    std::vector<std::size_t> stats(7, 0);

    util::WallTimer wTimer;
#pragma omp parallel for schedule(dynamic)
    for (std::size_t i = 0; i < mVerts.size(); ++i) {
        if (this->config.debugVertex) {
            if (this->config.vertIds.find(i) == this->config.vertIds.end())
                continue;

            std::cout << std::endl;
            std::cout << "//-------------------------" << std::endl;
            std::cout << "//     Vertex " << i << std::endl;
            std::cout << "//-------------------------" << std::endl;
        }

        if (i % 200 == 0) {
            std::size_t elapsed = wTimer.get_elapsed();
            float timePerVert = elapsed / (float)i;
            std::size_t remainingVerts = mVerts.size() - i;
            std::size_t remainingTime = remainingVerts * timePerVert;
            std::size_t percent = (i / (float)mVerts.size()) * 100.0f;
            std::cout << i << " of " << mVerts.size() << " (" << percent
                      << "%) ETA " << utils::getPrettyTimeString(remainingTime)
                      << std::endl;
        }

        // skip unreferenced vertices / faces
        if (this->meshInfo[i].verts.size() <= 2) {
            vertexConfidences[i] = 0.0f;
            continue;
        }

        //
        // Compute mean curvature at multiple radii
        //
        float initialRadius, radiusFactor;
        initialRadius = vertexConfidences[i];
        radiusFactor = this->config.radFactorInitial;

        if (this->config.debugVertex)
            std::cout << "Intial radius = " << initialRadius << std::endl;

        utils::Result origH;
        origH.vertexId = i;
        origH.name = "Mean curvature";
        origH.useCustomLabels = true;
        origH.xLabel = "radius";
        origH.yLabel = "mean curvature";

        utils::Result normH;
        normH.vertexId = i;
        normH.name = "Fitting result";
        normH.useCustomLabels = true;
        normH.xLabel = "radius";
        normH.yLabel = "normalized mean curvature";

        // precompute sample radii
        std::vector<float> sampleRadii;
        while (radiusFactor <= this->config.radFactorMax) {
            sampleRadii.push_back(initialRadius * radiusFactor);
            radiusFactor *= this->config.radFactorInc;
        }

        // preselect three radii at 20%, 50% and 100% and compute mean curvature
        // at those radii. If the curvature is planar at all three measures, we
        // can break early for this vertex and save some computations.
        std::size_t numRadii = sampleRadii.size();
        std::set<std::size_t> indices;
        indices.insert(numRadii * 0.2f);
        indices.insert(numRadii * 0.5f);
        indices.insert(numRadii - 1);

        if (indices.size() != 3)
            std::cout << "Warning: problem during presel. radii!" << std::endl;

        // collect mean curvature at preselected radii
        bool res = true;
        for (auto it = indices.begin(); (it != indices.end()) && res; ++it) {
            float const& radius = sampleRadii[*it];
            res = computeDataForSingleRadius(i, radius, origH, normH);
        }
        if (!res) {
            std::cout << "Skipping vertex: Cannot compute H for " << i << "!"
                      << std::endl;
            vertexValues[i] = 0.0f;
#pragma omp critical
            skippedVerts.push_back(i);
            continue;
        }

        // check if the surface is planar
        if (checkIfPlanarExact(normH)) {
            if (this->config.debugVertex)
                std::cout << "Planer exact! break early..." << std::endl;

            // assign largest radius
            vertexConfidences[i] = normH.data[normH.data.size() - 1][0];

            normH.name = "Classification: Planar";

            // DEBUGGING: write patch at final radius
            savePatchToDisk(i, vertexConfidences[i]);
            if (this->config.saveGraphs)
                utils::saveGraphToDisk(
                    origH, normH, vertexConfidences[i], this->config.outfile,
                    "planar_exact");

            stats[FittingResult::PLANAR_EARLY]++;
            continue;
        }

        // compute mean curvature for the remaining radii
        res = true;
        for (std::size_t k = 0; (k < sampleRadii.size()) && res; ++k) {
            // skip already computed radii
            if (indices.find(k) != indices.end()) continue;

            float const& radius = sampleRadii[k];
            res = computeDataForSingleRadius(i, radius, origH, normH);
        }
        if (!res) {
            std::cout << "Skipping vertex: Cannot compute H for " << i << "!"
                      << std::endl;
            vertexValues[i] = 0.0f;
#pragma omp critical
            skippedVerts.push_back(i);
            continue;
        }

        // sort data from smallest to largest radius
        std::sort(normH.data.begin(), normH.data.end(),
            utils::Result_less_than_key());

        // check planar again
        float planarRadius;
        if (checkIfPlanarApprox(normH, planarRadius)) {
            if (this->config.debugVertex)
                std::cout << "Planar approx!" << std::endl;
            vertexConfidences[i] = planarRadius;

            normH.name = "Class.: Planar";
            // DEBUGGING:
            savePatchToDisk(i, vertexConfidences[i]);
            if (this->config.saveGraphs)
                utils::saveGraphToDisk(
                    origH, normH, vertexConfidences[i],
                    this->config.outfile, "planar_appr");

            stats[FittingResult::PLANAR]++;
            continue;
        }

        // find optimal fit
        math::Vec4f coeffs;
        float fitError;
        findOptimalFit(i, origH, normH, coeffs, fitError);

        // select final radius based on the fit
        float radius;
        FittingResult fitRes =
            chooseFinalRadius(i, normH, coeffs, fitError, radius);
        stats[fitRes]++;
        vertexConfidences[i] = radius;

        // DEBUGGING: write patch at final radius
        savePatchToDisk(i, vertexConfidences[i]);

        // DEBUGGING: save graphs as gnuplot svg's
        if (this->config.saveGraphs) {
            normH.name = "Class.: Non-planar";
            utils::saveGraphToDisk(
                origH, normH, vertexConfidences[i], fitError, coeffs,
                this->config.outfile, "final");
        }
    }

    // DEBUGGING: print stats for each case
    std::size_t total = mVerts.size();
    const auto percent = [total](const std::size_t x)
                             -> float { return (x / (float)total) * 100.0f; };
    std::cout << "Fitting stats:" << std::endl;
    std::cout << "  NAN " << percent(stats[FittingResult::NO_EXTREMA_NAN])
              << "%" << std::endl;
    std::cout << "  One " << percent(stats[FittingResult::ONE_EXTREMA]) << "%"
              << std::endl;
    std::cout << "  Two " << percent(stats[FittingResult::TWO_EXTREMA]) << "%"
              << std::endl;
    std::cout << "  No " << percent(stats[FittingResult::NO_EXTREMA]) << "%"
              << std::endl;
    std::cout << "  Planar " << percent(stats[FittingResult::PLANAR]) << "%"
              << std::endl;
    std::cout << "  Planar early "
              << percent(stats[FittingResult::PLANAR_EARLY]) << "%"
              << std::endl;
    std::cout << "  else " << percent(stats[FittingResult::ALL_ELSE]) << "%"
              << std::endl;
}

void Estimator::findOptimalFit(std::size_t vertId, utils::Result& origH,
    utils::Result& normH, math::Vec4f& coeffs, float& fitError)
{
    // fit cubic polynomial
    float bestError;
    math::Vec4f bestCoeffs;
    utils::fitCubicPolynomial(normH.data, bestCoeffs, bestError);

    // DEBUGGING: save first fit to disk
    if (this->config.saveGraphs) {
        float radius;
        chooseFinalRadius(vertId, normH, coeffs, bestError, radius);
        utils::saveGraphToDisk(origH, normH, radius, bestError, bestCoeffs,
            this->config.outfile, "");
    }

    // TODO: classify fit to find planar region ??

    if (this->config.debugVertex)
        std::cout << "First fit error: " << bestError << std::endl;

    while (bestError > 0.02 && normH.data.size() > 5) {
        // remove largest radius
        std::size_t idxLast = normH.data.size() - 1;
        normH.data.erase(normH.data.begin() + idxLast);

        // refit
        float newError;
        math::Vec4f newCoeffs;
        utils::fitCubicPolynomial(normH.data, newCoeffs, newError);

        if (this->config.debugVertex)
            std::cout << "New error: " << newError << std::endl;

        // DEBUGGING: save first fit to disk
        if (this->config.saveGraphs) {
            float radius;
            chooseFinalRadius(vertId, normH, newCoeffs, newError, radius);
            utils::saveGraphToDisk(origH, normH, radius, newError, newCoeffs,
                this->config.outfile, "");
        }

        if (newError < bestError) {
            bestError = newError;
            bestCoeffs = newCoeffs;
        } else
            break;

        // TODO: check if we are planar now ??
    }

    coeffs = bestCoeffs;
    fitError = bestError;
}

bool Estimator::checkIfPlanarExact(utils::Result const& normH)
{
    for (std::size_t i = 0; i < normH.data.size(); ++i) {
        if (std::abs(normH.data[i][1]) > this->config.planarThres) return false;
    }
    return true;
}

bool Estimator::checkIfPlanarApprox(utils::Result const& normH, float& radius)
{
    float mean = 0.0f;
    for (std::size_t l = 0; l < normH.data.size(); ++l)
        mean += std::abs(normH.data[l][1]);
    mean /= normH.data.size();

    float meanRatio = mean / this->config.planarThres;
    if (mean <= this->config.planarThres && meanRatio < 0.9f) {
        // planar: use last radius below the threshold
        std::size_t num = normH.data.size();
        radius = normH.data[num - 1][0];
        if (std::abs(normH.data[num - 1][1]) > this->config.planarThres) {
            for (std::size_t i = num - 2; i >= 0; --i) {
                float const& normV = std::abs(normH.data[i][1]);
                if (normV <= this->config.planarThres) {
                    radius = normH.data[i][0];
                    break;
                }
            }
        }
        return true;
    }
    return false;
}

bool Estimator::computeDataForSingleRadius(unsigned int vertId,
    float const radius, utils::Result& origCurvatures, utils::Result& normH)
{
    // scale sphere with current radius
    mve::TriangleMesh::Ptr sphereCopy = this->getSphereCopy(radius);
    double intersectionVolume;
    bool res = this->computeIntersectionVolume(
        vertId, radius, sphereCopy, intersectionVolume);

    if (!res) {
        std::cout << "Warning: Volume computation failed!" << std::endl;
        return false;
    }

    // Calculate mean curvature H using the intersection volume
    double H = this->computeMeanCurvature(radius, intersectionVolume);

    // save values for current radius
    origCurvatures.data.push_back(math::Vec2d(radius, H));
    normH.data.push_back(math::Vec2f(radius, radius * H));
    return true;
}

bool Estimator::computeIntersectionVolume(unsigned int vertId, float radius,
    mve::TriangleMesh::Ptr currentSphere, double& volumeOut)
{
    mve::TriangleMesh::VertexList const& mVerts = this->mesh->get_vertices();

    RegionGrower regionGrower(
        vertId, this->mesh, &this->meshInfo, this->config);

    // get a surface patch of vertices inside the radius
    mve::TriangleMesh::Ptr patchMesh = mve::TriangleMesh::create();
    regionGrower.getSurfacePatch(radius, patchMesh);

    if (patchMesh->get_vertices().empty()) return false;

    // DEBUGGING: write out untranslated patch
    if (this->config.savePatchAll) {
        mve::TriangleMesh::ColorList& colors = patchMesh->get_vertex_colors();
        if (colors.empty()) {
            colors.resize(patchMesh->get_vertices().size(),
                math::Vec4f(1.0f, 1.0f, 1.0f, 1.0f));
        }
        mve::TriangleMesh::FaceList const& faces = patchMesh->get_faces();
        for (std::size_t k = 0; k < faces.size(); k = k + 3) {
            for (std::size_t l = 0; l < 3; ++l)
                colors[faces[k + l]] = math::Vec4f(1.0f, 0.0f, 0.0f, 1.0f);
        }
        // save mesh
        mve::geom::SavePLYOptions options;
        options.write_vertex_colors = true;
        options.write_face_normals = true;
        std::stringstream ss;
        ss << this->config.outfile << "_" << vertId << "_patch_r" << radius
           << ".ply";
        mve::geom::save_ply_mesh(patchMesh, ss.str(), options);
    }

    // translate patch to origin
    math::Vec3f const& centerVert = mVerts[vertId];
    math::Vec3f t = -centerVert;
    mve::TriangleMesh::VertexList& patchVerts = patchMesh->get_vertices();
    for (std::size_t i = 0; i < patchVerts.size(); ++i) patchVerts[i] += t;

    // calculate face normals of patch
    patchMesh->recalc_normals(true, true);

    // get list of vertices inside the intersection
    std::vector<std::size_t> faceIdsIntersection;
    this->getOptimizedIntersection(
        patchMesh, currentSphere, faceIdsIntersection);

    // DEBUGGING: write out volumes
    if (this->config.saveVolumeAll) {
        // for debugging, copy sphere and translate to centerVert
        mve::TriangleMesh::Ptr sphereCopy = currentSphere->duplicate();
        mve::TriangleMesh::VertexList& sVerts = sphereCopy->get_vertices();
        for (std::size_t i = 0; i < sVerts.size(); ++i) sVerts[i] -= t;

        // colorize vertices inside the intersection blue, all others white
        mve::TriangleMesh::ColorList& colors = sphereCopy->get_vertex_colors();
        math::Vec4f outColor(1.0f, 1.0f, 1.0f, 1.0f);
        if (colors.empty())
            colors.resize(sphereCopy->get_vertices().size(), outColor);
        math::Vec4f inColor(0.6f, 0.6705882352f, 0.0f, 1.0f);
        mve::TriangleMesh::FaceList const& faces = sphereCopy->get_faces();
        for (std::size_t k = 0; k < faceIdsIntersection.size(); ++k) {
            std::size_t faceId = faceIdsIntersection[k];
            for (std::size_t l = 0; l < 3; ++l)
                colors[faces[faceId + l]] = inColor;
        }

        // save mesh
        mve::geom::SavePLYOptions options;
        options.write_vertex_colors = true;
        std::stringstream ss;
        ss << this->config.outfile << "_" << vertId << "_sphere_r" << radius
           << ".ply";
        std::cout << "Writing " << ss.str() << std::endl;
        mve::geom::save_ply_mesh(sphereCopy, ss.str(), options);
    }

    // compute volume magnitude
    {
        // add terms for sphere intersection faces
        mve::TriangleMesh::VertexList const& sphereVerts =
            currentSphere->get_vertices();
        mve::TriangleMesh::FaceList const& sphereFaces =
            currentSphere->get_faces();
        double sum = 0.0;
        math::Vec3f a_k, b_k, c_k, n;
        for (std::size_t k = 0; k < faceIdsIntersection.size(); ++k) {
            std::size_t id = faceIdsIntersection[k];
            a_k = sphereVerts[sphereFaces[id]];
            b_k = sphereVerts[sphereFaces[id + 1]];
            c_k = sphereVerts[sphereFaces[id + 2]];
            n = (b_k - a_k).cross(c_k - a_k);

            // multiply by correction factor to minimize error
            sum += this->sphereVolumeCF * a_k.dot(n);
        }

        // add terms for patch faces
        mve::TriangleMesh::FaceList const& patchFaces = patchMesh->get_faces();
        for (std::size_t k = 0; k < patchFaces.size(); k += 3) {
            a_k = patchVerts[patchFaces[k]];
            c_k = patchVerts[patchFaces[k + 1]];
            b_k = patchVerts[patchFaces[k + 2]];
            n = (b_k - a_k).cross(c_k - a_k);
            sum += a_k.dot(n);
        }
        sum *= 1.0 / 6.0;
        volumeOut = sum;
    }
    return true;
}

void Estimator::getOptimizedIntersection(mve::TriangleMesh::ConstPtr patchMesh,
    mve::TriangleMesh::Ptr currentSphere, std::vector<std::size_t>& faceIdsOut)
{
    // kdtree of vertices in patch
    KDTree<3> kdTree(patchMesh->get_vertices());

    std::vector<std::size_t> sphereVertsInside;
    this->getIntersectionSphere(
        patchMesh, &kdTree, currentSphere, sphereVertsInside);

    // find faces on the intersection border
    std::vector<std::size_t> borderFaces;
    mve::TriangleMesh::FaceList& sphereFaces = currentSphere->get_faces();
    for (std::size_t i = 0; i < sphereFaces.size(); i += 3) {
        std::size_t const& idx0 = sphereFaces[i];
        std::size_t const& idx1 = sphereFaces[i + 1];
        std::size_t const& idx2 = sphereFaces[i + 2];

        int c = 0;
        c += (int)(sphereVertsInside[idx0] == 1);
        c += (int)(sphereVertsInside[idx1] == 1);
        c += (int)(sphereVertsInside[idx2] == 1);

        if (c >= 1 && c <= 2) borderFaces.push_back(i);
    }

    // refine border by iteratively subdiving border faces
    for (std::size_t iter = 0; iter < INTERSECTION_SUBDIV_ITERS; ++iter) {
        // divide each border face into two new faces
        std::vector<std::size_t> todoFaceIds;
        todoFaceIds.reserve(borderFaces.size() * 2);
        for (std::size_t i = 0; i < borderFaces.size(); i++) {
            // Split face into two. This will create one new vertex
            splitFace(currentSphere, todoFaceIds, borderFaces[i]);

            // add one position to inclusion vector for the new vertex
            sphereVertsInside.push_back(0);
        }

        // Check for each new face which vertices are within the intersection
        mve::TriangleMesh::VertexList& sVerts = currentSphere->get_vertices();
        // keep track which ones were checked already
        std::vector<std::size_t> checkedVerts(sVerts.size(), 0);
        for (std::size_t i = 0; i < todoFaceIds.size(); ++i) {
            std::size_t const& fId = todoFaceIds[i];

            for (std::size_t j = 0; j < 3; ++j) {
                std::size_t vIdx = sphereFaces[fId + j];
                if (!checkedVerts[vIdx]) {
                    sphereVertsInside[vIdx] = checkInsideIntersection(
                        patchMesh, &kdTree, sVerts[vIdx]);
                    checkedVerts[vIdx] = 1;
                }
            }
        }

        // Get border faces for next iteration. Only the faces which were
        // split have to be checked
        borderFaces.clear();
        for (std::size_t i = 0; i < todoFaceIds.size(); ++i) {
            std::size_t const& fId = todoFaceIds[i];

            std::size_t const& idx0 = sphereFaces[fId];
            std::size_t const& idx1 = sphereFaces[fId + 1];
            std::size_t const& idx2 = sphereFaces[fId + 2];

            int c = 0;
            c += (int)(sphereVertsInside[idx0] == 1);
            c += (int)(sphereVertsInside[idx1] == 1);
            c += (int)(sphereVertsInside[idx2] == 1);

            if (c >= 1 && c <= 2) borderFaces.push_back(fId);
        }
    }

    // write faceIdsOut
    for (std::size_t i = 0; i < sphereFaces.size(); i += 3) {
        std::size_t const& idx0 = sphereFaces[i];
        std::size_t const& idx1 = sphereFaces[i + 1];
        std::size_t const& idx2 = sphereFaces[i + 2];

        // A face is considered to be within the intersection, iff all three
        // vertices are within. This should not be too restrictive, since
        // we subdivided all faces which were on the border of intersection.
        if ((sphereVertsInside[idx0] && sphereVertsInside[idx1] &&
                sphereVertsInside[idx2]))
            faceIdsOut.push_back(i);
    }
}

void Estimator::splitFace(mve::TriangleMesh::Ptr mesh,
    std::vector<std::size_t>& facesTodo, std::size_t fId)
{
    mve::TriangleMesh::VertexList& sVerts = mesh->get_vertices();
    mve::TriangleMesh::FaceList& sFaces = mesh->get_faces();

    std::size_t aIdx = sFaces[fId];
    std::size_t bIdx = sFaces[fId + 1];
    std::size_t cIdx = sFaces[fId + 2];

    // get longest edge
    math::Vec3f const& a = sVerts[aIdx];
    math::Vec3f const& b = sVerts[bIdx];
    math::Vec3f const& c = sVerts[cIdx];

    float abl = (b - a).norm();
    float bcl = (c - b).norm();
    float cal = (a - c).norm();

    // determine indices to create two new triangles
    // idx0 will always be oposite the longest side
    std::size_t idx0, idx1, idx2;
    if (abl >= bcl) {
        if (abl >= cal) {
            // abl is longest side
            idx0 = cIdx;
            idx1 = aIdx;
            idx2 = bIdx;
        } else {
            // cal is longest
            idx0 = bIdx;
            idx1 = cIdx;
            idx2 = aIdx;
        }
    } else {
        if (bcl >= cal) {
            // bcl is longest
            idx0 = aIdx;
            idx1 = bIdx;
            idx2 = cIdx;
        } else {
            // cal is longest
            idx0 = bIdx;
            idx1 = cIdx;
            idx2 = aIdx;
        }
    }

    // get middle of edge to cut
    math::Vec3f const& v1 = sVerts[idx1];
    math::Vec3f const& v2 = sVerts[idx2];
    math::Vec3f middle = v1 + 0.5f * (v2 - v1);

    // Save the middle position as a new vertex
    std::size_t newIdx = sVerts.size();
    sVerts.push_back(middle);

    // Create two new faces without shifting face indices: One new face will
    // replace the old face.

    // write first one at position of old face
    sFaces[fId] = idx0;
    sFaces[fId + 1] = idx1;
    sFaces[fId + 2] = newIdx;
    facesTodo.push_back(fId);

    // append second face to the end
    facesTodo.push_back(sFaces.size());
    sFaces.push_back(idx0);
    sFaces.push_back(newIdx);
    sFaces.push_back(idx2);
}

bool Estimator::checkInsideIntersection(mve::TriangleMesh::ConstPtr mesh,
    KDTree<3>* kdTree, math::Vec3f const& vertex)
{
    mve::TriangleMesh::NormalList const& vertNormals =
        mesh->get_vertex_normals();
    mve::TriangleMesh::VertexList const& verts = mesh->get_vertices();

    // get mesh vertex which is closest to given vertex
    std::pair<std::size_t, float> res = kdTree->find_nn(vertex);
    math::Vec3f const& vertNormal = vertNormals[res.first];
    math::Vec3f const& nearestPoint = verts[res.first];
    math::Vec3f const a = nearestPoint - vertex;

    // face is behind surface
    if (a.dot(vertNormal) >= 0.0f) return true;

    return false;
}

void Estimator::getIntersectionSphere(mve::TriangleMesh::ConstPtr mesh,
    KDTree<3>* kdTree, mve::TriangleMesh::Ptr sphere,
    std::vector<std::size_t>& vertsInside)
{
    mve::TriangleMesh::VertexList const& sphereVerts = sphere->get_vertices();
    vertsInside.assign(sphereVerts.size(), 0);
    for (std::size_t i = 0; i < sphereVerts.size(); ++i) {
        if (checkInsideIntersection(mesh, kdTree, sphereVerts[i]))
            vertsInside[i] = 1;
    }
}

Estimator::FittingResult Estimator::chooseFinalRadius(std::size_t const vertId,
    utils::Result& normH, math::Vec4f coeffsCubic, float const& fitError,
    float& finalR)
{
    Estimator::FittingResult returnValue;

    // compute some radius values (x-axis) for further analysis
    float smallestR, largestR, middleR, closeToLargestR, closeToSmallestR;
    smallestR = normH.data[0][0];
    largestR = normH.data[normH.data.size() - 1][0];
    middleR = smallestR + (largestR - smallestR) / 2.0f;
    closeToLargestR = smallestR + 9.0f * ((largestR - smallestR) / 10.0f);
    closeToSmallestR = smallestR + 1.0f * ((largestR - smallestR) / 10.0f);

    if (this->config.debugVertex) {
        std::cout << "Id: " << vertId << "smallestR: " << smallestR
                  << " largest: " << largestR
                  << " closeToSmallest: " << closeToSmallestR
                  << " closeToLargest: " << closeToLargestR
                  << " middleR: " << middleR << std::endl;
    }

    // find extrema and inflection point of cubic function
    float x1, x2, inflPoint;
    utils::getInfoCubicPoly(coeffsCubic, x1, x2, inflPoint, this->config.debugVertex);
    if (this->config.debugVertex)
        std::cout << "x1 " << x1 << " x2 " << x2 << std::endl;

    bool inflWithin =
        (inflPoint >= closeToSmallestR && inflPoint <= closeToLargestR);

    if (IS_NAN(x1) || IS_NAN(x2)) {
        float point = middleR;
        if (inflWithin && fitError <= 0.01) point = inflPoint;
        finalR = smallestR + this->config.edgeSmooth * (point - smallestR);
        returnValue = FittingResult::NO_EXTREMA_NAN;
    } else {
        bool x1Within = (x1 >= smallestR && x1 <= largestR);
        bool x2Within = (x2 >= smallestR && x2 <= largestR);

        if ((x1Within && !x2Within) || (x2Within && !x1Within)) {
            if (this->config.debugVertex)
                std::cout << "only x1 XOR x2 inside" << std::endl;

            float x = x2Within ? x2 : x1;
            finalR = this->config.edgeSmooth * x +
                     (1.0f - this->config.edgeSmooth) * smallestR;

            returnValue = FittingResult::ONE_EXTREMA;
        } else if (!x1Within && !x2Within) {
            if (this->config.debugVertex)
                std::cout << "x1,x2 outside" << std::endl;

            float point = middleR;
            if (inflWithin && fitError <= 0.01) point = inflPoint;
            finalR = smallestR + this->config.edgeSmooth * (point - smallestR);
            returnValue = FittingResult::NO_EXTREMA;
        } else if (x1Within && x2Within && inflWithin) {
            if (this->config.debugVertex)
                std::cout << "all within... " << std::endl;
            finalR = smallestR + this->config.edgeSmooth * (x1 - smallestR);
            returnValue = FittingResult::TWO_EXTREMA;
        } else {
            if (this->config.debugVertex)
                std::cout << "interpolate, tend to first extrema; "
                          << std::endl;

            // tend to first extrema
            bool inflUseful = (inflWithin && inflPoint > closeToSmallestR &&
                               inflPoint < closeToLargestR);

            // select either infl point or middle radius
            float point = inflUseful ? inflPoint : middleR;
            float alpha = 1.0f - this->config.edgeSmooth;
            float interpol1, interpol2;
            if (x1 <= smallestR) {
                interpol1 = smallestR;
                interpol2 = point;
            } else {
                interpol1 = smallestR;
                interpol2 = x1;
            }

            if (this->config.debugVertex) {
                std::cout << "Interpolate btw.: " << interpol1 << " and "
                          << interpol2 << " alpha " << alpha << std::endl;
            }

            finalR = alpha * interpol1 + (1.0f - alpha) * interpol2;

            returnValue = FittingResult::ALL_ELSE;
        }
    }
    return returnValue;
}

void Estimator::savePatchToDisk(std::size_t const& vertId, float const& radius)
{
    // DEBUGGING: write patch at final radius
    if (!this->config.savePatchFinal) return;

    RegionGrower regionGrower(
        vertId, this->mesh, &this->meshInfo, this->config);
    mve::TriangleMesh::Ptr patchMesh = mve::TriangleMesh::create();
    regionGrower.getSurfacePatch(radius, patchMesh);

    mve::TriangleMesh::ColorList& colors = patchMesh->get_vertex_colors();
    if (colors.empty()) {
        colors.resize(patchMesh->get_vertices().size(),
            math::Vec4f(1.0f, 1.0f, 1.0f, 1.0f));
    }
    mve::TriangleMesh::FaceList const& faces = patchMesh->get_faces();
    for (std::size_t k = 0; k < faces.size(); k = k + 3) {
        for (std::size_t l = 0; l < 3; ++l)
            colors[faces[k + l]] = math::Vec4f(1.0f, 0.0f, 0.0f, 1.0f);
    }
    // save mesh
    mve::geom::SavePLYOptions options;
    options.write_vertex_colors = true;
    options.write_face_normals = true;
    std::stringstream ss;
    ss << this->config.outfile << "_" << vertId << "_patch_r" << radius
       << ".ply";
    mve::geom::save_ply_mesh(patchMesh, ss.str(), options);
}

void Estimator::computeFinalMeanCurvature()
{
    mve::TriangleMesh::VertexList const& mVerts = this->mesh->get_vertices();
    mve::TriangleMesh::ConfidenceList const& vertexConfidences =
        this->mesh->get_vertex_confidences();
    mve::TriangleMesh::ValueList& vertexValues =
        this->mesh->get_vertex_values();

    util::WallTimer finalMCTimer;
#pragma omp parallel for schedule(dynamic)
    for (std::size_t i = 0; i < mVerts.size(); ++i) {
        // debugging
        if (this->config.debugVertex) {
            if (this->config.vertIds.find(i) == this->config.vertIds.end())
                continue;
        }

        if (i % 200 == 0) {
            std::size_t elapsed = finalMCTimer.get_elapsed();
            float timePerVert = elapsed / (float)i;
            std::size_t remainingVerts = mVerts.size() - i;
            std::size_t remainingTime = remainingVerts * timePerVert;
            std::size_t percent = (i / (float)mVerts.size()) * 100.0f;
            std::cout << i << " of " << mVerts.size() << " (" << percent
                      << "%) ETA " << utils::getPrettyTimeString(remainingTime)
                      << std::endl;
        }

        // skip unreferenced vertices / faces
        if (this->meshInfo[i].verts.size() <= 2) {
            vertexValues[i] = 0.0f;
            continue;
        }

        float const& finalR = vertexConfidences[i];
        // scale sphere with current radius
        mve::TriangleMesh::Ptr sphereCopy = getSphereCopy(finalR);
        double intersectionVolume;

        bool res = this->computeIntersectionVolume(
            i, finalR, sphereCopy, intersectionVolume);

        if (!res) {
            vertexValues[i] = 0.0f;
            continue;
        }

        double H = this->computeMeanCurvature(finalR, intersectionVolume);
        vertexValues[i] = H;

        if (this->config.debugVertex)
            std::cout << "Final R: " << finalR << " H: " << H << std::endl;
    }
}

double Estimator::computeMeanCurvature(float radius, double intersectionVolume)
{
    double halfVolume = ((2.0 * MATH_PI * MATH_POW3(radius)) / 3.0);
    double H = (4.0 / (MATH_POW4(radius) * MATH_PI)) *
               (halfVolume - intersectionVolume);
    return H;
}

double Estimator::getSphereVolumeCorrectionFactor(
    mve::TriangleMesh::ConstPtr sphere)
{
    mve::MeshBase::VertexList const& vertices = sphere->get_vertices();
    mve::TriangleMesh::FaceList const& faces = sphere->get_faces();

    // assert sphere radius is 1.0
    assert(EPSILON_EQUAL(vertices[0].norm(), 1.0, 1e-6));

    // compute tessellated volume based on faces
    double sum = 0.0;
    for (std::size_t i = 0; i < faces.size(); i = i + 3) {
        math::Vec3f const& a_i = vertices[faces[i]];
        math::Vec3f const& b_i = vertices[faces[i + 1]];
        math::Vec3f const& c_i = vertices[faces[i + 2]];

        math::Vec3f n = (b_i - a_i).cross(c_i - a_i);
        sum += a_i.dot(n);
    }
    sum *= 1.0 / 6.0;

    // actual volume
    double volumeBall = (4.0 / 3.0) * MATH_PI;

    // correction factor
    double cf = volumeBall / sum;
    std::cout << "Sphere Volume Correction factor is: " << cf << std::endl;
    return cf;
}

mve::TriangleMesh::Ptr Estimator::getSphereCopy(float const& radius)
{
    mve::TriangleMesh::Ptr sphereCopy = this->sphere->duplicate();
    mve::TriangleMesh::VertexList& sphereVerts = sphereCopy->get_vertices();

    for (std::size_t i = 0; i < sphereVerts.size(); ++i)
        sphereVerts[i] = radius * this->dirSphereVerts[i];

    return sphereCopy;
}

void Estimator::smoothFinalRadii()
{
    mve::TriangleMesh::ConfidenceList& vertexConfidences =
        mesh->get_vertex_confidences();
    utils::smoothMeshProperty(
        this->mesh, vertexConfidences, FINAL_RADIUS_SMOOTHING_ITERS, 0.5f);
}

void Estimator::clampVertexValues()
{
    std::cout << "Clamping mean curvature to percentiles low: "
              << PERCENTILE_LOW << " high: " << PERCENTILE_HIGH << "..."
              << std::endl;
    mve::TriangleMesh::ConfidenceList& vertexValues =
        this->mesh->get_vertex_values();
    utils::percentileClamp(vertexValues, PERCENTILE_LOW, PERCENTILE_HIGH);
}

NAMESPACE_CURVEST_END
