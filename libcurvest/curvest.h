/*
 * Copyright (C) 2016, Patrick Seemann
 *
 * This software may be modified and distributed under the terms
 * of the BSD 3-Clause license. See the LICENSE.txt file for details.
 */

#pragma once

#include <set>

#include "defines.h"
#include "utils.h"
#include "KDTree.h"
#include "sphere.h"

#include "mve/mesh.h"
#include "mve/mesh_info.h"

NAMESPACE_CURVEST_BEGIN

// Multi-Scale Curvature Estimation: Given a triangle mesh, automatically
// and adaptively compute a radius at each vertex that is then used to compute
// the mean curvature using integral invariants by Yang et al [1].
//
// [1] Integral invariants (ball neighborhood) as described by
// "Yang, Yong-Liang, et al.: Robust principal curvatures on multiple
// scales. Symposium on Geometry Processing. 2006."
class Estimator
{
public:
    struct Config {
         //Number of Loop subdivision iterations for refining the sphere.
        std::size_t subdivIters;

        // The initial scale of each vertex is multiplied by this factor.
        // Increase this value to start off with a larger radius per vertex
        // (results in more smoothing of the final curvature).
        float radFactorInitial;

        // In each iteration, the current radius factor is multiplied by this.
        // This should be kept low to get enough samples for each vertex.
        float radFactorInc;

        // The radius factor is increased until this value is reached.
        float radFactorMax;

        // Number of iterations for smoothing the inital radius
        unsigned int radSmoothIters;

        // Controls how much edges are smoothed. Set this to 0.0 to perform
        // no smoothing. A value of 1.0 perform a lot of smoothing.
        // Note: A value of 0.0 can result in noisy results, so keeping this
        // value above 0.2 has proven to be useful.
        float edgeSmooth;

        // Controls up to when a vertex is considered to lay inside a planar
        // region. The range of this parameter is between 0.0 and inf.
        // A higher value results in smoothing of small geometric details.
        // In practice a value between 0.1 and 0.3 has shown best results.
        // Increasing this value might help when the mesh is really noisy or
        // when small geometry should be smoothed on purpose.
        float planarThres;

        // Compute the mean curvature at a single scale.
        bool usefixedRadius;

        // The radius to use for single-scale computation
        float fixedRadius;

        // Skip radius refinement step and use the confidences of the given
        // mesh to compute the final mean curvatures.
        bool useGivenConfidences;

        // Debugging options
        bool debugVertex;
        std::string outfile;
        std::set<unsigned int> vertIds;
        bool saveGraphs;
        bool saveVolumeAll;
        bool saveVolumeFinal;
        bool savePatchAll;
        bool savePatchFinal;
    };

    // DEBUGGIN: Different result types when fitting a polynomial
    enum FittingResult {
        NO_EXTREMA_NAN = 0,
        ONE_EXTREMA,
        TWO_EXTREMA,
        NO_EXTREMA,
        PLANAR,
        PLANAR_EARLY,
        ALL_ELSE
    };

    Estimator();

    // Sets the estimator config
    void setConfig(Config const& config);

    // Sets the mesh on which to perform computations
    void setMesh(mve::TriangleMesh::ConstPtr mesh);

    // Computes a multi-scale mean curvature field. The radius for each vertex
    // is saved in the output ply as its confidence and the mean curvature
    // computed at that radius is saved as its value.
    void compute();

    // Computes a single-scale mean curvature field given a radius. For every
    // vertex, the mean curvature is computed using the given radius. The
    // result is stored in the vertex's value field in the output ply.
    void computeFixedRadius(float const radius);

    // Returns a pointer to the resulting mesh
    mve::TriangleMesh::ConstPtr getResult();

    // Helper which computes the mean curvature for the given vertex at
    // multiple scales.
    bool computeDataForSingleRadius(unsigned int vertId, float const radius,
        utils::Result& origCurvatures, utils::Result& normalizedCurvatures);

    // Returns the number of vertices of the sphere
    std::size_t numSphereVerts();

private:
    // Resizes vertex values and vertex confidences vectors of the mesh
    void initialize();

    // Computes an initial radius for each vertex by calculating the average
    // edge lengths of its neighbors and then smoothing the radii over the
    // mesh.
    void computeInitialRadius();

    // Finds the final radius (aka the scale) for each vertex by computing its
    // mean curvature at multiple radii and then selecting the "best" one:
    // In planar regions the radius should be large where as in regions of
    // higher curvature, the radius should be small.
    // Note: The larger the radius the greater the smoothing of the mean
    // curvature at that point.
    //
    // The radius sampling can be controlled via the conifg.
    void refineRadii(std::vector<std::size_t>& skippedVerts);

    // Chooses the final radius given sampling data (radius and corresponding
    // normalized mean curvature) for a vertex. A cubic polynomial is used
    // to find interesting regions within the sampling.
    //
    // The final radius is written to the vertex confidence field of the vertex.
    FittingResult chooseFinalRadius(const std::size_t vertId,
        utils::Result& normH, math::Vec4f coeffsCubic, const float& fitError,
        float& finalR);

    // Calculates a correction factor that is used to minimize the error when
    // computing the intersection volume of a triangulated sphere and surface
    // patch.
    //
    // The factor is: f = (V_actual) / (V_approx)
    // where V_approx is the approximated volume of the sphere (centered at the
    // origin with radius 1.0) and V_actual is the actual volume of the sphere
    // computed analytically (using r = 1.0).
    double getSphereVolumeCorrectionFactor(mve::TriangleMesh::ConstPtr sphere);

    // Returns a copy of this->sphere where the vertices are scaled using
    // the given radius.
    mve::TriangleMesh::Ptr getSphereCopy(float const& radius);

    // Computes the intersection volume of a sphere and a surface patch (both
    // at the given radius) at the vertex with the given id.
    bool computeIntersectionVolume(unsigned int vertId, float radius,
        mve::TriangleMesh::Ptr currentSphere, double& volumeOut);

    // Returns a list of sphere face ids which form the intersection of the
    // sphere and the given patch of the mesh.
    void getIntersectionSphere(mve::TriangleMesh::ConstPtr mesh,
        KDTree<3>* kdTree, mve::TriangleMesh::Ptr sphere,
        std::vector<std::size_t>& vertsInside);

    // Computes the mean curvature.
    double computeMeanCurvature(float radius, double intersectionVolume);

    // Smoothes the final radii.
    void smoothFinalRadii();

    // Computes the final mean curvature for each vertex using its radius as
    // specified in its corresponding confidence field.
    void computeFinalMeanCurvature();

    // Checks whether all mean curvature measures are below the planar
    // threshold.
    bool checkIfPlanarExact(utils::Result const& normH);

    // Checks whether approximately all mean curvature measures are below the
    // planar threshold.
    bool checkIfPlanarApprox(utils::Result const& normH, float& radius);

    // Finds a cubic polynomial which represents the given result up to an
    // error.
    void findOptimalFit(std::size_t vertId, utils::Result& origH,
        utils::Result& normH, math::Vec4f& coeffs, float& fitError);

    // Helper to save a surface patch to disk.
    void savePatchToDisk(std::size_t const& vertId, float const& radius);

    // Calculates the intersection of the mesh with the given sphere. Returns
    // a list of face ids which are inside the intersection.
    void getOptimizedIntersection(mve::TriangleMesh::ConstPtr patchMesh,
        mve::TriangleMesh::Ptr currentSphere,
        std::vector<std::size_t>& faceIdsOut);

    // Splits a given face into two new faces. A new vertex is inserted on the
    // longest side of the face.
    //
    // One of the new faces is written into the position of the old face, so the
    // list of faces grows only by one.
    //
    // The face ids of the new faces will be added to the facesTodo vector.
    void splitFace(mve::TriangleMesh::Ptr mesh,
        std::vector<std::size_t>& facesTodo, std::size_t fId);

    // Checks for a single vertex if it is behind (an thus in the intersection)
    // or infront of the given mesh
    bool checkInsideIntersection(mve::TriangleMesh::ConstPtr mesh,
        KDTree<3>* kdTree, math::Vec3f const& vertex);

    // Initializes the sphere and precomputes direction vectors for each
    // vertex of the sphere.
    void initSphere();

    // Clamps the vertex values of the current mesh to PERCENTILE_LOW and
    // PERCENTILE_HIGH
    void clampVertexValues();

private:
    Config config;
    mve::TriangleMesh::Ptr mesh;
    mve::MeshInfo meshInfo;

    mve::TriangleMesh::ConstPtr sphere;
    std::vector<math::Vec3f> dirSphereVerts;
    float sphereVolumeCF;  // correction factor
};

/************************** Implementations **********************************/

inline Estimator::Estimator() : sphere(nullptr)
{
    this->config.radFactorInitial = DEFAULT_RADIUS_FACTOR_INITIAL;
    this->config.radFactorInc = DEFAULT_RADIUS_FACTOR_INC;
    this->config.radFactorMax = DEFAULT_RADIUS_FACTOR_MAX;
    this->config.edgeSmooth = DEFAULT_EDGE_SMOOTH;
    this->config.planarThres = DEFAULT_PLANAR_THRESHOLD;
    this->config.usefixedRadius = false;
    this->config.radSmoothIters = DEFAULT_INITIAL_RADIUS_SMOOTH_ITERS;
    this->config.debugVertex = false;
    this->config.saveGraphs = false;
    this->config.saveVolumeAll = false;
    this->config.saveVolumeFinal = false;
    this->config.savePatchAll = false;
    this->config.savePatchFinal = false;
}

inline void Estimator::setConfig(Config const& config)
{
    this->config = config;
}

inline void Estimator::setMesh(mve::TriangleMesh::ConstPtr mesh)
{
    this->mesh = mesh->duplicate();
    // make sure we have face and vertex normals
    this->mesh->recalc_normals(true, true);
    this->meshInfo.initialize(this->mesh);
}

inline void Estimator::initSphere()
{
    this->sphere = utils::getSphere(this->config.subdivIters);
    std::cout << "Created Sphere with a tessellation of "
              << this->sphere->get_vertices().size() << " vertices."
              << std::endl;

    // pre-compute sphere vertex direction vectors
    mve::TriangleMesh::VertexList const& sphereVerts =
        this->sphere->get_vertices();
    this->dirSphereVerts.resize(sphereVerts.size());
    for (std::size_t i = 0; i < sphereVerts.size(); ++i)
        this->dirSphereVerts[i] = sphereVerts[i] / sphereVerts[i].norm();

    // compute shpere volume correction factor
    this->sphereVolumeCF = getSphereVolumeCorrectionFactor(sphere);
}

inline std::size_t Estimator::numSphereVerts()
{
    if (this->sphere != nullptr) return this->sphere->get_vertices().size();
    return 0;
}

inline mve::TriangleMesh::ConstPtr Estimator::getResult()
{
    return this->mesh;
}
NAMESPACE_CURVEST_END
