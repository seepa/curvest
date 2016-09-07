/*
 * Copyright (C) 2016, Patrick Seemann
 *
 * This software may be modified and distributed under the terms
 * of the BSD 3-Clause license. See the LICENSE.txt file for details.
 */

#include <iostream>
#include <limits>
#include <sstream>

#include "util/arguments.h"
#include "util/tokenizer.h"

#include "mve/mesh_io_ply.h"
#include "mve/mesh.h"
#include "mve/mesh_info.h"

#include "libcurvest/curvest.h"
#include "libcurvest/utils.h"

int main(int argc, char** argv)
{
    // Setup argument parser.
    util::Arguments args;
    args.set_description("Multi-Scale Mean Curvature Estimation");
    args.set_exit_on_error(true);
    args.set_nonopt_maxnum(2);
    args.set_nonopt_minnum(2);

    // Help texts
    std::stringstream subdivIters, factorInitial, factorInc, factorMax,
        radSmooth, edgeSmooth, planarThres;
    subdivIters << "Number of subdivision iterations for refining the sphere ["
                << DEFAULT_SUBDIV_ITERS << "]";
    factorInitial << "Initial radius factor [" << DEFAULT_RADIUS_FACTOR_INITIAL
                  << "]";
    factorInc << "Radius factor increment [" << DEFAULT_RADIUS_FACTOR_INC
              << "]";
    factorMax << "Max radius factor [" << DEFAULT_RADIUS_FACTOR_MAX << "]";
    radSmooth << "Initial radius smoothing iterations ["
              << DEFAULT_INITIAL_RADIUS_SMOOTH_ITERS << "]";
    edgeSmooth << "Edge smoothing factor (0: no smooth, 1.0: full smooth) ["
               << DEFAULT_EDGE_SMOOTH << "]";
    planarThres << "Planar threshold [" << DEFAULT_PLANAR_THRESHOLD << "]";

    args.set_helptext_indent(22);
    args.set_usage(argv[0], "[ OPTIONS ] IN-PLY OUT");
    args.add_option('\0', "subdiv", true, subdivIters.str());
    args.add_option('\0', "factor-initial", true, factorInitial.str());
    args.add_option('\0', "factor-inc", true, factorInc.str());
    args.add_option('\0', "factor-max", true, factorMax.str());
    args.add_option('\0', "radius-smooth", true, radSmooth.str());
    args.add_option('e', "edge-smooth", true, edgeSmooth.str());
    args.add_option('p', "planar-thres", true, planarThres.str());
    args.add_option('d', "debug-verts", true,
        "(Debug) Comma separated list of vertices to debug");
    args.add_option(
        '\0', "save-graphs", false, "(Debug) Save a graph for each vertex");
    args.add_option(
        '\0', "save-volume-all", false, "(Debug) Save volume at each radius");
    args.add_option(
        '\0', "save-patch-all", false, "(Debug) Save patch at each radius");
    args.add_option(
        '\0', "save-patch-final", false, "(Debug) Save patch at final radius");
    args.add_option('\0', "fixed-radius", true,
        "Compute curvature field with a fixed ball radius.");
    args.add_option(
        '\0', "use-confidences", false, "Use mesh confidences as final radii.");

    args.parse(argc, argv);

    // Setup default options
    curvest::Estimator::Config conf;
    conf.subdivIters = DEFAULT_SUBDIV_ITERS;
    conf.radFactorInitial = DEFAULT_RADIUS_FACTOR_INITIAL;
    conf.radFactorInc = DEFAULT_RADIUS_FACTOR_INC;
    conf.radFactorMax = DEFAULT_RADIUS_FACTOR_MAX;
    conf.radSmoothIters = DEFAULT_INITIAL_RADIUS_SMOOTH_ITERS;
    conf.edgeSmooth = DEFAULT_EDGE_SMOOTH;
    conf.planarThres = DEFAULT_PLANAR_THRESHOLD;
    conf.debugVertex = false;
    conf.saveGraphs = false;
    conf.saveVolumeAll = false;
    conf.savePatchAll = false;
    conf.savePatchFinal = false;
    conf.usefixedRadius = false;
    conf.useGivenConfidences = false;

    // Parse arguments
    std::string inMeshPath = args.get_nth_nonopt(0);
    std::string outMeshPath = args.get_nth_nonopt(1);

    // DEBUGGING
    conf.outfile = outMeshPath;

    for (util::ArgResult const* arg = args.next_result(); arg != 0;
         arg = args.next_result()) {
        if (arg->opt == 0) continue;

        if (arg->opt->lopt == "subdiv")
            conf.subdivIters = arg->get_arg<std::size_t>();
        else if (arg->opt->lopt == "factor-initial")
            conf.radFactorInitial = arg->get_arg<float>();
        else if (arg->opt->lopt == "factor-inc")
            conf.radFactorInc = arg->get_arg<float>();
        else if (arg->opt->lopt == "factor-max")
            conf.radFactorMax = arg->get_arg<float>();
        else if (arg->opt->lopt == "radius-smooth")
            conf.radSmoothIters = arg->get_arg<unsigned int>();
        else if (arg->opt->lopt == "edge-smooth")
            conf.edgeSmooth = arg->get_arg<float>();
        else if (arg->opt->lopt == "planar-thres")
            conf.planarThres = arg->get_arg<float>();
        else if (arg->opt->lopt == "debug-verts") {
            conf.debugVertex = true;
            util::Tokenizer tk;
            tk.split(arg->get_arg<std::string>(), ',');
            for (std::size_t k = 0; k < tk.size(); k++) {
                unsigned int id =
                    util::string::convert<unsigned int>(std::string(tk[k]));
                conf.vertIds.insert(id);
            }
        } else if (arg->opt->lopt == "save-graphs")
            conf.saveGraphs = true;
        else if (arg->opt->lopt == "save-volume-all")
            conf.saveVolumeAll = true;
        else if (arg->opt->lopt == "save-patch-all")
            conf.savePatchAll = true;
        else if (arg->opt->lopt == "save-patch-final")
            conf.savePatchFinal = true;
        else if (arg->opt->lopt == "fixed-radius") {
            conf.usefixedRadius = true;
            conf.fixedRadius = arg->get_arg<float>();
        }
    }

    // make sure enough points will be sampled
    if (std::log(conf.radFactorMax) / std::log(conf.radFactorInc) < 7) {
        std::cout << "Warning: Less than 7 data samples will be calculated"
                  << " per vertex. You should consider increasing the"
                  << " maximum radius factor or the radius factor increment."
                  << std::endl;
        std::cout << "Continue anyway? (Enter or STRG+c)" << std::endl;
        std::cin.ignore();
    }

    // Warn the user if the sphere tessellation will be too high
    if (conf.subdivIters > 6) {
        std::cout << "Sphere tessellation will be very high, consider "
                     "lowering the subdiv argument!" << std::endl;
    }

    std::cout << "Running with config:      " << std::endl;
    std::cout << "  Subdivision Iters       " << conf.subdivIters << std::endl;
    if (conf.usefixedRadius)
        std::cout << "  Fixed radius            " << conf.fixedRadius
                  << std::endl;
    else {
        std::cout << "  RADIUS_FACTOR_INITIAL   " << conf.radFactorInitial
                  << std::endl;
        std::cout << "  RADIUS_FACTOR_INC       " << conf.radFactorInc
                  << std::endl;
        std::cout << "  RADIUS_FACTOR_MAX       " << conf.radFactorMax
                  << std::endl;
        std::cout << "  PLANAR_THRESHOLD        " << conf.planarThres
                  << std::endl;
        std::cout << "  EDGE_SMOOTH             " << conf.edgeSmooth
                  << std::endl;
        std::cout << "  VOLUME_APPROX_ITERS     " << INTERSECTION_SUBDIV_ITERS
                  << std::endl;
    }

    // Setup curvature estimator
    mve::TriangleMesh::Ptr mesh = mve::geom::load_ply_mesh(inMeshPath);
    curvest::Estimator curvest;
    curvest.setConfig(conf);
    curvest.setMesh(mesh);

    if (conf.usefixedRadius) {
        // compute single-scale curvature field
        curvest.computeFixedRadius(conf.fixedRadius);
    } else {
        // compute multi-scale curvature field
        curvest.compute();
    }
    mve::TriangleMesh::ConstPtr result = curvest.getResult();

    // write out resulting mesh:
    //   - mean curvatures are stored in the vertex value properties
    //   - radii are stored in the vertex confidence properties
    if (!conf.debugVertex) {
        mve::geom::SavePLYOptions options;
        options.write_vertex_confidences = true;
        options.write_vertex_values = true;
        options.write_vertex_normals = true;
        options.write_vertex_colors = true;
        std::stringstream ss;
        ss << outMeshPath;
        if (conf.debugVertex) ss << "_debv";

        ss << "_SPH-" << curvest.numSphereVerts();

        if (conf.usefixedRadius)
            ss << "_fixedRadius_" << conf.fixedRadius;
        else {
            ss << "_" << conf.radFactorInitial << "_" << conf.radFactorInc
               << "_" << conf.radFactorMax << "_pt" << conf.planarThres
               << "_esm" << conf.edgeSmooth;
        }

        if (!SMOOTH_FINAL_RADIUS) ss << "_notsmoothed";

        ss << "_adj_" << PERCENTILE_LOW << "_" << PERCENTILE_HIGH << ".ply";

        mve::geom::save_ply_mesh(result, ss.str(), options);
    }

    return 0;
}
