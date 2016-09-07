/*
 * Copyright (C) 2016, Patrick Seemann
 *
 * This software may be modified and distributed under the terms
 * of the BSD 3-Clause license. See the LICENSE.txt file for details.
 */

#include <iostream>
#include <limits>

#include "util/arguments.h"
#include "util/tokenizer.h"
#include "mve/mesh_io_ply.h"
#include "mve/mesh.h"
#include "mve/mesh_info.h"
#include "math/vector.h"

#include "libcurvest/utils.h"

#define DENSITY_MIN 0.0f
#define DENSITY_MAX 2.0f

struct AppSettings {
    bool linearMapping;
    bool polynomialMapping;
    bool sigmoid;
};

int main(int argc, char** argv)
{
    util::Arguments args;
    args.set_description("Creates a density field from the values stored in the 'value' field of a .ply file.");
    args.set_exit_on_error(true);
    args.set_nonopt_maxnum(2);
    args.set_nonopt_minnum(2);
    args.set_helptext_indent(22);
    args.set_usage(argv[0], "[ OPTIONS ] IN-PLY OUT");
    args.add_option(
        'l', "linear", false, "Use a linear mapping.");
    args.add_option(
        'p', "polynomial", false, "Use a polynomial mapping.");
    args.add_option(
        's', "sigmoid", false, "Use a sigmoid mapping.");
    args.parse(argc, argv);

    // Setup defaults
    AppSettings conf;
    conf.linearMapping = false;
    conf.polynomialMapping = true;
    conf.sigmoid = false;

    for (util::ArgResult const* i = args.next_result(); i != 0;
         i = args.next_result()) {
        if (i->opt == 0) continue;
        switch (i->opt->sopt) {
            case 'l':
                conf.linearMapping = true;
                conf.polynomialMapping = false;
                conf.sigmoid = false;
                break;
            case 'n':
                conf.linearMapping = false;
                conf.polynomialMapping = true;
                conf.sigmoid = false;
                break;
            case 's':
                conf.polynomialMapping = false;
                conf.linearMapping = false;
                conf.sigmoid = true;
                break;
            default:
                break;
        }
    }

    std::string infile = args.get_nth_nonopt(0);
    std::string outfile = args.get_nth_nonopt(1);
    mve::TriangleMesh::Ptr mesh = mve::geom::load_ply_mesh(infile);
    mve::TriangleMesh::ValueList& vertexValues = mesh->get_vertex_values();
    std::vector<float> density(vertexValues.size(), 0.0f);

    float min, max, absMin, absMax;
    min = std::numeric_limits<float>::max();
    max = -std::numeric_limits<float>::max();
    absMin = std::numeric_limits<float>::max();
    absMax = -std::numeric_limits<float>::max();

    for (std::size_t i = 0; i < vertexValues.size(); ++i) {
        float value = std::abs(vertexValues[i]);
        min = std::min(min, value);
        max = std::max(max, value);
        absMin = std::min(absMin, std::abs(value));
        absMax = std::max(absMax, std::abs(value));
    }

    std::cout << "Vertex values min: " << min << " max: " << max << std::endl;

    if (conf.linearMapping) {
        std::cout << "Using a linear mapping..." << std::endl;
        float range = max - min;
        for (std::size_t i = 0; i < vertexValues.size(); ++i) {
            float value =
                DENSITY_MAX * (std::abs(vertexValues[i]) - min) / range;
            density[i] = value;
        }
    } else if (conf.polynomialMapping) {
        std::cout << "Performing a polynomial mapping..." << std::endl;

        float dfl = absMin + 0.01f * (absMax - absMin); // left border
        float dfm = absMin + 0.9f * (absMax - absMin); // middle border
        float dfr = absMax; // right border
        float dfdmax = 3000.0; // max density factor
        float dfdmin = 0.1f; // min density factor
        // min density factor for values right of middle border
        float dFactorMinHighCurvature = 0.25f * dfdmax;

        std::cout << "MC abs min: " << absMin << " MC max: " << absMax
                  << std::endl;
        std::cout << "Left border: " << dfl << std::endl;
        std::cout << "Middle: " << dfm << std::endl;
        std::cout << "Right border: " << dfr << std::endl;
        std::cout << "Density fac. max: " << dfdmax << std::endl;
        std::cout << "Density fac. min: " << dfdmin << std::endl;
        std::cout << "Density fac. min (high curv.): "
                  << dFactorMinHighCurvature << std::endl;

        for (std::size_t i = 0; i < vertexValues.size(); ++i) {
            float value = std::abs(vertexValues[i]);

            if (value > dfl && value < dfm) {
                density[i] = (4.0 - 3.0 * (value - dfl) / (dfm - dfl))
                    * std::pow((value - dfl) / (dfm - dfl), 3.0) * (dfdmax-dfdmin) + dfdmin;
            } else if (value >= dfm) {
                // create exponential function where f(middle) is dFactorMax and
                // f(rightBorder) is dFactorMin
                double e, a, b, c, ymin;
                e = 2.0;
                b = -dfm;
                ymin = dFactorMinHighCurvature;
                double h = std::pow((dfr + b), e);
                c = dfdmax;
                a = (ymin - c) / h;
                density[i] = a * std::pow(value + b, e) + c;
            } else {
                density[i] = dfdmin;
            }
        }
    }else if (conf.sigmoid) {
        std::cout << "Performing a sigmoid mapping..." << std::endl;
        float dmin = 0.1f; // min density factor
        float dmax = 3000.0f; // max density factor
        float cmin = absMin + 0.01f * (absMax - absMin); // curvature min
        float cmax = absMin + 0.9f * (absMax - absMin); // curvature max
        for (std::size_t i = 0; i < vertexValues.size(); ++i) {
            float value = std::abs(vertexValues[i]);
            if (value <= cmin) {
                density[i] = dmin;
            } else if (value >= cmax) {
                density[i] = dmax;
            } else {
                float x = 2.0f * ((value - cmin) / (cmax - cmin)) - 1;
                float a = 4.0f;
                float scale = (dmax - dmin) / ( (1.0f / (1.0f + std::exp(-a))) - (1.0f / (1.0f + std::exp(a))) );
                density[i] = scale * ( (1.0f / (1.0f + std::exp(-a*x))) - (1.0f / (1.0f + std::exp(a))) ) + dmin;
            }
        }
    }

    // set density field as vertex values
    mesh->get_vertex_values() = density;

    // save mesh
    {
        mve::geom::SavePLYOptions options;
        options.write_vertex_values = true;
        options.write_vertex_confidences = false;
        std::stringstream ss;
        ss << outfile << "_density";
        if (conf.linearMapping)
            ss << "_linear";
        else if (conf.polynomialMapping)
            ss << "_polynomial";
        else if (conf.sigmoid)
            ss << "_sigmoid";
        ss << ".ply";
        mve::geom::save_ply_mesh(mesh, ss.str(), options);
    }
    return 0;
}
