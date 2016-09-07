/*
 * Copyright (C) 2016, Patrick Seemann
 *
 * This software may be modified and distributed under the terms
 * of the BSD 3-Clause license. See the LICENSE.txt file for details.
 */

#include <iostream>
#include <assert.h>

#include "util/arguments.h"
#include "util/tokenizer.h"
#include "mve/mesh_io_ply.h"
#include "mve/mesh.h"
#include "mve/mesh_info.h"
#include "mve/image_base.h"
#include "mve/image_io.h"

#include "colormap.h"
#include "math/vector.h"

#define COLORBAR_WIDTH 2000
#define COLORBAR_HEIGHT 200

struct AppSettings {
    bool useVertexValues;
    bool useVertexConfidences;
    float min, max;
    bool useCustomRange;
    bool saveColorbar;
    RGB low;
    RGB high;
    bool useCustomColors;
};

int main(int argc, char** argv)
{
    util::Arguments args;
    args.set_description("Colorize a mesh from its vertex values/confidences using diverging colormaps.");
    args.set_exit_on_error(true);
    args.set_nonopt_maxnum(2);
    args.set_nonopt_minnum(2);
    args.set_helptext_indent(22);
    args.set_usage(argv[0], "[ OPTIONS ] IN-PLY OUT");
    args.add_option('c', "confidences", false, "Use vertex confidences.");
    args.add_option('\0', "range", true, "Range of the colormap (comma-separated)");
    args.add_option('\0', "colors", true, "RGB colors for low and high ends of the colormap (e.g. 1.0,0.0,0.0,0.0,0.0,1.0 for red-blue)");
    args.add_option('\0', "colorbar", false, "Save the colorbar as a PNG.");
    args.parse(argc, argv);

    AppSettings conf;
    conf.useVertexValues = true;
    conf.useVertexConfidences = false;
    conf.min = 0.0f;
    conf.max = 0.0f;
    conf.useCustomRange = false;
    conf.saveColorbar = false;
    conf.useCustomColors = false;

    for (util::ArgResult const* arg = args.next_result(); arg != 0;
         arg = args.next_result()) {
        if (arg->opt == 0) continue;

        if (arg->opt->lopt == "confidences") {
            conf.useVertexConfidences = true;
            conf.useVertexValues = false;
        }
        else if (arg->opt->lopt == "range") {
            std::string str = arg->get_arg<std::string>();
            util::Tokenizer tk;
            tk.split(str, ',');
            if (tk.size() != 2) {
                std::cerr << "Need exaclty two comma separated values specifing the range!" << std::endl;
                std::exit(1);
            }
            conf.min = util::string::convert<float>(std::string(tk[0]));
            conf.max = util::string::convert<float>(std::string(tk[1]));
            conf.useCustomRange = true;

        }
        else if (arg->opt->lopt == "colorbar") {
            conf.saveColorbar = true;
        }
        else if (arg->opt->lopt == "colors") {
            std::string str = arg->get_arg<std::string>();
            util::Tokenizer tk;
            tk.split(str, ',');
            if (tk.size() != 6) {
                std::cerr << "Need exaclty six comma separated values specifing low and high rgb colors!" << std::endl;
                std::exit(1);
            }
            conf.low._r = util::string::convert<float>(std::string(tk[0]));
            conf.low._g = util::string::convert<float>(std::string(tk[1]));
            conf.low._b = util::string::convert<float>(std::string(tk[2]));

            conf.high._r = util::string::convert<float>(std::string(tk[3]));
            conf.high._g = util::string::convert<float>(std::string(tk[4]));
            conf.high._b = util::string::convert<float>(std::string(tk[5]));
            conf.useCustomColors = true;

        }
    }

    std::string infile = args.get_nth_nonopt(0);
    std::string outfile = args.get_nth_nonopt(1);
    mve::TriangleMesh::Ptr mesh = mve::geom::load_ply_mesh(infile);
    mve::TriangleMesh::VertexList const& verts = mesh->get_vertices();
    mve::TriangleMesh::ColorList& vcolor = mesh->get_vertex_colors();


    std::vector<float> & attrib = mesh->get_vertex_values();
    if (conf.useVertexConfidences)
        attrib = mesh->get_vertex_confidences();

    assert(verts.size() == attrib.size());
    vcolor.clear();
    vcolor.resize(attrib.size());


    float fmin, fmax;
    if (conf.useCustomRange) {
        fmin = conf.min;
        fmax = conf.max;
    } else {
        // Find min/max value of the attribute.
        fmin = std::numeric_limits<float>::max();
        fmax = -std::numeric_limits<float>::max();
        for (std::size_t i = 0; i < attrib.size(); ++i)
        {
            fmin = std::min(fmin, attrib[i]);
            fmax = std::max(fmax, attrib[i]);
        }
    }
    std::cout << "Mapping from " << fmin << " to " << fmax << std::endl;

    RGB low, high;
    Diverging cmap(fmin, fmax);
    if (fmin < 0.0f && fmax > 0.0f)
    {
        if (conf.useCustomColors) {
            low = conf.low;
            high = conf.high;
        } else {
            // RdBu color mapping
            MSH red(80.0, 1.08, 0.5);
            MSH blue(80.0, 1.08, -1.1);
            low = blue.toRGB();
            high = red.toRGB();
        }
        cmap.setLow(low._r, low._g, low._b);
        cmap.setHigh(high._r, high._g, high._b);
        cmap.setMidpoint(0.0);
    }
    else {
        if (conf.useCustomColors) {
            low = conf.low;
            high = conf.high;
        } else {
            // blue color mapping
            low = RGB(1.0, 1.0, 1.0);
            high = RGB(0.0, 0.0, 1.0);
        }
        cmap.setLow(low._r, low._g, low._b);
        cmap.setHigh(high._r, high._g, high._b);
        cmap.setMidpoint(fmin);
    }

    // Apply colormap
    for (std::size_t i = 0; i < vcolor.size(); ++i) {
        RGB rgb = cmap.colormap(attrib[i]);
        vcolor[i] = math::Vec4f(rgb._r, rgb._g, rgb._b);
    }

    // Save mesh
    {
        mve::geom::SavePLYOptions options;
        options.write_vertex_values = true;
        options.write_vertex_confidences = true;
        options.write_vertex_colors = true;
        std::stringstream ss;
        ss << outfile;
        if (conf.useCustomRange) {
            ss << "_range" << fmin << "to" << fmax;
        }
        ss << "_colored.ply";
        mve::geom::save_ply_mesh(mesh, ss.str(), options);
    }

    // Save colorbar
    if (conf.saveColorbar) {
        mve::ByteImage::Ptr im = mve::ByteImage::create(COLORBAR_WIDTH, COLORBAR_HEIGHT, 3);
        float step = (fmax - fmin) / im->width();
        float value = fmin;
        for (std::size_t x = 0; x < im->width(); ++x) {
            RGB rgb = cmap.colormap(value);
            for (std::size_t y = 0; y < im->height(); ++y) {
                im->at(x, y, 0) = (unsigned char) (rgb._r * 255.0f);
                im->at(x, y, 1) = (unsigned char) (rgb._g * 255.0f);
                im->at(x, y, 2) = (unsigned char) (rgb._b * 255.0f);
            }
            value += step;
        }
        std::stringstream name;
        name << outfile << "_range" << fmin << "to" << fmax << "_colorbar.png";
        mve::image::save_png_file(im, name.str());
    }


    return 0;
}
