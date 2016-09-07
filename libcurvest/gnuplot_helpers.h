/*
 * Copyright (C) 2016, Patrick Seemann
 *
 * This software may be modified and distributed under the terms
 * of the BSD 3-Clause license. See the LICENSE.txt file for details.
 */

#include <iostream>
#include <vector>
#include <sstream>

#include "defines.h"
#include "utils.h"

NAMESPACE_CURVEST_BEGIN
NAMESPACE_UTILS_BEGIN

//
// Saves a gnuplot multiplot of the given results
//
inline void saveMultiplot(std::string outPrefix, std::vector<Result> results)
{
    FILE* pipe = popen("gnuplot", "w");
    if (pipe == NULL) {
        std::cerr << "Could not open pipe for gnuplot!" << std::endl;
        return;
    }

    std::stringstream fpath;
    fpath << outPrefix << "_" << results[0].vertexId << "_multiplot";
    fpath << ".svg";
    fprintf(pipe, "set terminal svg size %d,%d enhanced font \"Verdana,8\"\n",
        results.size() * GRAPH_WIDTH, GRAPH_HEIGHT);
    fprintf(pipe, "set output \"%s\"\n", fpath.str().c_str());
    fprintf(pipe,
        "set multiplot layout 1, %d title \"Vertex %d\" font \",14\"\n",
        results.size(), results[0].vertexId);
    fprintf(pipe, "set xtics rotate\n");
    fprintf(pipe, "set size square\n");
    fprintf(pipe, "unset key\n");

    for (std::size_t i = 0; i < results.size(); ++i) {
        utils::Result r = results[i];
        if (r.fittedFunction.empty())
            fprintf(pipe, "set title \"%s\"\n", r.name.c_str());
        else
            fprintf(pipe, "set title sprintf(\"%s \\n [%s]\")\n",
                r.name.c_str(), r.fittedFunction.c_str());

        if (r.useCustomLabels) {
            fprintf(pipe, "set xlabel \"%s\"\n", r.xLabel.c_str());
            fprintf(pipe, "set ylabel \"%s\"\n", r.yLabel.c_str());
        } else {
            fprintf(pipe, "set xlabel \"x\"\n");
            fprintf(pipe, "set ylabel \"y\"\n");
        }

        // set options
        if (!r.setOptions.empty()) {
            for (std::size_t k = 0; k < r.setOptions.size(); ++k) {
                std::string opt = r.setOptions[k].first;
                fprintf(pipe, "%s\n", opt.c_str());
            }
        }

        if (r.fittedFunction.empty()) {
            fprintf(pipe, "plot '-' with points lt rgb \"blue\" \n");
            for (std::size_t l = 0; l < r.data.size(); ++l) {
                fprintf(pipe, "%f %f\n", r.data[l][0], r.data[l][1]);
            }
            fprintf(pipe, "e\n");
        } else {
            // plot points and fitted function
            fprintf(pipe, "plot %s lt rgb \"red\", ", r.fittedFunction.c_str());
            fprintf(pipe, "'-' with points lt rgb \"blue\"\n");
            for (std::size_t l = 0; l < r.data.size(); ++l) {
                fprintf(pipe, "%f %f\n", r.data[l][0], r.data[l][1]);
            }
            fprintf(pipe, "e\n");
        }

        // unset options
        if (!r.setOptions.empty()) {
            for (std::size_t k = 0; k < r.setOptions.size(); ++k) {
                std::string opt = r.setOptions[k].second;
                fprintf(pipe, "%s\n", opt.c_str());
            }
        }
    }
    fprintf(pipe, "unset multiplot\n");
    fflush(pipe);
    pclose(pipe);
}

//
// Saves a gnuplot multiplot graph (including fitted polynomial) to disk
//
inline void saveGraphToDisk(utils::Result const& origH,
    utils::Result const& normH, float const& radius, float const& error,
    math::Vec4f const& coeffs, std::string outfile, std::string extName)
{
    utils::Result tmpNormH = normH;
    utils::Result tmpOrigH = origH;

    std::stringstream nameSS;
    nameSS << tmpNormH.name << " [e = " << error << "]";
    tmpNormH.name = nameSS.str();

    std::stringstream functionSS;
    functionSS << coeffs[0];
    if (coeffs[1] > 0.0f)
        functionSS << " + " << coeffs[1] << "*x";
    else
        functionSS << " " << coeffs[1] << "*x";
    if (coeffs[2] > 0.0f)
        functionSS << " + " << coeffs[2] << "*x*x";
    else
        functionSS << " " << coeffs[2] << "*x*x";
    if (coeffs[3] > 0.0f)
        functionSS << " + " << coeffs[3] << "*x*x*x";
    else
        functionSS << " " << coeffs[3] << "*x*x*x";

    tmpNormH.fittedFunction = functionSS.str();

    // save gnuplot "set" option to display final radius
    // set arrow from 0.35,graph(0,0) to 0.35,graph(1,1) nohead
    std::stringstream ssRadiusArrow, ssUnsetArrow;
    ssRadiusArrow << "set arrow from " << radius << ",graph(0,0) to " << radius
                  << ",graph(1,1) nohead lt 6";
    ssUnsetArrow << "unset arrow";

    std::pair<std::string, std::string> showFinalRadiusOpt =
        std::make_pair(ssRadiusArrow.str(), ssUnsetArrow.str());

    tmpOrigH.setOptions.push_back(showFinalRadiusOpt);
    tmpNormH.setOptions.push_back(showFinalRadiusOpt);

    std::vector<utils::Result> multiPlotResults;
    multiPlotResults.push_back(tmpOrigH);
    multiPlotResults.push_back(tmpNormH);

    std::stringstream extSS;
    extSS << "_" << extName << "_err_" << error;
    utils::saveMultiplot(outfile + extSS.str(), multiPlotResults);
}


//
// Saves a gnuplot multiplot graph to disk
//
inline void saveGraphToDisk(utils::Result const& origH,
    utils::Result const& normH, float const& radius, std::string outfile,
    std::string extName)
{
    utils::Result tmpNormH = normH;
    utils::Result tmpOrigH = origH;

    // save gnuplot "set" option to display final radius
    // set arrow from 0.35,graph(0,0) to 0.35,graph(1,1) nohead
    std::stringstream ssRadiusArrow, ssUnsetArrow;
    ssRadiusArrow << "set arrow from " << radius << ",graph(0,0) to " << radius
                  << ",graph(1,1) nohead lt 6";
    ssUnsetArrow << "unset arrow";

    std::pair<std::string, std::string> showFinalRadiusOpt =
        std::make_pair(ssRadiusArrow.str(), ssUnsetArrow.str());

    tmpOrigH.setOptions.push_back(showFinalRadiusOpt);
    tmpNormH.setOptions.push_back(showFinalRadiusOpt);

    std::vector<utils::Result> multiPlotResults;
    multiPlotResults.push_back(tmpOrigH);
    multiPlotResults.push_back(tmpNormH);

    std::stringstream extSS;
    extSS << "_" << extName;
    utils::saveMultiplot(outfile + extSS.str(), multiPlotResults);
}

NAMESPACE_UTILS_END
NAMESPACE_CURVEST_END
