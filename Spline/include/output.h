/**
 * @file output.h
 * @author OuyangShangke
 * @brief Funcs in this file are all for output. 
 * @details The aim of those funcs is to make the output coordinate with the gramma of Octave code, 
 * so that one can immediately execute the output file by Octave. There's no need to open the output file and edit it again.
 * @version 0.1
 * @date 2021-11-23
 * @copyright Copyright (c) 2021
 * 
 */


#ifndef __OUTPUT_H__
#define __OUTPUT_H__

#include "Config.h"
#include <fstream>
#include <iostream>

/**
 * @brief Output a std::map<double, Vec<double,_Dim>> object.
 * @details The output code can plot a curve, which is produced by the data in the std::map object.
 * Besides, one can change the appearance of the curve(include legend, linetype) by simply introducing different param of this func.
 * 
 * @tparam _Dim the dimention of the Vec object.
 * @param _path the path of the output file
 * @param _t the std::map object
 * @param _index index to avoid the same varname
 * @param _legend legend of this curve
 * @param _linetype(dispensable) the line type of the curve
 * @return 0 for OK
 */
template<u_int _Dim>
int out4map(const STRING& _path, const std::map<double, Vec<double,_Dim>>& _t,
            int _index, const STRING& _legend, const STRING& _linetype = "")
{
    std::ofstream outfile;
    outfile.open(_path, std::ios::app);
    vector2_double p(_Dim);
    for(auto it = _t.begin() ; it != _t.end() ; it ++)
    {
        Vec<double, _Dim> tmp = it -> second;
        for(int i = 0 ; i < _Dim ; i ++)
            p[i].push_back(tmp[i]);
    } 
    for(int i = 0 ; i < _Dim ; i ++)
    {
        outfile << "x" << i << "_" << _index << "=[ ";
        for(auto it = p[i].begin() ; it < p[i].end() ; it ++)
            outfile << *it << " "; 
        outfile << "];" << std::endl;
    }
    outfile << "plot(";
    for(int i = 0 ; i < _Dim ; i ++)
        outfile << "x" << i << "_" << _index << ",";
    outfile << "'"<<_linetype<<";"<<_legend<<";')" << std::endl;
    return 0;
}

/**
 * @brief A specialization of the above template func out4map.
 * @details This specialization is for the Vec which is one-dimention.
 * Every element in the map can be regarded as a point in R^2.
 * The output code of this func can plot a curve which connecting all these points.
 * 
 * @param _path the path of the output file
 * @param _t the std::map object
 * @param _index index to avoid the same varname
 * @param _legend legend of this curve
 * @param _linetype(dispensable) the line type of the curve
 * @return 0 for OK
 */
template<>
int out4map<1>(const STRING& _path, const std::map<double, Vec<double,1>>& _t, 
               int _index, const STRING& _legend, const STRING& _linetype)
{
    std::ofstream outfile;
    outfile.open(_path, std::ios::app);
    outfile << "x_" << _index << "=[ ";
    for(auto it = _t.begin() ; it != _t.end() ; it ++)
        outfile << it -> first << " ";    
    outfile << "];" << std::endl;
    outfile << "y_" << _index << "=[ ";
    for(auto it = _t.begin() ; it != _t.end() ; it ++)
        outfile << (it -> second)[0] << " ";    
    outfile << "];" << std::endl;
    outfile << "plot(x_"<<_index<<",y_"<<_index<<",'"<<_linetype<<";"<<_legend<<";')" << std::endl;

    outfile.close();
    return 0;
}

/**
 * @brief Output a std::map<double, double> object.
 * @details Every element in the map can be regarded as a point in R^2.
 * The output code of this func can plot a curve which connecting all these points.
 * 
 * @param _path the path of the output file
 * @param _t the std::map object
 * @param _index index to avoid the same varname
 * @param _legend legend of this curve
 * @param _linetype(dispensable) the line type of the curve 
 * @return 0 for OK
 */
int out4map(const STRING& _path, const M_DOUBLE_DOUBLE& _t, int _index, const STRING& _legend,
            const STRING& _linetype = "")
{
    std::ofstream outfile;
    outfile.open(_path, std::ios::app);
    outfile << "x_" << _index << "=[ ";
    for(auto it = _t.begin() ; it != _t.end() ; it ++)
        outfile << it -> first << " ";    
    outfile << "];" << std::endl;
    outfile << "y_" << _index << "=[ ";
    for(auto it = _t.begin() ; it != _t.end() ; it ++)
        outfile << it -> second << " ";    
    outfile << "];" << std::endl;
    outfile << "plot(x_"<<_index<<",y_"<<_index<<",'"<<_linetype<<";"<<_legend<<";')" << std::endl;

    outfile.close();
    return 0;
}

/**
 * @brief create the file if it dosen't exist, or clean the file if it already exists.
 * 
 * @param _path path of the output file
 * @return
 */
int clean_file(const STRING& _path)
{
    /// clean the file
    std::fstream clean(_path, std::ios::out);
    clean.close();
    return 0;
}

/**
 * @brief Output Octave code for creating a new figure.
 * 
 * @param _path path of the output file 
 * @param _index index of figure
 * @param _hold(dispensible) a bool variable to show whether output the code "hold on". \n
 * This param is dispensible. It's default value is true which means output the "hold on".
 * @param _tilte(dispensible) title of the plot.
 * @return
 */
int out_new_fig(const STRING& _path, int _index = 1, bool _hold = true, STRING _title = STRING(""))
{
    std::ofstream outfile;
    outfile.open(_path, std::ios::app);
    outfile << "figure(" << _index << ")" << std::endl;
    if(_hold)
        outfile << "hold on" << std::endl;
    if(!_title.empty())
        outfile << "title('" << _title << "');" << std::endl;
    outfile.close();
    return 0;
}

/**
 * @brief Output some other Octave codes.
 * 
 * @param _path path of the output file
 * @param _settings It's a vector of std::string, and every element in the vector is a line of Octave code.
 * @return
 */
int out_other_settings(const STRING& _path, const vector_string& _settings)
{
    std::ofstream outfile;
    outfile.open(_path, std::ios::app);
    for(auto it = _settings.begin() ; it < _settings.end() ; it ++)
        outfile << *it << std::endl;
    outfile.close();
    return 0;
}

/**
 * @brief Output the code for saving plot.
 * 
 * @param _path path of the output file
 * @param _save_path path of the saving plots.
 * @return
 */
int out_save(const STRING& _path, const STRING& _save_path)
{
    std::ofstream outfile;
    outfile.open(_path, std::ios::app);
    outfile << "legend()" << std::endl;
    outfile << "saveas(gcf, '"<<_save_path << "')" << std::endl;
    outfile.close();
    return 0;
}

#endif