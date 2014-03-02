/*
 * =====================================================================================
 *
 *       Filename:  algorithm.h
 *
 *    Description:  some algorithms
 *
 *        Version:  1.0
 *        Created:  02/24/2014 11:03:31 PM
 *       Revision:  none
 *       Compiler:  g++
 *
 *         Author:  Xinyun Ma
 *   Organization:  Tsinghua University
 *
 * =====================================================================================
 */


#ifndef ALGORITHM2_H_INCLUDED
#define ALGORITHM2_H_INCLUDED

#include <cmath>
#include <fstream>
#include <vector>

#include <boost/unordered_map.hpp>

#include "class.h"
#include "common.h"
#include "const.h"

using namespace std;
using namespace boost;

void get_anno_info(ifstream& in_anno, const unsigned int anno_choice,
    unordered_map<string, gene_info> & map_g_anno);

// get the read count of each exon
size_t get_exon_rd_cnt(unordered_map<string, gene_info> & map_g_info, ifstream & in_rdmap,
    size_t & tot_valid_rd_cnt, vector<double> & GBC);

void express_estimate(unordered_map<string, gene_info> & map_g_info,
    size_t tot_valid_rd_cnt, double alpha, const vector<double> & GBC, ofstream& out);

#endif // ALGORITHM_H_INCLUDED
