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

#include <fstream>
#include <vector>
#include <iostream>
#include <sstream>
#include <math.h>

#include "class.h"
#include "const.h"
#include "common.h"

using namespace std;

void get_anno_info(ifstream& in_anno, const int anno_choice, map<string, gene_info> & map_g_anno);

// get the read count of each exon
int get_exon_rd_cnt(map<string, gene_info> & map_g_info, ifstream & in_rdmap,
    size_t & tot_valid_rd_cnt, vector<double> & GBC);


void max_isoform_bisearch(gene_info& g,int k);

double max_likelihood_given_C(gene_info& g);

double max_likelihood(gene_info& g, double alpha, const vector<double> & GBC);

void calcuAllTheGenes(map<string, gene_info> & map_g_info,
    size_t tot_valid_rd_cnt, double alpha, const vector<double> & GBC, ofstream& out);

double get_log_likelihood(const gene_info& g);

double get_gradient_of_log_likelihood(const gene_info& g, int i);

//return the first line of urd file
void get_GBC_bin(ifstream& infile);

void get_LBC_curve(gene_info& g, double* LBC);

void get_LBC_matrix(gene_info& g);

void get_GBC_matrix(gene_info& g, const vector<double> & GBC);

//N:number of bin
//exon_N:number of exon
void get_Curve_from_bin(const vector<double>& bin, int N, const vector<int>& length,
    int exon_N, vector<double>& Area);


#endif // ALGORITHM_H_INCLUDED
