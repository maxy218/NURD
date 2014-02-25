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

//return -1 if fail
//else, return the index of x
//vec is a vector of boundary of each interval
//  n+1 elements represent n intervals, the first elem is 0
//the interval is left close and right open
template <class T >
int bin_search(const vector<T> & vec, const T & x);

// two vector version. It's totally different with the one vector version.
// one vector version can easily transform into two vector version.
// vec1: starts     vec2:ends
// reference: introduction to the Design and analysis of algorithms(second edition)
//  related chapter: chapter4. Chinese version, P104
template <class T >
int bin_search(const vector<T>& vec1, const vector<T>& vec2, const T & x);

// This version is for the array that the element is sorted from large element to small element. The reverse of above
template <class T >
int bin_search_reverse(const vector<T>& vec1, const vector<T>& vec2, const T & x);

// two vector and multiple return value version.
// It's similar with the above binary search version
// If there are multiple hit, return a list of recode. The list record the index of two vector.
// If there's no hit, return the null list, whose length is 0.
template <class T >
list<int> bin_search_multi(const vector<T>& vec1, const vector<T>& vec2, const T & x);

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
