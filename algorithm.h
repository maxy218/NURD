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
#include "myData.h"
#include "someGlobal.h"
#include "common.h"

using namespace std;

void get_anno_info(ifstream& in_anno, const int anno_choice, map<string, gene_info> & map_g_anno);

// get the read count of each exon
int get_exon_rd_cnt(map<string, gene_info> & map_g_info, ifstream & in_rdmap,
    ofstream & out_nurd, map<string, vector<int> > & gene_rd_cnt);

void max_isoform_bisearch(gene_info& g,int k);

double max_likelihood_given_C(gene_info& g);

double max_likelihood(gene_info& g, double alpha);

void calcuAllTheGenes(ifstream& infile, ofstream& out, double alpha);

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

void output_nurd_file_old(ofstream& out_nurd, const vector<double>& GBC, const map<string, gene_info>& map_g_anno, const map<string,vector<int> >& gene_read_count, int total_valid_read_count);


#endif // ALGORITHM_H_INCLUDED
