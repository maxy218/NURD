#ifndef ALGORITHM_H_INCLUDED
#define ALGORITHM_H_INCLUDED

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


//Linear search with Fabonacci Series
void max_isoform(gene& ps,int k);

//binary search
void max_isoform_bisearch(gene& ps,int k);

double max_likelihood_given_C(gene& ps);

//choice represent the models: 0=>URD, 1=>GN-URD, 2=>MN-URD, 3=>LN-URD, 4=>1-M, 5=>5-M
double max_likelihood(gene& ps, int choice);

void calcuAllTheGenes(ifstream& infile, ofstream& out, double alpha);



//return -1 if fail
//else, return the index of x
//vec is a vector of boundary of each interval
//  n+1 elements represent n intervals, the first elem is 0
//the interval is left close and right open
template <class T >
int my_bin_search(const vector<T>& vec, T x);

//linear search
int my_lin_search(const vector<int>& vec, int x);

// two vector version. It's totally different with the one vector version.
// one vector version can easily transform into two vector version.
// vec1: starts     vec2:ends
// reference: introduction to the Design and analysis of algorithms(second edition)
//  related chapter: chapter4. Chinese version, P104
template <class T >
int my_bin_search(const vector<T>& vec1, const vector<T>& vec2, T x);

// This version is for the array that the element is sorted from large element to small element. The reverse of above
template <class T >
int my_bin_search_reverse(const vector<T>& vec1, const vector<T>& vec2, T x);

// two vector and multiple return value version.
// It's similar with the above binary search version
// If there are multiple hit, return a list of recode. The list record the index of two vector.
// If there's no hit, return the null list, whose length is 0.
template <class T >
list<int> my_bin_search_multi(const vector<T>& vec1, const vector<T>& vec2, T x);


//// haven't judge whether reads have same length.
int get_read_len(ifstream& in_sam);

// annotype: the type of annotation. 1 -> refflat, 2 -> GTF
int get_exon_read_count(ifstream& in_refFlat, ifstream& in_sam, ofstream& out_nurd, bool if_chr, int annotype);

#endif // ALGORITHM_H_INCLUDED
