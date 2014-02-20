/*
 * =====================================================================================
 *
 *       Filename:  class.h
 *
 *    Description:  the definitions of some classes.
 *
 *        Version:  1.0
 *        Created:  02/19/2014 10:21:42 PM
 *       Revision:  none
 *       Compiler:  g++
 *
 *         Author:  Xinyun Ma
 *   Organization:  Tsinghua University
 *
 * =====================================================================================
 */


#ifndef CLASS_H_INCLUDED
#define CLASS_H_INCLUDED

#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <list>
#include <stdlib.h>
#include <algorithm>
#include <map>
#include "common.h"
#include "someGlobal.h"
using namespace std;

//gene class represent the real gene information
class gene{
public:
    string gene_name;
    int N;
    int M;
    int w;//total read number
    vector<string> iso_name;//maybe there are multiple isoforms from one gene
    vector<int> l;
    vector<int> x;
    vector<double> a;//0-1 matrix
    vector<double> c;//float matrix. used in calculation, can be 0-1 matrix or GBC or LBC or the mixture of GBC or LBC

    vector<double> GBC;
    vector<double> LBC;
    vector<int> iso_len;

    vector<double> theta;

    //0=>valid
    //1=>#read=0
    int ifValid;

    gene(
		string _gene_name,
		int _N, // # exon
		int _M, // # isoform
		int _w,//total read number
		const vector<string>& _isoform_name,//maybe there are multiple isoforms from one gene
		const vector<int>& _l,
		const vector<int>& _x,
		const vector<double>& _a,//0-1 matrix
		const vector<double>& _c,//float matrix. used in calculation, can be 0-1 matrix or GBC or LBC or the mixture of GBC or LBC

		const vector<double>& _GBC,
		const vector<double>& _LBC,
		const vector<int>& _isoform_length,

		const vector<double>& _theta,
		int _ifValid
	);

    gene(const gene& g1);
	gene();

    bool getGeneData(ifstream& infile);
    int judgeIfValid();	//judge whether the gene is valid. If one gene is invalid, may be it is because:
						//1: there are 0 reads mapped to the gene

};

class isoform_anno{
public:
    string g_name;
    string name;
    string chrom;
    string strand;
    _chr_coor txStart;
    _chr_coor txEnd;
    _chr_coor cdsStart;
    _chr_coor cdsEnd;
    _chr_coor exonCount;
    vector<_chr_coor> exonStarts;
    vector<_chr_coor> exonEnds;

    void output();
    isoform_anno();
    isoform_anno(int anno_type, string line);

};

class gene_anno{
public:
    string g_name;

	int iso_num;
	int exon_num;

    //here, one elem in vector represent one isoform
    vector<string> iso_name;
    vector<string> iso_chrom;
	string chrom;
    vector<string> iso_strand;
	string strand;

	//it's nothing with strand, only represent the pos on chromosome.
	// If necessary, it can be easily got from exon_start_g and exon_end_g
	_chr_coor g_start;
	_chr_coor g_end;

	vector<vector<_chr_coor> > exon_iso_idx;

    //here, one elem in vector represent one exon
    vector<_chr_coor> exon_start_g;	// start position, coordinate on gene
    vector<_chr_coor> exon_end_g;	// end position, coordinate on gene

    vector<_chr_coor> exon_len;
    vector<_chr_coor> exon_g_start_l;	//original name: exonGeneLocalStarts_fromLeft. Count from 0, record the total length from start before this exon starts!
    vector<_chr_coor> exon_g_start_r;	//original name: exonGeneLocalStarts_fromRight.Count from 0, record the total length from end before this exon starts!

	//this boundary is for: binary search when dealing with read cnt.
	vector<_chr_coor> exon_g_bound;

    _chr_coor g_len;

	// function
    gene_anno();
    gene_anno(list<isoform_anno>& g_set);
};

bool if_some_state_on(vector<bool>& state);

//the first two para are const, the last four para are the results we want
void exon_len_split(
		const vector<vector<int> >& exon_iso_start,
		const vector<vector<int> >& exon_iso_end,
		vector<int>& exon_len,
		vector<vector<int> >& exon_iso_idx,
		list<int>& list_split_s,
		list<int>& list_split_e);

vector<double> getGBC(
	map<string, list<_chr_coor> > gene_read_pos_map,
	map<string,gene_anno> gene_anno_map,
	int num_bins
	);

//the following work is to change the return type from bool to int, different return value represent different error type
bool if_gene_valid(const gene_anno& gene);

#endif // CLASS_H_INCLUDED
