/*
 * =====================================================================================
 *
 *     Filename:  class.h
 *
 *  Description:  the definitions of some classes.
 *
 *    Version:  1.0
 *    Created:  02/19/2014 10:21:42 PM
 *     Revision:  none
 *     Compiler:  g++
 *
 *     Author:  Xinyun Ma
 *   Organization:  Tsinghua University
 *
 * =====================================================================================
 */


#ifndef CLASS_H_INCLUDED
#define CLASS_H_INCLUDED

#include <string>
#include <vector>

#include "const.h"
using namespace std;

struct isoform_anno{
public:
  string gene_name;
  string name;
  string chrom;
  string strand;
  _chr_coor tx_start;
  _chr_coor tx_end;
  _chr_coor cds_start;
  _chr_coor cds_end;
  _chr_coor exon_cnt;
  vector<_chr_coor> exon_starts;
  vector<_chr_coor> exon_ends;

  isoform_anno();
  isoform_anno(int anno_type, string line);
};

struct gene_info{
public:
  int iso_num;
  int exon_num;

  //here, one elem in vector represent one isoform
  string gene_name;
  vector<string> iso_name; //maybe there are multiple isoforms from one gene
  string chrom;
  vector<string> iso_chrom; //which chromosomes are the isoforms from. 
                            //maybe different isoform comes from different chromosomes, which is invalid here.
  string strand;
  vector<string> iso_strand; //which chromosomes are the isoforms from.
                             //maybe different isoform comes from different chromosomes, which is invalid here.

  //it's nothing with strand, only represent the pos on chromosome.
  // If necessary, it can be easily got from exon_start_g and exon_end_g
  _chr_coor g_start;
  _chr_coor g_end;

  vector<vector<int> > exon_iso_idx;

  //here, one elem in vector represent one exon
  vector<_chr_coor> exon_start_g;  // start position, coordinate on gene
  vector<_chr_coor> exon_end_g;  // end position, coordinate on gene

  vector<_chr_coor> exon_len;
  vector<_chr_coor> exon_g_start_l;  //original name: exonGeneLocalStarts_fromLeft. Count from 0, record the total length from start before this exon starts!
  vector<_chr_coor> exon_g_start_r;  //original name: exonGeneLocalStarts_fromRight.Count from 0, record the total length from end before this exon starts!

  //this boundary is for: binary search when dealing with read cnt.
  vector<_chr_coor> exon_g_bound;

  _chr_coor g_len;
  int tot_rd_cnt; //total read counts
  vector<int> rd_cnt;
  vector<double> GBC;
  vector<double> LBC;
  vector<double> a;//0-1 matrix
  vector<double> c;//float matrix. used in expression estimation
                   //it can be 0-1 matrix or GBC or LBC or the mixture of GBC or LBC

  vector<_chr_coor> iso_len;
  vector<double> theta;
  int is_valid;

  //following are some constructor functions.
  gene_info();
  gene_info(const gene_info& g);
  gene_info(const vector<isoform_anno>& iso_vec);

  //following are some other member function
  int if_valid();  //judge whether the gene is valid. If one gene is invalid, may be it is because:
            //1: there are 0 reads mapped to the gene
};

//the following work is to change the return type from bool to int, different return value represent different error type
bool if_gene_anno_valid(const gene_info& gene);

#endif // CLASS_H_INCLUDED
