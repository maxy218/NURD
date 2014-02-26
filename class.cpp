/*
 * =====================================================================================
 *
 *       Filename:  class.cpp
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


#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <list>
#include <stdlib.h>
#include <algorithm>
#include <map>
#include "common.h"
#include "const.h"
#include "class.h"
using namespace std;

//gene_info: constructor functions
gene_info::gene_info(){} // do nothing. Don't caculate based on the object initialized by default constructor.

//assignment member-by-member
//  when we are assigning the container-type members, maybe smart-pointer is a better choice. 
gene_info::gene_info(const gene_info& g){
  iso_num = g.iso_num;
  exon_num = g.exon_num;

  gene_name = g.gene_name;
  iso_name = g.iso_name; // smart pointer maybe is a better choice. 
  chrom = g.chrom;
  iso_chrom = g.iso_chrom; // smart pointer maybe is a better choice. 
  strand = g.strand;
  iso_strand = g.iso_strand; 

  g_start = g.g_start;
  g_end = g.g_end;

  exon_iso_idx = g.exon_iso_idx; // smart pointer maybe is a better choice.
  exon_start_g = g.exon_start_g; // smart pointer maybe is a better choice.
  exon_end_g = g.exon_end_g; // smart pointer maybe is a better choice.
  exon_len = g.exon_len; // smart pointer maybe is a better choice.
  exon_g_start_l = g.exon_g_start_l; // smart pointer maybe is a better choice. 
  exon_g_start_r = g.exon_g_start_r; // smart pointer maybe is a better choice.
  exon_g_bound = g.exon_g_bound; // smart pointer maybe is a better choice.

  g_len = g.g_len;
  tot_rd_cnt = g.tot_rd_cnt;
  rd_cnt = g.rd_cnt; // smart pointer maybe is a better choice.
  GBC = g.GBC; // smart pointer maybe is a better choice.
  LBC = g.LBC; // smart pointer maybe is a better choice.
  a = g.a; // smart pointer maybe is a better choice.
  c = g.c; // smart pointer maybe is a better choice.

  iso_len = g.iso_len; // smart pointer maybe is a better choice.
  theta = g.theta; // smart pointer maybe is a better choice.
  is_valid = g.is_valid;
}

//0=>valid
//1=>#read=0
int gene_info::if_valid(){
    if(tot_rd_cnt <= 0){
        return 1;
    }
    return 0;
}



/////////////////////////////////////////////////////////
/////////    new class definition    /////////////
/////////////////////////////////////////////////////////

//default constructor
isoform_anno::isoform_anno():
  gene_name (""),
  name (""),
  chrom (""),
  strand (""),
  tx_start (0),
  tx_end (0),
  cds_start (0),
  cds_end (0),
  exon_cnt (0),
  exon_starts ( vector<_chr_coor>(0) ),
  exon_ends ( vector<_chr_coor>(0) )
{}

isoform_anno::isoform_anno(int anno_type, string line){
  if(anno_type == 1)// refflat annotation format
  {
    vector<string> str_vec = delimiter(line,'\t');

    gene_name = str_vec[0];
    name = str_vec[1];
    chrom = str_vec[2];
    strand = str_vec[3];
    tx_start = atoi((str_vec[4]).c_str());
    tx_end = atoi((str_vec[5]).c_str());
    cds_start = atoi((str_vec[6]).c_str());
    cds_end = atoi((str_vec[7]).c_str());
    exon_cnt = atoi((str_vec[8]).c_str());

    string tmp;
    vector<string>::iterator iter;

    exon_starts = vector<_chr_coor>(exon_cnt);
    exon_ends = vector<_chr_coor>(exon_cnt);

    vector<string> pos;
    pos = delimiter(str_vec[9], ',', exon_cnt);
    int i = 0;
    for(iter=pos.begin(); iter != pos.end(); iter++){
      exon_starts[i++] = atoi((*iter).c_str());
    }

    pos = delimiter(str_vec[10], ',', exon_cnt);
    i = 0;
    for(iter=pos.begin(); iter != pos.end(); iter++){
      exon_ends[i++] = atoi((*iter).c_str());
    }
  }
}

bool if_some_state_on(const vector<bool>& state){
  int iso_num = state.size();
  for(int i = 0; i < iso_num; i++){
    if(state[i]){
      return true;
    }
  }
  return false;
}

//the first two para are const, the last four para are the results we want
void exon_len_split(
    const vector<vector<_chr_coor> >& exon_iso_start,
    const vector<vector<_chr_coor> >& exon_iso_end,
    vector<_chr_coor>& exon_len,
    vector<vector<int> >& exon_iso_idx,
    vector<_chr_coor>& vec_splt_s,
    vector<_chr_coor>& vec_splt_e)
{
  list<int> list_s;
  list<int> list_e;

  int iso_num = exon_iso_start.size();
  for(int i = 0; i < iso_num; i++){
    int tmp_exon_num = exon_iso_start[i].size();
    for(int j = 0; j < tmp_exon_num; j++){
      list_s.push_back(exon_iso_start[i][j]);
    }
  }
  for(int i = 0; i < iso_num; i++){
    int tmp_exon_num = exon_iso_end[i].size();
    for(int j = 0; j < tmp_exon_num; j++){
      list_e.push_back(exon_iso_end[i][j]);
    }
  }

  list_s.sort();
  list_s.unique();
  list_e.sort();
  list_e.unique();

  //split begin

  //state: on or off(true or false)
  //on: current postion for i-th isoform is exon
  //off: ... is intron
  vector<bool> state(iso_num);
  for(int i = 0; i < iso_num; i++){
    state[i] = false;
  }

  //index: each elem is from 0 to exon_num-1
  //it represent the current exon ordinal for each isoform
  //until the i-th end has been dealt with, the index++, and at the same time, the state turn off
  vector<int> index(iso_num);
  for(int i = 0; i < iso_num; i++){
    index[i] = 0;
  }

  //isoform_to_deal: represent the isoform ordinal that should deal with in some operation.
  //it's mostly used when seeking the first start and so on.
  vector<int> isoform_to_deal;

  int cur_pos=0;
  while(list_s.size()>0 && list_e.size()> 0){
    if(if_some_state_on(state)){
      //start come first when there are isoform on
      //at this time, there are no isoform turn off from on, because the first end comes later than the first start
      if(list_s.front() < list_e.front() ){
        cur_pos = list_s.front();
        list_s.pop_front();
        vec_splt_s.push_back(cur_pos);
        vec_splt_e.push_back(cur_pos);

        //turn the related isoform on
        for(int i = 0; i < iso_num; i++){
          if(index[i] < exon_iso_start[i].size()){
            if(exon_iso_start[i][ index[i] ] == cur_pos){
              state[i] = true;
            }
          }
        }
      }
      //end come first when there are isoform on
      else if(list_e.front() < list_s.front() ){
        //turn the related isoform off
        //because the end is the first event of all the starts and ends, so the related isoform ending with this end must be on
        cur_pos = list_e.front();
        vec_splt_e.push_back(cur_pos);
        for(int i = 0; i < iso_num; i++){
          if(index[i] < exon_iso_start[i].size()){
            if(exon_iso_end[i][ index[i] ] == cur_pos){
              state[i] = false;
              index[i]++;
            }
          }
        }
        //need following step. If there are isoforms still on, then need add a start to vec_splt_s
        if(if_some_state_on(state)){
          vec_splt_s.push_back(cur_pos);
        }
        list_e.pop_front();
      }
      //start and end are on the same position
      else{
        cur_pos = list_s.front();
        vec_splt_s.push_back(cur_pos);
        vec_splt_e.push_back(cur_pos);
        for(int i = 0; i < iso_num; i++){
          if(index[i] < exon_iso_start[i].size()){
            if(exon_iso_end[i][ index[i] ] == cur_pos){
              state[i] = false;
              index[i]++;
            }
            else if(exon_iso_start[i][ index[i] ] == cur_pos){
              state[i] = true;
            }
          }
        }
        list_s.pop_front();
        list_e.pop_front();
      }
    }
    //if no isoform on, then find the first start as the new start point
    else{
      cur_pos=list_s.front();
      vec_splt_s.push_back(cur_pos);
      for(int i = 0; i < iso_num; i++){
        if(index[i] < exon_iso_start[i].size()){
          if(exon_iso_start[i][ index[i] ] == cur_pos){
            state[i] = true;
          }
        }
      }
      list_s.pop_front();
    }
  }

  //deal with the left ends
  while(list_e.size()> 0){
    //turn the related isoform off
    cur_pos = list_e.front();
    vec_splt_e.push_back(cur_pos);
    for(int i = 0; i < iso_num; i++){
      if(exon_iso_end[i][ index[i] ] == cur_pos){
        state[i] = false;
        index[i]++;
      }
    }
    //need following step. If there are isoforms still on, then need add a start to vec_splt_s
    if(if_some_state_on(state)){
      vec_splt_s.push_back(cur_pos);
    }
    list_e.pop_front();
  }
  //split done!

  //exon indicator begin!
  for(int i = 0; i < iso_num; i++){
    vector<int> tmp_vec;
    exon_iso_idx.push_back(tmp_vec);
  }
  int total_exon_num = vec_splt_s.size();

  vector<_chr_coor>::iterator iter_vec_s = vec_splt_s.begin();
  vector<_chr_coor>::iterator iter_vec_e = vec_splt_e.begin();
  while(iter_vec_s != vec_splt_s.end()){
    exon_len.push_back((*iter_vec_e)-(*iter_vec_s));
    iter_vec_s++;
    iter_vec_e++;
  }
  for(int i = 0; i < iso_num; i++){
    int iso_exon_index = 0;
    bool if_on = false;
    for(iter_vec_s = vec_splt_s.begin(), iter_vec_e = vec_splt_e.begin();
      iter_vec_s != vec_splt_s.end(), iter_vec_e != vec_splt_e.end(), iso_exon_index < exon_iso_start[i].size();
      iter_vec_s++,iter_vec_e++)
    {
      if((*iter_vec_s) < exon_iso_start[i][iso_exon_index]){
        exon_iso_idx[i].push_back(0);
      }
      else if((*iter_vec_s) >= exon_iso_start[i][iso_exon_index] && (*iter_vec_e) < exon_iso_end[i][iso_exon_index]){
        exon_iso_idx[i].push_back(1);
      }
      else if((*iter_vec_e) == exon_iso_end[i][iso_exon_index]){
        exon_iso_idx[i].push_back(1);
        iso_exon_index++;
      }
    }
    for(;iter_vec_s != vec_splt_s.end(),iter_vec_e != vec_splt_e.end(); iter_vec_s++, iter_vec_e++){
      exon_iso_idx[i].push_back(0);
    }
  }
  //exon indicator done!
}



gene_info::gene_info(const vector<isoform_anno>& iso_vec){
  vector<isoform_anno>::const_iterator iter_gene = iso_vec.begin();
  vector<_chr_coor> vec_exon_splt_s;
  vector<_chr_coor> vec_exon_splt_e;

  //use tmp, because there exists exon spliting event
  vector<vector<_chr_coor> > tmp_exon_iso_start;
  vector<vector<_chr_coor> > tmp_exon_iso_end;

  gene_name = (*iter_gene).gene_name;

  iter_gene = iso_vec.begin();
  while(iter_gene != iso_vec.end()){
    vector<_chr_coor>::iterator tmp_iter;
    iso_name.push_back((*iter_gene).name);
    chrom = (*iter_gene).chrom;
    iso_chrom.push_back(chrom);
    strand = (*iter_gene).strand;
    iso_strand.push_back(strand);

    tmp_exon_iso_start.push_back((*iter_gene).exon_starts);
    tmp_exon_iso_end.push_back((*iter_gene).exon_ends);
    iter_gene++;
  }
  iso_num = iso_name.size();

  vector<_chr_coor> tmp_exon_g_len;
  vector<vector<_chr_coor> > tmp_exon_iso_idx;
  ///////////////////////////////////////////////////////////////////////////
  exon_len_split(
    tmp_exon_iso_start,
    tmp_exon_iso_end,
    tmp_exon_g_len,
    tmp_exon_iso_idx,
    vec_exon_splt_s,
    vec_exon_splt_e);
  ///////////////////////////////////////////////////////////////////////////

  exon_num = tmp_exon_g_len.size();

  size_t idx = 0;
  size_t sz = vec_exon_splt_s.size();
  for(idx = 0; idx < sz; ++idx){
    exon_start_g.push_back(vec_exon_splt_s[idx]);
    exon_end_g.push_back(vec_exon_splt_e[idx]);
  }

  g_start = exon_start_g[0];
  g_end = exon_end_g[exon_end_g.size()-1];

  if(strand == "+"){
    exon_len = tmp_exon_g_len;
    exon_iso_idx = tmp_exon_iso_idx;
  }
  else{
    exon_len = vector<_chr_coor>(exon_num);
    for(int i = exon_num-1; i >= 0; i--){
      exon_len[i] = tmp_exon_g_len[exon_num - 1 - i];
    }
    for(int j = 0; j < iso_num; j++){
      vector<_chr_coor> tmp_iso_ind = vector<_chr_coor>(exon_num);
      for(int i = exon_num - 1; i >= 0; i--){
        tmp_iso_ind[i] = tmp_exon_iso_idx[j][exon_num - 1 - i];
      }
      exon_iso_idx.push_back(tmp_iso_ind);
    }

    // if neg strand, reverse the exon boundary
    reverse(exon_start_g.begin(), exon_start_g.end());
    reverse(exon_end_g.begin(), exon_end_g.end());
  }

  g_len = 0;
  exon_g_bound = vector<_chr_coor>(exon_num+1);
  exon_g_bound[0] = 0;
  for(int i = 0; i < exon_num; i++){
    g_len += exon_len[i];
    exon_g_bound[i+1] = exon_g_bound[i] + exon_len[i];
  }

  exon_g_start_l = vector<_chr_coor>(exon_num);
  exon_g_start_r = vector<_chr_coor>(exon_num);
  exon_g_start_l[0] = 0;
  exon_g_start_r[0] = 0;
  if(exon_num >= 2){
    for(int i = 1; i < exon_num; i++){
      exon_g_start_l[i] = exon_g_start_l[i-1] + exon_len[i-1];
      exon_g_start_r[i] = exon_g_start_r[i-1] + exon_len[exon_num-i];
    }
  }

  // initialize the vector members of gene_info
  rd_cnt = vector<int>(exon_num, 0);

  iso_len = vector<_chr_coor>(iso_num, 0);
  for(size_t i = 0; i < iso_num; ++i){
    for(size_t j = 0; j < exon_num; ++j){
      iso_len[i] += exon_iso_idx[i][j]*exon_len[j];
    }
  }

  a = vector<double>(iso_num * exon_num, 0.0);
  for(size_t i = 0; i < iso_num; ++i){
    for(size_t j = 0; j < exon_num; ++j){
      a[i* exon_num + j] = exon_iso_idx[i][j];
    }
  }
  c = GBC = LBC = a;

  int min_iso_len = iso_len[0];
  for(size_t i = 1; i < iso_num; ++i){
    if(iso_len[i] < min_iso_len){
      min_iso_len = iso_len[i];
    }
  }

  theta = vector<double>(iso_num, 0.0);
  for(size_t i = 0; i < iso_num; ++i){
     theta[i] = 1.0/(double)min_iso_len;
  }

  // the member of is_valid should be set after get the read counts.
}

//the following work is to change the return type from bool to int, different return value represent different error type
bool if_gene_anno_valid(const gene_info& gene){
  bool if_valid = true;
  //  if there're duplicate isoforms
  //  sorted list may be not efficient. Maybe hash is faster!
  //  hash: just judge whether every isoform name is new, having no duplicate
  if(gene.iso_num > 1){
    list<string> list_iso_name;

    for(int i = 0; i < gene.iso_num; i++){
      list_iso_name.push_back(gene.iso_name[i]);
    }
    list_iso_name.sort();
    list_iso_name.unique();
    if(list_iso_name.size() != gene.iso_num){
      return false;
    }
  }

  //if having different orient
  if(gene.iso_num > 1){
    for(int i = 0; i < gene.iso_num - 1; i++){
      if(gene.iso_strand[i] != gene.iso_strand[i+1]){
        return false;
      }
    }
  }

  //if from different chrome
  if(gene.iso_num > 1){
    for(int i = 0; i < gene.iso_num - 1; i++){
      if(gene.iso_chrom[i] != gene.iso_chrom[i+1]){
        return false;
      }
    }
  }

  //if isoform comes from chromosome whose name has underline "_"
  //delete the isoforms on such chrom: chr6_qbl_hap2
  //method: delimit the chrom column, if there exists "_", then delete it
  size_t found;
  for(int i = 0; i < gene.iso_num; i++){
    found = gene.iso_chrom[i].find('_');
    if(found != string::npos){
      return false;
    }
  }
  return true;
}
