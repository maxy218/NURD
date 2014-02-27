/*
 * =====================================================================================
 *
 *       Filename:  algorithm.cpp
 *
 *    Description:  some algorithms
 *
 *        Version:  1.0
 *        Created:  02/24/2014 11:02:13 PM
 *       Revision:  none
 *       Compiler:  g++
 *
 *         Author:  Xinyun Ma
 *   Organization:  Tsinghua University
 *
 * =====================================================================================
 */



#include <fstream>
#include <vector>
#include <iostream>
#include <sstream>
#include <algorithm>
#include <math.h>

#include "class.h"
#include "const.h"
#include "common.h"
#include "algorithm.h"

using namespace std;

//deal with refflat format annotation file.
void get_anno_refflat(ifstream& in_anno, map<string, vector<isoform_anno> >& g_iso_anno_map){
  map<string, gene_info>::iterator iter_map_g_anno;

  string temp_line;
  while(getline(in_anno, temp_line)){
    isoform_anno g(1, temp_line); //use refflat format annotation file.
    g_iso_anno_map[g.gene_name].push_back(g);
  }
}

//only deal with the exon annotation. CDS and start/end_codon is ignored.
// 0 -> chr; 1 -> data source; 2 -> function; 3 -> start; 4 -> end; 5 -> score; 6 -> strand; 7 -> phase; 8 -> gene id and trans id
void get_anno_GTF(ifstream& in_anno, map<string, vector<isoform_anno> >& g_iso_anno_map){
  map<string, gene_info>::iterator iter_map_g_anno;

  string temp_line;
  map<string, isoform_anno> iso_anno_map; // key: transcript name, value: transcript annotation.
  map<string, isoform_anno>::iterator iter_map_iso_anno; // key: transcript name, value: transcript annotation.
  while(getline(in_anno, temp_line)){
    vector<string> str_vec = delimiter(temp_line,'\t');
    vector<string> gene_trans_id_vec = delimiter(str_vec[8],'\"');

    _chr_coor start_pos, end_pos;

    if(str_vec[2] == "exon") // only deal with exon.
    {
      start_pos = atoi(str_vec[3].c_str()) - 1; // "-1" is because the GTF starts from 1, not 0 (which is refFlat style.)
      end_pos = atoi(str_vec[4].c_str()); // the end position is the same with refflat.

      // trans name and chr name are needed to identify a transcript.
      // some trans may come from different chromosome.
      string iso_chr_combined = gene_trans_id_vec[3] + "\t" + str_vec[0]; 

      // if has been dealt before.
      if( (iter_map_iso_anno = iso_anno_map.find(iso_chr_combined)) != iso_anno_map.end()) 
      {
        iter_map_iso_anno -> second.exon_starts.push_back(start_pos);
        iter_map_iso_anno -> second.exon_ends.push_back(end_pos);
      }
      else // this trans has not been dealt.
      {
        iso_anno_map[ iso_chr_combined ] = isoform_anno();

        isoform_anno & iso_anno = iso_anno_map[ iso_chr_combined ];
        iso_anno.gene_name = gene_trans_id_vec[1];
        iso_anno.name = gene_trans_id_vec[3];
        iso_anno.chrom = str_vec[0];
        iso_anno.strand = str_vec[6];

        iso_anno.exon_starts.push_back(start_pos);
        iso_anno.exon_ends.push_back(end_pos);
      }
    }
  }

  for(iter_map_iso_anno = iso_anno_map.begin(); iter_map_iso_anno != iso_anno_map.end(); iter_map_iso_anno++)
  {
    isoform_anno & iso_anno = iter_map_iso_anno->second;
    sort(iso_anno.exon_starts.begin(), iso_anno.exon_starts.end());
    sort(iso_anno.exon_ends.begin(), iso_anno.exon_ends.end());
    unique(iso_anno.exon_starts.begin(), iso_anno.exon_starts.end());
    unique(iso_anno.exon_ends.begin(), iso_anno.exon_ends.end());

    iso_anno.exon_cnt= iso_anno.exon_starts.size();
    iso_anno.tx_start= iso_anno.exon_starts[0];
    iso_anno.tx_end = iso_anno.exon_ends[ iso_anno.exon_cnt - 1 ];

    g_iso_anno_map[iso_anno.gene_name].push_back( iso_anno );
  }
}

void get_anno_info(ifstream & in_anno, const unsigned int anno_choice, 
    map<string, gene_info> & map_g_anno){
  map<string, gene_info>::iterator iter_map_g_anno;
  map<string, vector<isoform_anno> > g_iso_anno_map;
  map<string, vector<isoform_anno> >::iterator iter_g_iso;

  if(anno_choice == 1){
    get_anno_refflat(in_anno, g_iso_anno_map);
  }
  else{
    get_anno_GTF(in_anno, g_iso_anno_map);
  }

  for(iter_g_iso = g_iso_anno_map.begin(); iter_g_iso != g_iso_anno_map.end(); iter_g_iso++){
    //map_g_anno[(*iter_g_iso).first] = gene_info( g_iso_anno_map[(*iter_g_iso).first] );
    map_g_anno[(*iter_g_iso).first] = gene_info( iter_g_iso -> second );
  }

  iter_map_g_anno = map_g_anno.begin();
  while(iter_map_g_anno != map_g_anno.end()){
    if(!if_gene_anno_valid( (*iter_map_g_anno).second) ){
      map_g_anno.erase(iter_map_g_anno++);
    }
    else{
      iter_map_g_anno++;
    }
  }
}

// output the data to the nurd file, which is a temporary file.
void output_nurd_file(const vector<double>& GBC, const map<string, gene_info>& map_g_anno, 
    size_t tot_valid_rd_cnt, ofstream& out_nurd){
  // output GBC curve
  output_vector<double>(GBC, out_nurd, '\t');

  // output valid read number and gene number
  out_nurd << tot_valid_rd_cnt << '\t';
  out_nurd << map_g_anno.size() << endl;

  // output each gene's detail information
  map<string, gene_info>::const_iterator iter_map_g_anno; // const iterator
  for(iter_map_g_anno = map_g_anno.begin(); iter_map_g_anno != map_g_anno.end(); ++iter_map_g_anno){
    string gene_name = (*iter_map_g_anno).second.gene_name;

    //output as "nurd" format
    //first line: basic information
    out_nurd << gene_name << '\t';
    out_nurd << (*iter_map_g_anno).second.exon_num << '\t';
    out_nurd << (*iter_map_g_anno).second.iso_num << endl;

    //second line: isoform name
    output_vector((*iter_map_g_anno).second.iso_name, out_nurd, '\t');

    //third line: exon length
    output_vector((*iter_map_g_anno).second.exon_len, out_nurd, '\t');

    //forth line: read count in exon
    output_vector((*iter_map_g_anno).second.rd_cnt, out_nurd, '\t');

    //other lines: gene structure
    output_2D_vector((*iter_map_g_anno).second.exon_iso_idx, out_nurd, '\t');
  }
}

// get the map relation of (chr -> (pos -> gene))
void get_map_chr_pos_gene(map<string, gene_info> & map_g_info, 
    map<string, map<_chr_coor, string> > & map_chr_pos_gene){
  map<string, gene_info>::iterator iter_map_g_info;

  for(iter_map_g_info = map_g_info.begin(); iter_map_g_info != map_g_info.end(); ++iter_map_g_info){
    gene_info & g = iter_map_g_info -> second;

    //frist map
    if(map_chr_pos_gene.find(g.chrom) == map_chr_pos_gene.end()){
      map_chr_pos_gene[g.chrom] = map<_chr_coor, string>();
      map_chr_pos_gene[g.chrom][g.g_start] = g.gene_name;
    }
    else{
      map_chr_pos_gene[g.chrom][g.g_start] = g.gene_name;
    }
  }
}

// get the map relation of (chr -> all the start position on this chromosome)
void get_map_chr_start_pos(map<string, gene_info> & map_g_info,
    map<string, vector<_chr_coor> > & map_chr_start_pos){
  map<string, vector<_chr_coor> >::iterator iter_map_chr_pos;
  map<string, gene_info>::iterator iter_map_g_info;

  for(iter_map_g_info = map_g_info.begin(); iter_map_g_info != map_g_info.end(); ++iter_map_g_info){
    gene_info & g = iter_map_g_info -> second;

    //second map
    if(map_chr_start_pos.find(g.chrom) == map_chr_start_pos.end()){
      map_chr_start_pos[g.chrom] = vector<_chr_coor>(1, g.g_start);
    }
    else{
      map_chr_start_pos[g.chrom].push_back(g.g_start);
    }
  }

  //sort the vector before binary search
  for(iter_map_chr_pos = map_chr_start_pos.begin(); iter_map_chr_pos != map_chr_start_pos.end(); 
      ++iter_map_chr_pos){
    sort(iter_map_chr_pos -> second.begin(), iter_map_chr_pos -> second.end());
  }
}

// get the map relation of (chr -> all the end position on this chromosome)
void get_map_chr_end_pos(map<string, gene_info> & map_g_info, map<string, 
    map<_chr_coor,string> > & map_chr_pos_gene,
    map<string, vector<_chr_coor> > & map_chr_start_pos,
    map<string, vector<_chr_coor> > & map_chr_end_pos){
  map<string, vector<_chr_coor> >::iterator iter_map_chr_pos;
  // get map_chr_end_pos
  // because the start pos is sorted, so the correspond end pos can't be sure sorted. So the corresponding vector should be got based on the start pos
  for(iter_map_chr_pos = map_chr_start_pos.begin(); iter_map_chr_pos != map_chr_start_pos.end();
      ++iter_map_chr_pos){
    size_t chr_s_pos_size = iter_map_chr_pos -> second.size();
    map_chr_end_pos[iter_map_chr_pos -> first] = vector<_chr_coor>(chr_s_pos_size);
    vector<_chr_coor> & ref_end_pos = map_chr_end_pos[iter_map_chr_pos -> first];
    vector<_chr_coor> & ref_start_pos = map_chr_end_pos[iter_map_chr_pos -> first];
    map<_chr_coor, string> & ref_map_pos_gene = map_chr_pos_gene[iter_map_chr_pos -> first];
    for(size_t i = 0; i < chr_s_pos_size; i++){
      string & gene_name = ref_map_pos_gene[ iter_map_chr_pos -> second[i] ];
      ref_end_pos[i] = map_g_info[gene_name].g_end;
    }
  }
}

// get the read count of each exon
ofstream out_nurd("tmp_nurd.txt");
size_t get_exon_rd_cnt(map<string, gene_info> & map_g_info, ifstream & in_rdmap, 
    size_t & tot_valid_rd_cnt, vector<double> & GBC){
  clock_t start_time,end_time;
  time_t cur_time;
  start_time = clock();

  clock_t tmp_start, tmp_end;
  stringstream ss (stringstream::in | stringstream::out);

  map<string, gene_info>::iterator iter_map_g_info;

  end_time = clock();
  ss << "gene annotation time: " << ((double)end_time-start_time)/CLOCKS_PER_SEC << " seconds.\n";
  std_output_with_time(ss.str());
  ss.str("");
  start_time = end_time;

  ss << "gene annotation done!" << endl;
  std_output_with_time(ss.str());
  ss.str("");

  //gene read count begins
  ////////////////////////////////////////
  /////// read count starts!
  tot_valid_rd_cnt = 0;

  // bin number in GBC. GBC is calculated at the same time of reads counting.
  vector<int> int_GBC = vector<int>(GBC_BIN_NUM, 0);

  size_t outlier_read_cnt = 0;

  //// first key is chr name, second key is gene start pos, second value is gene name. The first value is a map container.
  //// This is map nest definition.
  map<string, map<_chr_coor, string> > map_chr_pos_gene;

  // the key is chr name, the value is a vector, in which are the gene start positions on this chr.
  // no gene name. Gene names can be easily got from the former hash map.
  // this vector should be sorted before binary search.
  map<string, vector<_chr_coor> > map_chr_start_pos;
  map<string, vector<_chr_coor> > map_chr_end_pos;
  map<string, vector<_chr_coor> >::iterator iter_map_chr_pos_vec;

  get_map_chr_pos_gene(map_g_info, map_chr_pos_gene);
  get_map_chr_start_pos(map_g_info, map_chr_start_pos);
  get_map_chr_end_pos(map_g_info, map_chr_pos_gene, map_chr_start_pos, map_chr_end_pos);

  string temp_line;

  //deal with header
  in_rdmap.seekg(0, ios::beg);
  while(getline(in_rdmap, temp_line)){
    if(temp_line[0] != '@'){
      break;
    }
  }
  //deal with mapped reads
  unsigned int RD_FLAG_MASK_REVERSE_MAP = 0x10;
  // vector<string> sam_column = vector<string>(4);
  vector<string> sam_column = vector<string>(10);

  do{
    unsigned int invalid_type = 0;//0: valid, 1: not gene region, 2: not exon region 3: multi gene 4: invalid chrome

    //extract the map information from sam file
    // delimiter(sam_column,temp_line,'\t',4,true); // because only the first 4 fields are used. reduce the time that was
    delimiter(sam_column, temp_line, '\t', 10, true); // because only the first 10 fields are used. reduce the time that was
    string & chrName = sam_column[2];
    unsigned int read_Flag = atoi(sam_column[1].c_str());

    if(chrName == "*"){
      continue;
    }

    iter_map_chr_pos_vec = map_chr_start_pos.find(chrName);
    if(iter_map_chr_pos_vec == map_chr_start_pos.end()){
      continue;
    }

    _chr_coor read_pos = atoi(sam_column[3].c_str());
    if( (read_Flag & RD_FLAG_MASK_REVERSE_MAP) != 0 ){ // reverse read
      read_pos += sam_column[9].size();
    }

    vector<int> vec_gene_idx = bin_search_multi(map_chr_start_pos[chrName],map_chr_end_pos[chrName],read_pos);

    if(vec_gene_idx.size() == 0){ // if no gene cover the read, deal with the next read
      continue;
    }

    bool if_map_to_gene = false;
    bool if_multi_map = false;
    int exon_idx = -1;
    string g_name;
    for(size_t idx = 0; idx < vec_gene_idx.size(); ++idx){
      g_name = map_chr_pos_gene[chrName][ map_chr_start_pos[chrName][ vec_gene_idx[idx] ] ];

      gene_info & g_info = map_g_info[g_name];
      if(g_info.strand == "+"){
        exon_idx = bin_search(g_info.exon_start_g, g_info.exon_end_g, read_pos);
      }
      ///// if on the negtive strand, should use the reverse search
      else{
        exon_idx = bin_search_reverse(g_info.exon_start_g, g_info.exon_end_g, read_pos);
      }

      if(exon_idx != -1){
        if(!if_map_to_gene){
          if_map_to_gene = true;
        }
        else{
          if_multi_map = true;
          break;
        }
      }
    }
    if(if_multi_map){
      continue;
    }

    gene_info & g_info = map_g_info[g_name];
    if(exon_idx != -1){
      tot_valid_rd_cnt++;
      int gene_read_pos = -1; // initialized as an invalid position.
      if(g_info.strand == "+"){
        g_info.rd_cnt[exon_idx]++;

        //for GBC, read pos
        gene_read_pos = g_info.exon_g_start_l[exon_idx] + read_pos - g_info.exon_start_g[exon_idx];
      }
      else{
        g_info.rd_cnt[exon_idx]++;

        //for GBC, read pos
        // still exon_g_start_l, if you figure out, it's obvious.
        gene_read_pos = g_info.exon_g_start_l[exon_idx] + g_info.exon_end_g[exon_idx] - read_pos;
      }

      //For GBC
      //only use the genes with single isoform
      if(g_info.iso_name.size() == 1)
      {
        _chr_coor gene_len = g_info.g_len;
        if( gene_read_pos < gene_len && gene_read_pos >= 0 ){
          int_GBC[ (int)(gene_read_pos*GBC_BIN_NUM)/gene_len ]++;
        }
        else if(gene_read_pos == gene_len){
          int_GBC[ GBC_BIN_NUM-1 ]++;
        }
        else{
          outlier_read_cnt++;
        }
      }
    }
  }while(getline(in_rdmap,temp_line));
  //////////  read count ends!
  //////////////////////////////////////////

  //////////////////////////////////
  //////// get GBC
  double total_GBC = 0;
  for(size_t i = 0; i < GBC_BIN_NUM; i++){
    total_GBC += int_GBC[i];
  }
  for(size_t i = 0; i < GBC_BIN_NUM; i++){
    GBC[i] = ((double)int_GBC[i]*GBC_BIN_NUM)/total_GBC;
  }
  //////////////////////////////////

  end_time = clock();
    ss << "read count time: " << ((double)end_time-start_time)/CLOCKS_PER_SEC << " seconds.\n";
    std_output_with_time(ss.str());
    ss.str("");
    start_time = end_time;

  // update the is_valid information of each gene.
  for(iter_map_g_info = map_g_info.begin(); iter_map_g_info != map_g_info.end(); ++iter_map_g_info){
    iter_map_g_info->second.tot_rd_cnt = sum_vector(iter_map_g_info->second.rd_cnt);
  }

  //output the information to nurd file.
//  output_nurd_file(GBC, map_g_info, tot_valid_rd_cnt, out_nurd);

  return 0;
}

double get_log_likelihood(const gene_info& g){
  size_t N = g.exon_num;
  size_t M = g.iso_num;

  double likeli = 0.0;
  for(size_t j = 0; j < N; j++){
    double tempSum1 = 0.0;
    double tempSum2 = 0.0;
    for(size_t i = 0; i < M; i++){
      double tmp_double = 0.0;
      tmp_double = g.c[i*N+j] * g.theta[i];
      tempSum1 += g.exon_len[j] * tmp_double;
      tempSum2 += tmp_double;
    }
    likeli = likeli - g.tot_rd_cnt*tempSum1 + g.rd_cnt[j]*log(g.exon_len[j]*g.tot_rd_cnt*tempSum2+EPSILON);
  }
  return likeli;
}

double get_gradient_of_log_likelihood(const gene_info& g, size_t i){
  size_t N = g.exon_num;
  size_t M = g.iso_num;

  double gradient = 0.0;
  for(size_t j = 0; j < N; j++){
    gradient -= g.tot_rd_cnt*g.exon_len[j]*g.c[i*N+j];
    if(g.rd_cnt[j]*g.c[i*N+j] != 0)
    {
      double tmp_sum = EPSILON;
      for(size_t k = 0; k < M; k++)
      {
        tmp_sum += g.c[k*N+j]*g.theta[k];
      }
      gradient += g.rd_cnt[j]*g.c[i*N+j]/tmp_sum;
    }
  }
  return gradient;
}

void max_isoform_bisearch(gene_info& g, size_t k){
  const double LOCAL_EPSILON = 1e-8;
  const double LOCAL_EPSILON_GRADIENT = 1e-8;
  double left, right;

  double effective_iso_len = 0.0;
  size_t N = g.exon_num;
  for(size_t j = 0; j < N; j++)
  {
    effective_iso_len += g.c[k*N+j]*g.exon_len[j];
  }
  left = 0;
  right = 1/effective_iso_len;

  ///////// quick check of the boundary. If gradient near zero < 0 => 0; if gradient near right boundary > 0 => right boundary.
  //check left boundary
  g.theta[k] = left + LOCAL_EPSILON;
  if( get_gradient_of_log_likelihood(g, k) < 0)
  {
    g.theta[k] = left;
    return;
  }
  //check right boundary
  g.theta[k] = right - LOCAL_EPSILON;
  if( get_gradient_of_log_likelihood(g, k) > 0)
  {
    g.theta[k] = right;
    return;
  }
  ///////// quick check end

  double tmp = (right+left)/2;
  double gradient = 1.0;
  // while(right - left > LOCAL_EPSILON && (gradient > LOCAL_EPSILON_GRADIENT || gradient < -LOCAL_EPSILON_GRADIENT ))
  while(right - left > LOCAL_EPSILON)
  {
    g.theta[k] = tmp;
    gradient = get_gradient_of_log_likelihood(g, k);
    if( gradient < 0 )
    {
      right = tmp;
    }
    else if(gradient > 0)
    {
      left = tmp;
    }
    else
    {
      break;
    }
    tmp = (left+right)/2;
  }
}

double max_likelihood_given_C(gene_info& g){
  size_t M = g.iso_num;

  double f_value = get_log_likelihood(g);
  size_t iteration = 0;

  while(true){
    for(size_t i = 0; i < M; i++){
      max_isoform_bisearch(g, i);
    }

    double new_f_value = get_log_likelihood(g);
    if(fabs(new_f_value - f_value) < EPSILON){
      return new_f_value;
    }
    f_value = new_f_value;

    if(iteration < MAX_ITER_NUM){
      iteration++;
    }
    else{
      break;
    }
  }
  return f_value;
}

void get_GBC_matrix(gene_info& g, const vector<double> & GBC);
void get_LBC_matrix(gene_info& g);

double max_likelihood(gene_info& g, double alpha, const vector<double> & GBC){
  size_t M = g.iso_num;
  size_t N = g.exon_num;

  get_GBC_matrix(g, GBC);
  get_LBC_matrix(g);
  double rel_len = 0; // relative length
  for(size_t i = 0; i < M; i++){
    for(size_t j = 0; j < N; j++){
      rel_len = (double)g.exon_len[j]/g.iso_len[i];
      g.c[i*N+j] = ( g.GBC[i*N+j]*alpha + g.LBC[i*N+j]*(1-alpha) )/rel_len/GBC_BIN_NUM;
    }
  }
  return max_likelihood_given_C(g);
}

void get_LBC_curve(gene_info& g, vector<double>& LBC){
  size_t M = g.iso_num;
  size_t N = g.exon_num;

  double tempSum = 0.0;
  double LBC_sum = 0.0;
  for(size_t j = 0; j < N; j++){
    tempSum = 0.0;
    for(size_t i = 0; i < M; i++){
      tempSum += g.theta[i]*g.a[i*N+j];
    }
    if(tempSum < EPSILON){
      tempSum += EPSILON;
    }
    LBC[j] = g.rd_cnt[j]/(g.exon_len[j]*tempSum);
    LBC_sum += LBC[j];
  }

  for(size_t j = 0; j < N; j++){
    LBC[j] = LBC[j]*N/LBC_sum;
  }
}

void get_curve_from_hist(const vector<double> & hist_h, const vector<double> & hist_l,
    const vector<double> & len, vector<double>& area){

  size_t hist_size = hist_h.size();
  size_t len_size = len.size();
  size_t tot_hist_len = sum_vector(hist_l);
  size_t tot_len = sum_vector(len);

  for(size_t i = 0; i < len_size; i++){
    area[i] = 0.0;
  }

  size_t idx1 = 0, idx2 = 0; // idx1: index of hist, idx2: index of len
  double cur_hist_len = hist_l[idx1];
  double cur_len; 
  for(; idx2 < len_size; ++idx2){
    area[idx2] = 0.0;
    cur_len = (double)len[idx2] / tot_len * tot_hist_len;
    while(idx1 < hist_size && cur_len > cur_hist_len){
      area[idx2] += cur_hist_len * hist_h[idx1];
      cur_len -= cur_hist_len;
      cur_hist_len = hist_l[++idx1]; 
    }
    if(idx1 >= hist_size){
      return;
    }
    area[idx2] += cur_len * hist_h[idx1];
    cur_hist_len -= cur_len;
  }
}

void get_GBC_matrix(gene_info& g, const vector<double> & GBC){
  size_t M = g.iso_num;
  size_t N = g.exon_num;
  vector<double> area = vector<double>(N, 0);
  vector<double> len = vector<double>(N, 0);
  vector<double> GBC_l = vector<double>(GBC_BIN_NUM, 1);
 
  for(size_t i = 0; i < M; i++){
    for(size_t j = 0; j < N; j++){
      len[j] = g.exon_len[j]*g.exon_iso_idx[i][j];
    }
    get_curve_from_hist(GBC, GBC_l, len, area);
    double tot_area = sum_vector(area);
    for(size_t j = 0; j < N; j++){
      g.GBC[i*N+j] = area[j];
    }
  }
}

void get_LBC_matrix(gene_info& g){
  size_t M = g.iso_num;
  size_t N = g.exon_num;

  vector<double> LBC_h = vector<double>(N, 0.0);
  vector<double> LBC_l = vector<double>(N, 0.0);
  get_LBC_curve(g, LBC_h);

  vector<double> area= vector<double>(N, 0.0);
  vector<double> len = vector<double>(N, 0);

  for(size_t j = 0; j < N; j++){
    for(size_t i = 0; i < M; i++){
      LBC_l[j] += g.exon_iso_idx[i][j];
    }
    LBC_l[j] *= g.exon_len[j];
  }

  for(size_t i = 0; i < M; i++){
    for(size_t j = 0; j < N; j++){
      len[j] = (double)g.exon_len[j]*g.exon_iso_idx[i][j];
    }
    get_curve_from_hist(LBC_h, LBC_l, len, area);
    double tot_area = sum_vector(area);
    for(size_t j = 0; j < N; j++){
      g.LBC[i*N+j] = area[j] * GBC_BIN_NUM / tot_area;
    }
  }
}

void calcuAllTheGenes(map<string, gene_info> & map_g_info, 
    size_t tot_valid_rd_cnt, double alpha, const vector<double> & GBC, ofstream& out){
  map<string, gene_info>::iterator iter_map_g_info = map_g_info.begin();
  for(; iter_map_g_info != map_g_info.end(); ++iter_map_g_info){
    gene_info & g = iter_map_g_info -> second;

    g.is_valid = g.if_valid();
    if(g.is_valid == 1){
      out << g.gene_name << "\t" << g.iso_num << "\t" << g.tot_rd_cnt << "\t";
      for(size_t i = 0; i < g.iso_num; i++){
        out << g.iso_name[i] <<",";
      }
      out << "\t";
      for(size_t j = 0; j < g.iso_num; j++){
      out << 0 << ",";
      }
      out << "\t" << 0; // total Theta, total expression.
      out << "\n";
    }
    else if(g.is_valid == 0){
      out << g.gene_name << "\t" << g.iso_num << "\t" << g.tot_rd_cnt << "\t";
      for(size_t i = 0; i < g.iso_num; i++){
        out << g.iso_name[i] << ",";
      }
      out << "\t";

    //////  expression estimation
      max_likelihood(g, alpha, GBC);
      double totalTheta = 0.0;
      for(size_t ii = 0; ii < g.iso_num; ii++){
        out << g.theta[ii]*g.tot_rd_cnt/tot_valid_rd_cnt*1e9 << ",";
        totalTheta += g.theta[ii];
      }
      out << "\t" << totalTheta*g.tot_rd_cnt/tot_valid_rd_cnt*1e9;
      out << "\n";
    }
  }
}

