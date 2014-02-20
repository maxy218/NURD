#include <fstream>
#include <vector>
#include <iostream>
#include <sstream>
#include <algorithm>
#include <math.h>

#include "class.h"
#include "myData.h"
#include "someGlobal.h"
#include "common.h"
#include "algorithm2.h"

using namespace std;

//deal with refflat format annotation file.
void get_anno_refflat(ifstream& in_anno, map<string, gene_anno> & map_g_anno, map<string,list<isoform_anno> >& isoform_anno_map){
  map<string,gene_anno>::iterator iter_map_g_anno;
  map<string,list<isoform_anno> >::iterator iter_gene_anno_map;

  string temp_line;
  while(getline(in_anno, temp_line)){
    isoform_anno g(1,temp_line); //use refflat format annotation file.
    isoform_anno_map[g.g_name].push_back(g);
  }
}

//only deal with the exon annotation. CDS and start/end_codon is ignored.
void get_anno_GTF(ifstream& in_anno, map<string, gene_anno> & map_g_anno, map<string,list<isoform_anno> >& isoform_anno_map){
  map<string,gene_anno>::iterator iter_map_g_anno;
  map<string,list<isoform_anno> >::iterator iter_gene_anno_map;

  string temp_line;
  map<string, isoform_anno> dealt_trans_anno_map; // key: transcript name, value: transcript annotation.
  map<string, isoform_anno>::iterator trans_anno_iter; // key: transcript name, value: transcript annotation.
  while(getline(in_anno, temp_line)){
    // 0 -> chr; 1 -> data source; 2 -> function; 3 -> start; 4 -> end; 5 -> score; 6 -> strand; 7 -> phase; 8 -> gene id and trans id
    vector<string> str_vec = delimiter(temp_line,'\t');
    vector<string> gene_trans_id_vec = delimiter(str_vec[8],'\"');

    _chr_coor start_pos, end_pos;

    if(str_vec[2] == "exon") // only deal with exon.
    {
      start_pos = atoi(str_vec[3].c_str()) - 1; // "-1" is because the GTF starts from 1, not 0 (which is refFlat style.)
      end_pos = atoi(str_vec[4].c_str()); // the end position is the same with refflat.

      string key_trans_name_chr = gene_trans_id_vec[3] + "\t" + str_vec[0]; // trans name and chr name are needed to identify a transcript. some trans may come from different chromosome.

      if(dealt_trans_anno_map.find( key_trans_name_chr ) != dealt_trans_anno_map.end()) // if has been dealt before.
      {
        dealt_trans_anno_map[ key_trans_name_chr ].exonStarts.push_back(start_pos);
        dealt_trans_anno_map[ key_trans_name_chr ].exonEnds.push_back(end_pos);
      }
      else // this trans has not been dealt.
      {
        dealt_trans_anno_map[ key_trans_name_chr ] = isoform_anno();
        dealt_trans_anno_map[ key_trans_name_chr ].g_name = gene_trans_id_vec[1];
        dealt_trans_anno_map[ key_trans_name_chr ].name = gene_trans_id_vec[3];
        dealt_trans_anno_map[ key_trans_name_chr ].chrom = str_vec[0];
        dealt_trans_anno_map[ key_trans_name_chr ].strand = str_vec[6];

        dealt_trans_anno_map[ key_trans_name_chr ].exonStarts.push_back(start_pos);
        dealt_trans_anno_map[ key_trans_name_chr ].exonEnds.push_back(end_pos);
      }
    }
  }

  for(trans_anno_iter = dealt_trans_anno_map.begin(); trans_anno_iter != dealt_trans_anno_map.end(); trans_anno_iter++)
  {
    sort(trans_anno_iter->second.exonStarts.begin(),trans_anno_iter->second.exonStarts.end());
    sort(trans_anno_iter->second.exonEnds.begin(),trans_anno_iter->second.exonEnds.end());
    unique(trans_anno_iter->second.exonStarts.begin(),trans_anno_iter->second.exonStarts.end());
    unique(trans_anno_iter->second.exonEnds.begin(),trans_anno_iter->second.exonEnds.end());

    trans_anno_iter->second.exonCount = trans_anno_iter->second.exonStarts.size();
    trans_anno_iter->second.txStart = trans_anno_iter->second.exonStarts[0];
    trans_anno_iter->second.txEnd = trans_anno_iter->second.exonEnds[ trans_anno_iter->second.exonCount - 1 ];

    isoform_anno_map[trans_anno_iter->second.g_name].push_back( trans_anno_iter->second );
  }
}

void get_anno_info(ifstream& in_anno, const int anno_choice, map<string, gene_anno> & map_g_anno){
  map<string,gene_anno>::iterator iter_map_g_anno;
  map<string,list<isoform_anno> > isoform_anno_map;
  map<string,list<isoform_anno> >::iterator iter_gene_anno_map;

  if(anno_choice == 1){
    get_anno_refflat(in_anno, map_g_anno, isoform_anno_map);
  }
  else{
    get_anno_GTF(in_anno, map_g_anno, isoform_anno_map);
  }

  for(iter_gene_anno_map = isoform_anno_map.begin(); iter_gene_anno_map != isoform_anno_map.end(); iter_gene_anno_map++){
    map_g_anno[(*iter_gene_anno_map).first] = gene_anno( isoform_anno_map[(*iter_gene_anno_map).first] );
  }

  iter_map_g_anno = map_g_anno.begin();
  while(iter_map_g_anno != map_g_anno.end()){
    if(!if_gene_valid( (*iter_map_g_anno).second) ){
      map_g_anno.erase(iter_map_g_anno++);
    }
    else{
      iter_map_g_anno++;
    }
  }
}

