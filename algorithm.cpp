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
#include "algorithm.h"

using namespace std;

bool get_anno_info(ifstream& in_anno, const int anno_choice, map<string, gene_anno> & map_g_anno);

void max_isoform_bisearch(gene& ps,int k){
	double local_epsilon = 1e-8;
	double local_epsilon_gradient = 1e-8;
	double left, right;
	
	double effective_iso_len = 0.0;
	int N = ps.N;
	for(int j = 0; j < N; j++)
	{
		effective_iso_len += ps.c[k*N+j]*ps.l[j];
	}
	left = 0; 
	right = 1/effective_iso_len;
	
	///////// quick check of the boundary. If gradient near zero < 0 => 0; if gradient near right boundary > 0 => right boundary.
	//check left boundary
	ps.theta[k] = left + local_epsilon;
	if( get_gradient_of_log_likelihood(ps,k) < 0)
	{
		ps.theta[k] = left;
		return;
	}
	//check right boundary
	ps.theta[k] = right - local_epsilon;
	if( get_gradient_of_log_likelihood(ps,k) > 0)
	{
		ps.theta[k] = right;
		return;
	}
	///////// quick check end
	
	double tmp = (right+left)/2;
	double gradient = 1.0;
	// while(right - left > local_epsilon && (gradient > local_epsilon_gradient || gradient < -local_epsilon_gradient ))
	while(right - left > local_epsilon)
	{
		ps.theta[k]=tmp;
		gradient = get_gradient_of_log_likelihood(ps,k);
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

double max_likelihood_given_C(gene& ps){
    int M = ps.M;

    double f_value = get_log_likelihood(ps);
    int iteration = 0;

    while(true){
        for(int i = 0; i < M; i++){
            max_isoform_bisearch(ps, i);
        }
        double new_f_value = get_log_likelihood(ps);
        if(fabs(new_f_value - f_value) < epsilon){
            return new_f_value;
        }
        f_value = new_f_value;

        if(iteration < max_iteration){
			iteration++;
		}
        else{
            break;
        }
    }
    return f_value;
}

double max_likelihood(gene& ps, double alpha){
    bool if_nurd = true;//defualt is true, unless the choice is 0

    int iteration = 0;//means the total number of iterate;
    int iter = 0;     //represent the current number of iterate;

    double f_value = get_log_likelihood(ps);
    double new_f_value = 0.0;

	get_GBC_matrix(ps);

	do{
		get_LBC_matrix(ps);
		int M = ps.M;
		int N = ps.N;

		vector<double> tempLength = vector<double>(M*N);
		for(int i = 0; i < M; i++){
			for(int j = 0; j < N; j++){
				tempLength[i*N+j] = (double)ps.l[j]/ps.iso_len[i];
				ps.c[i*N+j] = ( ps.GBC[i*N+j]*alpha + ps.LBC[i*N+j]*(1-alpha) )/tempLength[i*N+j]/bin_N;
			}
		}

		new_f_value = max_likelihood_given_C(ps);
		iter++;
	}while(iter <= iteration);
    return f_value;
}

void calcuAllTheGenes(ifstream& infile, ofstream& out, double alpha){

	int buf_size = 4096;
    char buf[buf_size];
	
    int totalRead;
    infile >> totalRead;
    int totalGeneNum;
    infile >> totalGeneNum;

    string s;
    while(true){
        gene g; 
		
		if( !g.getGeneData(infile) ){ // get the data for gene g.
			break;
		}

        if(g.ifValid == 1){
            out << g.gene_name << "\t" << g.M << "\t" << g.w << "\t";
            for(int i = 0; i < g.M; i++){
                out << g.iso_name[i] <<",";
            }
            out << "\t";
			for(int j = 0; j < g.M; j++){
				out << 0 << ",";
			}
			out << "\t" << 0; // total Theta, total expression.
			out << "\n";
        }
        else if(g.ifValid == 0){
            out << g.gene_name << "\t" << g.M << "\t" << g.w << "\t";
            for(int i = 0; i < g.M; i++){
                out << g.iso_name[i] << ",";
            }
            out << "\t";

			//////	expression estimation
			max_likelihood(g,alpha);
			double totalTheta=0.0;
			for(int ii = 0; ii < g.M; ii++){
				out << g.theta[ii]*g.w/totalRead*1e9 << ",";
				totalTheta += g.theta[ii];
			}
			out << "\t" << totalTheta*g.w/totalRead*1e9;
			out << "\n";
        }
    }
}

//return -1 if fail
//else, return the index of x
//vec is a vector of boundary of each interval
//  n+1 elements represent n intervals, the first elem is 0
//the interval is left close and right open
template <class T >
int my_bin_search(const vector<T>& vec, T x){
    int left = 0;
    int right = vec.size()-1;
    int mid = left + (right - left)/2;//find the middle interval, size-1 for there are n+1 elems in vector

    if(x >= vec[right] || x < vec[left]){
        return -1;
    }
    else{
        while(right - left>1){
            if(x < vec[mid]){
                right = mid;
                mid = left + (right - left)/2;
            }
            else if(x >= vec[mid]){
                left = mid;
                mid = left + (right - left)/2;
            }
            else{
                return mid; 
            }
        }
        return left;
    }
}

//linear search
int my_lin_search(const vector<int>& vec, int x){
    int right=1;
    if(x>=vec[vec.size()-1]||x<vec[0]){
        return -1;
    }
    while(x>=vec[right]){
        right++;
    }
    return right-1;
}


// two vector version. It's totally different with the one vector version.
// one vector version can easily transform into two vector version.
// vec1: starts     vec2:ends
// reference: introduction to the Design and analysis of algorithms(second edition)
//  related chapter: chapter4. Chinese version, P104
template <class T >
int my_bin_search(const vector<T>& vec1, const vector<T>& vec2, T x){
    int l=0;
    int r=vec1.size()-1;
    int m;

    //if x is out of bound, return -1. false
    if(x>=vec2[r]||x<vec1[l]){
        return -1;
    }
    else{
        while(l<=r){
            m = (l+r)/2;
            if(x>=vec1[m]&&x<vec2[m]){
                return m;
            }
            else if(x<vec1[m]){
                r=m-1;
            }
            else{
                l=m+1;
            }
        }
        return -1;
    }
}

// This version is for the array that the element is sorted from large element to small element. The reverse of above
template <class T >
int my_bin_search_reverse(const vector<T>& vec1, const vector<T>& vec2, T x){
    int l=0;
    int r=vec1.size()-1;
    int m;

    //if x is out of bound, return -1. false
    if(x>=vec2[l]||x<vec1[r]){
        return -1;
    }
    else{
        while(l<=r){
            m = (l+r)/2;
            if(x>=vec1[m]&&x<vec2[m]){
                return m;
            }
            else if(x<vec1[m]){
                l=m+1;
            }
            else{
                r=m-1;
            }
        }
        return -1;
    }
}

// two vector and multiple return value version.
// It's similar with the above binary search version
// If there are multiple hit, return a list of recode. The list record the index of two vector.
// If there's no hit, return the null list, whose length is 0.
template <class T >
list<int> my_bin_search_multi(const vector<T>& vec1, const vector<T>& vec2, T x){
    int l = 0;
    int r = vec1.size()-1;
    int m;
    list<int> result;

    //if x is out of bound, return -1. false
    if(x >= vec2[r] || x < vec1[l]){
        return result;
    }
    else{
        while(l <= r){
            m = (l+r)/2;
            if(x >= vec1[m] && x < vec2[m]){
                result.push_back(m);
                break;
            }
            else if(x < vec1[m]){
                r = m-1;
            }
            else{
                l = m+1;
            }
        }

        int m_bak = m;
		int flank_gene = 10; // allowing for flank 10 genes to search the covered genes.

        // search the region before m. Stop when the x is larger than the right bound
        // not perfect. Because only starts are sorted, so it can't guarantee to find all the proper interval.
        m--;
		int l_flank = 0;
        while(m>=0 && l_flank < flank_gene){
			l_flank++;
            if(x >= vec1[m] && x < vec2[m]){
                result.push_back(m);
            }
			m--;
        }
        // search the region after m. Stop when the x is smaller than the left bound
        // This half should be perfect.
        m = m_bak + 1;
        int total_size = vec1.size();
		int r_flank = 0;
        while(m < total_size && r_flank < flank_gene){
			r_flank++;
            if(x >= vec1[m] && x < vec2[m]){
                result.push_back(m);
            }
			m++;
        }

        result.sort();
        return result;
    }
}

void output_refflat_format_anno(ofstream& out, isoform_anno& iso)
{
	out<<iso.g_name<<"\t"<<iso.name<<"\t"<<iso.chrom<<"\t"<<iso.strand<<"\t"<<iso.txStart<<"\t"<<iso.txEnd<<"\t"<<iso.cdsStart<<"\t"<<iso.cdsEnd<<"\t"<<iso.exonCount<<"\t";;
	
	for(int i = 0; i < iso.exonCount; i++){
		out<<iso.exonStarts[i]<<",";
	}
	out<<"\t";
	
	for(int i = 0; i < iso.exonCount; i++){
		out<<iso.exonEnds[i]<<",";
	}
	out<<"\n";
}

void get_annotation(map<string,gene_anno>& map_g_anno, ifstream& in_annotation, int annotype)
{
	map<string,gene_anno>::iterator iter_map_g_anno;
    map<string,list<isoform_anno> > isoform_anno_map;
    map<string,list<isoform_anno> >::iterator iter_gene_anno_map;

	string temp_line;

	if(annotype == 1) // deal with refflat format annotation file.
	{
		while(getline(in_annotation, temp_line)){
			isoform_anno g(1,temp_line); //use refflat format annotation file.
			isoform_anno_map[g.g_name].push_back(g);
		}
	}
	else if(annotype == 2) // deal with GTF format annotation file. only deal with the exon annotation. CDS and start/end_codon is ignored.
	{
		map<string, isoform_anno> dealt_trans_anno_map; // key: transcript name, value: transcript annotation.
		map<string, isoform_anno>::iterator trans_anno_iter; // key: transcript name, value: transcript annotation.
		while(getline(in_annotation, temp_line)){
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
	
	
	for(iter_gene_anno_map = isoform_anno_map.begin(); iter_gene_anno_map != isoform_anno_map.end(); iter_gene_anno_map++){
		map_g_anno[(*iter_gene_anno_map).first] = gene_anno( isoform_anno_map[(*iter_gene_anno_map).first] );
    }

	iter_map_g_anno=map_g_anno.begin();
	while(iter_map_g_anno!=map_g_anno.end()){
        if(!if_gene_valid( (*iter_map_g_anno).second) ){
            map_g_anno.erase(iter_map_g_anno++);
        }
		else{
			iter_map_g_anno++;
		}
	}
}

void output_nurd_file(ofstream& out_nurd, const vector<double>& GBC, const map<string,gene_anno>& map_g_anno, const map<string,vector<int> >& gene_read_count, int total_valid_read_count){
	// output GBC curve
	for(int i = 0; i < GBC.size(); i++){
		out_nurd << GBC[i] << "\t";
	}
	out_nurd << "\n";

	// output valid read number and gene number
	out_nurd << total_valid_read_count << "\t";
	out_nurd << map_g_anno.size() << "\n";

	// output each gene's detail information
	map<string,gene_anno>::const_iterator iter_map_g_anno; // const iterator
    for(iter_map_g_anno = map_g_anno.begin(); iter_map_g_anno != map_g_anno.end(); iter_map_g_anno++){
        string tmp_geneName = (*iter_map_g_anno).second.g_name;

        //output as "nurd" format
        //first line: basic information
        out_nurd << tmp_geneName << "\t";
        out_nurd << (*iter_map_g_anno).second.exon_len.size() << "\t";
        out_nurd << (*iter_map_g_anno).second.iso_name.size() << "\n";

        //second line: isoform name
        for(int i = 0; i < (*iter_map_g_anno).second.iso_name.size(); i++){
            out_nurd << (*iter_map_g_anno).second.iso_name[i]<<"\t";
        }
        out_nurd << "\n";

        //third line: exon length
        const vector<_chr_coor>& tmp_exonLen = (*iter_map_g_anno).second.exon_len;
        int tmp_exon_num = tmp_exonLen.size();
        for(int i = 0; i < tmp_exon_num; i++){
            out_nurd << tmp_exonLen[i] << "\t";
        }
        out_nurd << "\n";
        //forth line: read count in exon
		map<string,vector<int> >::const_iterator iter_map_g_rdcnt = gene_read_count.find(tmp_geneName); // rdcnt: short for read count.
		if(iter_map_g_rdcnt != gene_read_count.end()){
			for(int i = 0; i < tmp_exon_num; i++){
				out_nurd << iter_map_g_rdcnt->second[i] << "\t";
			}
		}
        out_nurd << "\n";
        //other line: gene structure
        for(int i = 0; i < (*iter_map_g_anno).second.iso_name.size(); i++){
            for(int j = 0; j < tmp_exon_num; j++){
                out_nurd << (*iter_map_g_anno).second.exon_iso_idx[i][j] << "\t";
            }
            out_nurd << "\n";
        }
    }
}

int get_exon_read_count(ifstream& in_annotation, ifstream& in_sam, ofstream& out_nurd, bool if_chr, int annotype){
    clock_t start_time,end_time;
    time_t cur_time;
    start_time = clock();

	clock_t tmp_start, tmp_end;
	stringstream ss (stringstream::in | stringstream::out);
		
    map<string,gene_anno> map_g_anno;
    map<string,gene_anno>::iterator iter_map_g_anno;

	get_annotation(map_g_anno, in_annotation, annotype);
	
	end_time = clock();
    ss<<"gene annotation time: "<<((double)end_time-start_time)/CLOCKS_PER_SEC<<" seconds.\n";
    std_output_with_time(ss.str());
    ss.str("");
    start_time=end_time;
	
	ss<<"gene annotation done!"<<endl;
    std_output_with_time(ss.str());
    ss.str("");
	
    //gene read count begins
	////////////////////////////////////////
	/////// read count starts!
	map<string,vector<int> > gene_read_count;
	for(iter_map_g_anno = map_g_anno.begin(); iter_map_g_anno != map_g_anno.end(); iter_map_g_anno++){
		gene_read_count[(*iter_map_g_anno).first] = vector<int>((*iter_map_g_anno).second.exon_num);
	}

	map<string, list<_chr_coor> > gene_read_pos_map;
	for(iter_map_g_anno = map_g_anno.begin(); iter_map_g_anno != map_g_anno.end(); iter_map_g_anno++){
		gene_read_pos_map[(*iter_map_g_anno).first] = list<_chr_coor>();
	}

	int total_read_count = 0;
	int total_valid_read_count = 0;

	// bin number in GBC. GBC is calculated at the same time of reads counting.
	int num_bins = 10;
	vector<double> GBC(num_bins);
	for(int i = 0; i < num_bins; i++){
		GBC[i] = 0.0;
	}
	int outlier_read_cnt = 0;

	//// first key is chr name, second key is gene start pos, second value is gene name. The first value is a map container.
	//// This is map nest definition.
	map<string, map<int,string> > map_chr_pos_gene;
	map<string, map<int,string> >::iterator iter_map_chr_pos_gene; // to help judge if some element exists

	// the key is chr name, the value is a vector, in which are the gene start positions on this chr.
	// no gene name. Gene names can be easily got from the former hash map.
	// this vector should be sorted before binary search.
	map<string, vector<int> > map_chr_start_pos_vec;
	map<string, vector<int> > map_chr_end_pos_vec;
	map<string, vector<int> >::iterator iter_map_chr_pos_vec;

	int start_pos, end_pos;
	for(iter_map_g_anno = map_g_anno.begin(); iter_map_g_anno != map_g_anno.end(); iter_map_g_anno++)
	{
		//frist map
		iter_map_chr_pos_gene = map_chr_pos_gene.find((*iter_map_g_anno).second.chrom);
		if(iter_map_chr_pos_gene == map_chr_pos_gene.end()){
			map<int,string> tmp_map;
			start_pos = (*iter_map_g_anno).second.g_start;
			tmp_map[start_pos] = (*iter_map_g_anno).first;

			map_chr_pos_gene[(*iter_map_g_anno).second.chrom] = tmp_map;
		}
		else{
			start_pos = (*iter_map_g_anno).second.g_start;
			map_chr_pos_gene[(*iter_map_g_anno).second.chrom][start_pos] = (*iter_map_g_anno).first;
		}

		//second map
		iter_map_chr_pos_vec = map_chr_start_pos_vec.find((*iter_map_g_anno).second.chrom);
		if(iter_map_chr_pos_vec == map_chr_start_pos_vec.end()){
			vector<int> tmp_vec;
			start_pos = (*iter_map_g_anno).second.g_start;
			tmp_vec.push_back(start_pos);

			map_chr_start_pos_vec[(*iter_map_g_anno).second.chrom] = tmp_vec;
		}
		else{
			start_pos = (*iter_map_g_anno).second.g_start;
			map_chr_start_pos_vec[(*iter_map_g_anno).second.chrom].push_back(start_pos);
		}
	}

	//sort the vector before binary search
	for(iter_map_chr_pos_vec = map_chr_start_pos_vec.begin(); iter_map_chr_pos_vec != map_chr_start_pos_vec.end(); iter_map_chr_pos_vec++){
		sort((*iter_map_chr_pos_vec).second.begin(),(*iter_map_chr_pos_vec).second.end());
	}

	// get map_chr_end_pos_vec
	// because the start pos is sorted, so the correspond end pos can't be sure sorted. So the corresponding vector should be got based on the start pos
	for(iter_map_chr_pos_vec = map_chr_start_pos_vec.begin(); iter_map_chr_pos_vec != map_chr_start_pos_vec.end(); iter_map_chr_pos_vec++){
		int chr_s_pos_size = (*iter_map_chr_pos_vec).second.size();
		vector<int> tmp_ends(chr_s_pos_size);
		for(int i = 0; i < chr_s_pos_size; i++){
			string tmp_geneName = map_chr_pos_gene[(*iter_map_chr_pos_vec).first][ (*iter_map_chr_pos_vec).second[i] ];
			tmp_ends[i] = map_g_anno[tmp_geneName].g_end;
		}
		map_chr_end_pos_vec[(*iter_map_chr_pos_vec).first] = tmp_ends;
	}

	string temp_line;
	
	//deal with header
	in_sam.seekg(0,ios::beg);
	while(getline(in_sam,temp_line)){
		if(temp_line[0] != '@'){
			break;
		}
	}
	//deal with mapped reads
	int read_Flag_mask_reverse = 0x10;
	// vector<string> sam_column = vector<string>(4);
	vector<string> sam_column = vector<string>(10);
	
	do{
		total_read_count++;
		int invalid_type = 0;//0: valid, 1: not gene region, 2: not exon region 3: multi gene 4: invalid chrome

		//extract the map information from sam file
		// delimiter(sam_column,temp_line,'\t',4,true); // because only the first 4 fields are used. reduce the time that was
		delimiter(sam_column, temp_line, '\t', 10, true); // because only the first 10 fields are used. reduce the time that was
		string chrName = sam_column[2];

		int read_Flag = atoi(sam_column[1].c_str());

		string read_name = sam_column[0];

		if(chrName == "*"){
			continue;
		}
		else{
			iter_map_chr_pos_vec = map_chr_start_pos_vec.find(chrName);
			if(iter_map_chr_pos_vec == map_chr_start_pos_vec.end()){
				continue;
			}
			
			_chr_coor read_pos = atoi(sam_column[3].c_str());
			if( (read_Flag & read_Flag_mask_reverse) !=0 ){ // reverse read
				read_pos += sam_column[9].size();
			}
			
			list<int> l_gene_index = my_bin_search_multi(map_chr_start_pos_vec[chrName],map_chr_end_pos_vec[chrName],read_pos);

			if(l_gene_index.size() == 0){ // if no gene cover the read, deal with the next read
				continue;
			}
			else{
				list<int>::iterator tmp_iter_list;
				bool if_map_to_gene = false;
				bool if_multi_map = false;
				int exon_idx = -1;
				string g_name;
				for(tmp_iter_list = l_gene_index.begin(); tmp_iter_list != l_gene_index.end(); tmp_iter_list++){
					g_name = map_chr_pos_gene[chrName][ map_chr_start_pos_vec[chrName][(*tmp_iter_list)] ];

					if(map_g_anno[g_name].strand == "+"){
						exon_idx = my_bin_search(map_g_anno[g_name].exon_start_g, map_g_anno[g_name].exon_end_g, read_pos);
					}
					///// if on the negtive strand, should use the reverse search
					else{
						exon_idx = my_bin_search_reverse(map_g_anno[g_name].exon_start_g, map_g_anno[g_name].exon_end_g, read_pos);
					}

					if(exon_idx!=-1){
						if(!if_map_to_gene){
							if_map_to_gene = true;
						}
						else{
							if_multi_map = true;
							break;
						}
					}
				}
				if(!if_multi_map){
					if(exon_idx != -1){
						total_valid_read_count++;
						int gene_read_pos = -1; // initialized as an invalid position.
						if(map_g_anno[g_name].strand == "+"){
							gene_read_count[g_name][exon_idx]++;

							//for GBC, read pos
							gene_read_pos = map_g_anno[g_name].exon_g_start_l[exon_idx] + read_pos - map_g_anno[g_name].exon_start_g[exon_idx];
						}
						else{
							gene_read_count[g_name][exon_idx]++;

							//for GBC, read pos
							// still exon_g_start_l, if you figure out, it's obvious.
							gene_read_pos = map_g_anno[g_name].exon_g_start_l[exon_idx] + map_g_anno[g_name].exon_end_g[exon_idx] - read_pos;
							// gene_read_pos_map[g_name].push_back(gene_read_pos);
						}

						//For GBC
						//only use the genes with single isoform
						if(map_g_anno[g_name].iso_name.size() == 1)
						{
							int gene_len = map_g_anno[g_name].g_len;
							if( gene_read_pos < gene_len && gene_read_pos >= 0 ){
								GBC[ (int)(gene_read_pos*num_bins)/gene_len ]++;
							}
							else if(gene_read_pos == gene_len){
								GBC[ num_bins-1 ]++;
							}
							else{
								outlier_read_cnt++;
							}
						}
					}
				}
			}
			continue;
		}

	}while(getline(in_sam,temp_line));	
	//////////  read count ends!
	//////////////////////////////////////////

	//////////////////////////////////
	//////// get GBC
	double total_GBC = 0;
	for(int i = 0; i < num_bins; i++){
		total_GBC += GBC[i];
	}
	for(int i = 0; i < num_bins; i++){
		GBC[i] = ((double)GBC[i]*num_bins)/total_GBC;
	}
	//////////////////////////////////

	end_time = clock();
    ss << "read count time: " << ((double)end_time-start_time)/CLOCKS_PER_SEC << " seconds.\n";
    std_output_with_time(ss.str());
    ss.str("");
    start_time = end_time;

	//output the information to nurd file.
	output_nurd_file(out_nurd, GBC, map_g_anno, gene_read_count, total_valid_read_count);
	
	return 0;
}
