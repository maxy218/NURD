/*
 * =====================================================================================
 *
 *       Filename:  main.cpp
 *
 *    Description:  the framework of NURD. 
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


#include <algorithm>
#include <cerrno>   // error information
#include <cstdlib>
#include <cstring>  // strerror
#include <ctime>
#include <fstream>
#include <iostream>
#include <list>
#include <map>
#include <sstream>
#include <stdexcept>
#include <vector>

#include <dirent.h> // directory operation
#include <unistd.h> // parsing argument
#include <sys/stat.h>   // mkdir and access

#include "class.h"
#include "common.h"
#include "algorithm.h"

using namespace std;

void usage()
{
  cout << "NURD is a tool to estimate isoform expression with RNA-Seq data." << endl;
  cout << "Version: 1.0.8" << endl;
  cout << "========================================" << endl;
  cout << "Usage:\t" << "NURD [options] <-G annotation.gtf>|<-R annotation.refflat> <-S mapping_file.sam>" << endl;
  cout << "\t" << "-G: annotation.gtf: gene annotation in gtf format." << endl;
  cout << "\t" << "-R: annotation.refflat: gene annotation in refflat format." << endl;
  cout << "\t" << "-S: mapping_result.sam: reads mapping result in sam format." << endl;
  cout << "========================================" << endl;
  cout << "\t" << "Optional:" << endl;
  cout << "\t" << "-A: the weight of GBC when mixturing the GBC and LBC into one gene structure matrix. It's a float number between 0-1. Default: 0.5" << endl;
  cout << "\t" << "-O: output_dir: the directory to output the estimation result. Default: current directory" << endl;
  cout << "========================================" << endl;
}

// get the name of sam file.
// take both linux and windows into consider.
string get_file_name(const string& dir, int& flag){
  flag = 0;
  if(dir.size() > 0){
    if( dir[dir.size() - 1] == '\\' || dir[dir.size() - 1] == '/' ){
      flag = 1;
      cerr << "fatal error!" << endl; // should throw some exception.
      return "";
    }
    size_t length = 0;
    int i = dir.size() - 1;
    for(; i >= 0 ; i--){
      if( dir[i] == '\\' || dir[i] == '/' ){
        break;
      }
      length++;
    }
    return dir.substr(i + 1 , length);
  }
  else{
    flag = 2;
    return "";
  }
}

// create directory.
int create_dir(const string&  dirName)
{
  string dir = dirName + "/";

  for( int i=1; i < dir.length(); i++)
  {
    if(dir[i] == '/')
    {
      dir[i] = 0;
      if( access(dir.c_str(), F_OK) != 0 )
      {
      if(mkdir(dir.c_str(), 0755) == -1) // mode 755 means: rwxr-xr-x
        {
          return -1;
        }
      }
      dir[i] = '/';
    }
  }
  return 0;
}

int main(int argc, char**argv)
{
  // parsing argument
  map<string,string> argu_parse_result = map<string,string>();
  map<string,string>::iterator argu_parse_iterator;

  ifstream in_anno; // annotation file.
  ifstream in_rdmap; // reads mapping file.
  ofstream out_nurd; // tmp file. nurd file.
  ifstream in_nurd; // tmp file. nurd file.
  ofstream out_expr; // expression estimation file.

  string in_anno_name;
  string in_rdmap_name;
  string out_nurd_name;
  string out_expr_name;

  int anno_choice = 2; // choice of annotation: 1: refflat, 2: GTF
  double alpha = 0.5; // weight of GBC and LBC. 0.5 as default.

  try
  {
    int ch;
    opterr = 0;
    while((ch = getopt(argc,argv,"O:R:G:S:"))!= -1)
    {
      switch(ch)
      {
        case 'O': argu_parse_result["O"] = optarg; break;
        case 'R': argu_parse_result["R"] = optarg; break;
        case 'G': argu_parse_result["G"] = optarg; break;
        case 'S': argu_parse_result["S"] = optarg; break;
        case 'A': argu_parse_result["A"] = optarg; break;
      }
    }

    //post process of argument.
    //situation 1: -O is not specified.
    if( (argu_parse_iterator = argu_parse_result.find("O")) == argu_parse_result.end() )
    {
      argu_parse_result["O"] = "./";
    }
    else
    {
      argu_parse_result["O"] += "/";
      if( create_dir(argu_parse_result["O"]) == -1 )
      {
        throw runtime_error("mkdir error!\t");
      }
    }

    //situation 2: neithor of  -G or -R is specified: error
    if( (argu_parse_result.find("R") == argu_parse_result.end()) && (argu_parse_result.find("G") == argu_parse_result.end()) )
    {
      throw runtime_error("no annotation is specified.");
    }

    //situation 3: both of  -G and -R are specified: error
    if( (argu_parse_result.find("R") != argu_parse_result.end()) && (argu_parse_result.find("G") != argu_parse_result.end()) )
    {
      throw runtime_error("there are multiple annotation formats.");
    }
    //situation 4: -R is specified, but no refflat file.
    if( argu_parse_result.find("R") != argu_parse_result.end() )
    {
      in_anno_name = argu_parse_result["R"];
      in_anno.open(in_anno_name.c_str());
      if(!in_anno.is_open()){
        throw runtime_error("can not open refflat file! Please check whether there exists the refflat file you specified.");
      }
      anno_choice = 1;
    }
    //situation 5: -G is specified, but no gtf file.
    if( argu_parse_result.find("G") != argu_parse_result.end() )
    {
      in_anno_name = argu_parse_result["G"];
      in_anno.open(in_anno_name.c_str());
      if(!in_anno.is_open()){
        throw runtime_error("can not open gtf file! Please check whether there exists the gtf file you specified.");
      }
      anno_choice = 0;
    }

    //situation 6: -S is not specified: error
    if( argu_parse_result.find("S") == argu_parse_result.end() )
    {
      throw runtime_error("no mapping file is specified.");
    }
    else
    {
      in_rdmap_name = argu_parse_result["S"];
      in_rdmap.open(in_rdmap_name.c_str());
      if(!in_rdmap.is_open()){
        throw runtime_error("can not open sam file! Please check whether there exists the sam file you specified.");
      }
    }
    
    // situation 7: -A is out of range.
    if( argu_parse_result.find("A") == argu_parse_result.end() )
    {
      argu_parse_result["A"] = "0.5"; // default: 0.5
    }
    else{
      alpha = atof(argu_parse_result["A"].c_str());
      if(alpha > 1.0 || alpha < 0.0){
        throw runtime_error("illegal alpha. (alpha should be in [0, 1])");
      }
    }
  }
  catch(runtime_error err)
  {
    cerr << "Exception catched:" << "\t";
    cerr << err.what() << endl << endl;
    usage();
    return 1;
  }

  stringstream ss (stringstream::in | stringstream::out);
  std_output_with_time(string("expression estimation start!\n"));

  clock_t start_time, end_time;
  start_time=clock();

  //get the file names.
  int flag = 0;
  string in_rdmap_name_no_dir = get_file_name(in_rdmap_name, flag);
  if(flag != 0){
    cerr << "Invalid sam file name!" << endl;
    return 1;
  }

  out_nurd_name = argu_parse_result["O"] + in_rdmap_name_no_dir + ".nurd";
  std_output_with_time("sam file:\t" + out_nurd_name + "\n");
  out_expr_name = out_nurd_name + ".all_expr";
  std_output_with_time("expression file:\t" + out_expr_name + "\n");


  try{
    out_nurd.open(out_nurd_name.c_str());
    if(!out_nurd.is_open()){
      throw runtime_error("can not open nurd file! Please check whether there exists the nurd file.");
    }
  }
  catch(runtime_error err){
    cerr << "Exception catched:" << "\t";
    cerr << err.what() << endl << endl;
    return 1;
  }

  map<string, gene_info> map_g_info;
  get_anno_info(in_anno, anno_choice, map_g_info);

  cout << "number of genes with anno:" << map_g_info.size() << endl;

  vector<double> GBC = vector<double>(GBC_bin_num, 0.0);

  size_t tot_valid_rd_cnt = 0;
  get_exon_rd_cnt(map_g_info, in_rdmap, tot_valid_rd_cnt, GBC);




        out_nurd.close();

        //expression estimation
        try{
                in_nurd.open(out_nurd_name.c_str());
                if(! in_nurd){
                        throw runtime_error("can not open nurd file! Please check whether there exists the nurd file.");
                }
        }
        catch(runtime_error err){
                cerr<<"Exception catched:"<<"\t";
                cerr<<err.what()<<endl<<endl;
                return 1;
        }


        out_expr.open(out_expr_name.c_str());

        time_t esti_start, esti_end;
        esti_start = clock();
        calcuAllTheGenes(map_g_info, tot_valid_rd_cnt, alpha, GBC, out_expr);
        esti_end = clock();

        ss<<"expression estimation time:\t"<<((double)esti_end-esti_start)/CLOCKS_PER_SEC<<" seconds.\n";
    std_output_with_time(ss.str());
    ss.str("");

    ss<<"expression done!"<<endl;
    std_output_with_time(ss.str());
    ss.str("");

        end_time = clock();
    ss << "total time:\t" << ((double)end_time-start_time)/CLOCKS_PER_SEC << " seconds." <<endl;
    std_output_with_time(ss.str());
    ss.str("");

    return 0;




















}
