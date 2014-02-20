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


#include <iostream>
#include <vector>
#include <list>
#include <map>
#include <fstream>
#include <cstdlib>
#include <algorithm>
#include <sstream>

#include <unistd.h> // parsing argument
#include <sys/stat.h>   // mkdir and access
#include <dirent.h> // directory operation
#include <cerrno>   // error information
#include <cstring>  // strerror
#include <stdexcept>
#include <ctime>

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
  ofstream out_expr; // expression estimation file.

  string in_anno_name;
  string in_rdmap_name;
  string out_nurd_name;
  string out_expr_name;

  int anno_choice = 0; // choice of annotation: 0: gtf, 1: refflat.
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
      in_anno.open(argu_parse_result["R"].c_str());
      if(!in_anno.is_open()){
        throw runtime_error("can not open refflat file! Please check whether there exists the refflat file you specified.");
      }
      anno_choice = 1;
    }
    //situation 5: -G is specified, but no gtf file.
    if( argu_parse_result.find("G") != argu_parse_result.end() )
    {
      in_anno.open(argu_parse_result["G"].c_str());
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
      in_rdmap.open(argu_parse_result["S"].c_str());
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
}
