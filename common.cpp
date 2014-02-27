/*
 * =====================================================================================
 *
 *       Filename:  common.cpp 
 *
 *    Description:  some useful functions
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


#include <cstdio>
#include <ctime>
#include <fstream>
#include <string>
#include <vector>

#include <unistd.h> // parsing argument
#include <sys/stat.h> 	// mkdir and access 

#include "common.h"

using namespace std;

bool std_output_with_time(const string & s)
{
  time_t cur_time;
  time(&cur_time);
  printf("[%15.15s]:  %s", ctime(&cur_time)+4, s.c_str());
}

//delimeter
vector<string> delimiter(const string & str, const char deli){
  vector<string> t;
  size_t l_ind = 0;
  size_t r_ind = 0;
  while(str[r_ind] != '\0'){
    if(str[r_ind] != deli){
      r_ind++;
    }
    else{
      t.push_back(str.substr(l_ind, r_ind - l_ind));
      r_ind++;
      l_ind = r_ind;
    }
  }
  if(r_ind != l_ind){
    t.push_back(str.substr(l_ind, r_ind - l_ind));
  }
  return t;
}

// if user know that there are N fields, then the vector can be initialized by size N. Maybe could speed up.
// if there are more than N fields, it will return the first N fields
// return by reference, to speed up
void delimiter_ret_ref(const string & str, const char deli, const int N, vector<string> & t){
  size_t l_ind = 0;
  size_t r_ind = 0;
  size_t cur_field = 0;
  while(str[r_ind] != '\0'){
    if(str[r_ind] != deli){
      r_ind++;
    }
    else{
      if(cur_field >= N){
         return; 
      }
      t[cur_field++] = str.substr(l_ind,r_ind-l_ind);
      r_ind++;
      l_ind = r_ind;
    }
  }
  if(r_ind != l_ind){
    if(cur_field >= N){
      return;
    }
    t[cur_field++] = str.substr(l_ind,r_ind-l_ind);
  }
  return;
}

