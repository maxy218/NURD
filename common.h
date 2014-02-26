/*
 * =====================================================================================
 *
 *       Filename:  common.h
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


#ifndef COMMON_H_INCLUDED
#define COMMON_H_INCLUDED

#include <algorithm>
#include <ctime>
#include <cstdio>
#include <string>
#include <vector>
#include <iostream>
#include <fstream>

#include <unistd.h> // parsing argument
#include <sys/stat.h> 	// mkdir and access 
#include <cerrno> 	// error information
#include <cstring>	// strerror

using namespace std;

bool std_output_with_time(string s);

vector<string> delimiter(string str, char deli);

// if user know that there are N fields, then the vector can be initialized by size N. Maybe could speed up.
// if there are more than N fields, it will return the first N fields
////////        can add a parameter to the function, such as flag which is boolean variable, can tell whether there are more than N fields.
vector<string> delimiter(string str, char deli, int N);

////// return the first N fields, if there are more than N fields. If less than N, the extra field are filled by Null value.
vector<string> delimiter(string str, char deli, int N, bool if_head_N);

////// return the first N fields, if there are more than N fields. If less than N, the extra field are filled by Null value.
////// use the reference to reduce the time consuming.
void delimiter(vector<string>& t, string str, char deli, int N, bool if_head_N);

template<typename T>
void output_vector(const vector<T> & vec, ostream& out, const char deli){
  const size_t sz = vec.size();
  for(size_t idx = 0; idx < sz; ++idx){
    out << vec[idx] << deli;
  }
  out << endl;
}

template<typename T>
void output_2D_vector(const vector< vector<T> > & vec, ostream& out, const char deli){
  const size_t sz1 = vec.size();
  if(sz1 == 0){
    return;
  }
  const size_t sz2 = vec[0].size();
  for(size_t idx1 = 0; idx1 < sz1; ++idx1){
    for(size_t idx2 = 0; idx2 < sz2; ++idx2){
      out << vec[idx1][idx2] << deli;
    }
    out << endl;
  }
  out << endl;
}

template<typename T>
T sum_vector(const vector<T> & vec){
  T sum = 0;
  const size_t sz = vec.size();
  for(size_t idx = 0; idx < sz; ++idx){
    sum += vec[idx];
  }
  return sum;
}

//return -1 if fail
//else, return the index of x
//vec is a vector of boundary of each interval
//  n+1 elements represent n intervals, the first elem is 0
//the interval is left close and right open
// [a, b)
template <class T >
int bin_search(const vector<T> & vec, const T & x){
  // invalid boundaries.
  if(vec.size() < 2){
    return -1;
  }

  int left = 0;
  int right = vec.size()-2;
  int mid;//find the middle interval, size-1 for there are n+1 elems in vector

  if(x >= vec[right + 1] || x < vec[left]){
    return -1;
  }
  while(left <= right){
    mid = left + (right - left)/2;
    if(x >= vec[mid] && x < vec[mid+1]){
      return mid;
    }
    else if(x < vec[mid]){
      right = mid - 1;
    }
    else{
      left = mid + 1;
    }
  }
  return left;
}

// two vector version. It's totally different with the one vector version.
// one vector version can easily transform into two vector version.
// vec1: starts     vec2:ends
// reference: introduction to the Design and analysis of algorithms(second edition)
//  related chapter: chapter4. Chinese version, P104
template <class T >
int bin_search(const vector<T>& vec1, const vector<T>& vec2, const T & x){
  int left = 0;
  int right = vec1.size() - 1;
  int mid;

  //if x is out of bound, return -1. false
  if(x >= vec2[right] || x < vec1[left]){
    return -1;
  }
  while(left <= right){
    mid = left + (right - left)/2;
    if(x >= vec1[mid] && x < vec2[mid]){
      return mid;
    }
    else if(x < vec1[mid]){
      right = mid - 1;
    }
    else{
      left = mid + 1;
    }
  }
  return -1;
}

// This version is for the array that the element is sorted from large element to small element. The reverse of above
template <class T >
int bin_search_reverse(const vector<T>& vec1, const vector<T>& vec2, const T& x){
  int left = 0;
  int right = vec1.size() - 1;
  int mid;

  //if x is out of bound, return -1. false
  if(x >= vec2[left] || x < vec1[right]){
    return -1;
  }
  while(left <= right){
    mid = left + (right - left)/2;
    if(x >= vec1[mid] && x < vec2[mid]){
      return mid;
    }
    else if(x < vec1[mid]){
      left = mid + 1;
    }
    else{
      right = mid - 1;
    }
  }
  return -1;
}

// two vector and multiple return value version.
// It's similar with the above binary search version
// If there are multiple hit, return a list of recode. The list record the index of two vector.
// If there's no hit, return the null list, whose length is 0.
template <class T >
vector<int> bin_search_multi(const vector<T>& vec1, const vector<T>& vec2, const T & x){
  int left = 0;
  int right = vec1.size()-1;
  int mid;
  vector<int> result;

  //if x is out of bound, return -1. false
  if(x >= vec2[right] || x < vec1[left]){
    return result;
  }
  while(left <= right){
    mid = left + (right - left)/2;
    if(x >= vec1[mid] && x < vec2[mid]){
      result.push_back(mid);
      break;
    }
    else if(x < vec1[mid]){
      right = mid - 1;
    }
    else{
      left = mid + 1;
    }
  }

  int mid_bak = mid;
  int flank_gene = 10; // allowing for flank 10 genes to search the covered genes.

  // search the region before m. Stop when the x is larger than the right bound
  // not perfect. Because only starts are sorted, so it can't guarantee to find all the proper interval.
  mid--;
  int left_flank = 0;
  while(mid >= 0 && left_flank < flank_gene){
    left_flank++;
    if(x >= vec1[mid] && x < vec2[mid]){
      result.push_back(mid);
    }
    mid--;
  }
  // search the region after m. Stop when the x is smaller than the left bound
  // This half should be perfect.
  mid = mid_bak + 1;
  int total_size = vec1.size();
  int right_flank = 0;
  while(mid < total_size && right_flank < flank_gene){
    right_flank++;
    if(x >= vec1[mid] && x < vec2[mid]){
      result.push_back(mid);
    }
    mid++;
  }
  sort(result.begin(), result.end());
  return result;
}

#endif // COMMON_H_INCLUDED 
