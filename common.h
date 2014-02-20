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

#include <ctime>
#include <cstdio>
#include <string>
#include <iostream>

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

#endif // COMMON_H_INCLUDED 
