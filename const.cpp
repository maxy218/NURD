/*
 * =====================================================================================
 *
 *       Filename:  const.cpp
 *
 *    Description:  the definitions of some consts.
 *
 *        Version:  1.0
 *        Created:  02/25/2014 23:48:09 PM
 *       Revision:  none
 *       Compiler:  g++
 *
 *         Author:  Xinyun Ma
 *   Organization:  Tsinghua University
 *
 * =====================================================================================
 */


#include <vector>
#include <string>
#include "const.h"
using namespace std;

const double EPSILON = 1e-10;
const double maxFloat = 1e100;
const int max_iteration = (int)1e5;

int bin_N = 10; // bin number got from nurd file. not in exon read counts part.
vector<double> bin;

const int GBC_bin_num = 10;

// typedef unsigned int _chr_coor;
typedef int _chr_coor;
