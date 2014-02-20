#ifndef MYDATA_H_INCLUDED
#define MYDATA_H_INCLUDED

#include <fstream>
#include <string>
#include <vector>
#include <time.h>
#include <math.h>

#include "class.h"
#include "someGlobal.h"
#include "algorithm.h"
using namespace std;

double get_log_likelihood(const gene& ps);

double get_gradient_of_log_likelihood(const gene& ps, int i);

//return the first line of urd file
void get_GBC_bin(ifstream& infile);

void get_LBC_curve(gene& g, double* LBC);

void get_LBC_matrix(gene& g);

void get_GBC_matrix(gene& g);

//N:number of bin
//exon_N:number of exon
void get_Curve_from_bin(const vector<double>& bin, int N, const vector<int>& length, int exon_N, vector<double>& Area);

#endif // MYDATA_H_INCLUDED
