#include <vector>
#include <string>
#include "someGlobal.h"
using namespace std;

const double epsilon = 1e-10;
const double maxFloat = 1e100;
const int max_iteration = (int)1e5;

int bin_N; // bin number got from nurd file. not in exon read counts part.
vector<double> bin;

// typedef unsigned int _chr_coor;
typedef int _chr_coor;
