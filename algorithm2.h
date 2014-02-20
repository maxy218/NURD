#ifndef ALGORITHM2_H_INCLUDED
#define ALGORITHM2_H_INCLUDED

#include <fstream>
#include <vector>
#include <iostream>
#include <sstream>
#include <math.h>

#include "class.h"
#include "myData.h"
#include "someGlobal.h"
#include "common.h"

using namespace std;

void get_anno_info(ifstream& in_anno, const int anno_choice, map<string, gene_anno> & map_g_anno);

#endif // ALGORITHM_H_INCLUDED
