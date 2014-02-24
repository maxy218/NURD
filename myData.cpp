#include <fstream>
#include <string>
#include <vector>
#include <time.h>
#include <math.h>

#include "myData.h"
#include "class.h"
#include "someGlobal.h"
#include "algorithm.h"
using namespace std;

double get_log_likelihood(const gene_info& g){
  int N = g.exon_num;
  int M = g.iso_num;

  double likeli = 0.0;
  
  for(int j = 0; j < N; j++){
    double tempSum1 = 0.0;
    double tempSum2 = 0.0;
    for(int i = 0; i < M; i++){
      double tmp_double = 0.0;
      tmp_double = g.c[i*N+j] * g.theta[i];
      tempSum1 += g.exon_len[j] * tmp_double;
      tempSum2 += tmp_double;
    }
    likeli = likeli - g.tot_rd_cnt*tempSum1 + g.rd_cnt[j]*log(g.exon_len[j]*g.tot_rd_cnt*tempSum2+EPSILON);
  }

  return likeli;
}

double get_gradient_of_log_likelihood(const gene_info& g, int i){
  int N = g.exon_num;
  int M = g.iso_num;

  double gradient = 0.0;
  for(int j = 0; j < N; j++){
  gradient -= g.tot_rd_cnt*g.exon_len[j]*g.c[i*N+j];
  if(g.rd_cnt[j]*g.c[i*N+j] != 0)
  {
      double tmp_sum = EPSILON;
      for(int k = 0; k < M; k++)
      {
        tmp_sum += g.c[k*N+j]*g.theta[k];
      }
      gradient += g.rd_cnt[j]*g.c[i*N+j]/tmp_sum;
    }
  }
  return gradient;
}

//return the first line of urd file
void get_GBC_bin(ifstream& infile){
  int i = 0;
  double x;
  vector<double> temp_bin = vector<double>(100);
  while(infile.peek() != EOF){
    if(isdigit(infile.peek())){
      infile >> x;
      temp_bin[i] = x;
      i++;
    }
    if(infile.peek() == '\n'){
      break;
    }
    infile.ignore();
  }
  bin_N = i;
  bin = vector<double>(bin_N);
  for(int i = 0; i < bin_N; i++){
    bin[i] = temp_bin[i];
  }
}

void get_LBC_curve(gene_info& g, vector<double>& LBC){
  int M = g.iso_num;
  int N = g.exon_num;

  double tempSum = 0.0;
  double LBC_sum = 0.0;
  for(int j = 0; j < N; j++){
    tempSum = 0.0;
    for(int i = 0; i < M; i++){
      tempSum += g.theta[i]*g.a[i*N+j];
    }
    if(tempSum < EPSILON){
      tempSum += EPSILON;
    }
    LBC[j] = g.rd_cnt[j]/(g.exon_len[j]*tempSum);
    LBC_sum += LBC[j];
  }

  for(int j = 0; j < N; j++){
    LBC[j] = LBC[j]*N/LBC_sum;
  }
}

void get_LBC_matrix(gene_info& g){
  int M = g.iso_num;
  int N = g.exon_num;

  vector<double> LBC = vector<double>(N);
  get_LBC_curve(g, LBC);

  vector<int> new_len = vector<int>(N);
  int new_tot_len = 0;

  for(int j = 0; j < N; j++){
    new_len[j] = 0;
    for(int i = 0; i < M; i++){
      new_len[j] += (int)g.a[i*N+j];
    }
    new_len[j] *= g.exon_len[j];
    new_tot_len += new_len[j];
  }

  vector<double> newLBC = vector<double>(new_tot_len);
  int curPos = 0;
  for(int j = 0; j < N; j++){
    for(int subj = 0; subj < new_len[j]; subj++){
      newLBC[curPos+subj] = LBC[j];
    }
    curPos += new_len[j];
  }

  double tempSum = 0.0;
  for(int j = 0; j < new_tot_len; j++){
    tempSum += newLBC[j];
  }
  for(int j = 0; j < new_tot_len; j++){
    newLBC[j] = newLBC[j]*new_tot_len/tempSum;
  }

  vector<double> temp = vector<double>(N);
  vector<int> temp_length = vector<int>(N);
  int gene_length = 0;
  for(int j = 0; j < N; j++){
    gene_length += g.exon_len[j];
  }

  for(int i = 0; i < M; i++){
    for(int j = 0; j < N; j++){
      temp_length[j] = g.exon_len[j]*(int)g.a[i*N+j];
    }
    get_Curve_from_bin(newLBC, new_tot_len, temp_length, N, temp);
    for(int j = 0; j < N; j++){
      g.LBC[i*N+j] = temp[j]*bin_N/new_tot_len;
    }
  }
}

void get_GBC_matrix(gene_info& g){
  int M = g.iso_num;
  int N = g.exon_num;
  vector<double> temp = vector<double>(N);
  vector<int> temp_length = vector<int>(N);

  for(int i = 0; i < M; i++){
    for(int j = 0; j < N; j++){
      temp_length[j] = g.exon_len[j]*(int)g.a[i*N+j];
    }
    get_Curve_from_bin(bin, bin_N, temp_length, N, temp);
    for(int j = 0; j < N; j++){
      g.GBC[i*N+j] = temp[j];
    }
  }
}

//N:number of bin
//exon_N:number of exon
void get_Curve_from_bin(const vector<double>& bin, int N, const vector<int>& length, 
    int exon_N, vector<double>& Area){
  vector<double> new_len = vector<double>(exon_N);

  int total_Length = 0;
  for(int i = 0; i < exon_N; i++){
    total_Length += length[i];
  }
  for(int i = 0; i < exon_N; i++){
    new_len[i] = (double)length[i]/total_Length*N;
  }

  double left_coor = 0.0;
  double right_coor = new_len[0];

  double temp_pos = left_coor;
  double delta = 0;
  for(int i = 0; i < exon_N; i++){
    Area[i] = 0.0;

    while(right_coor > temp_pos){
      delta = min(((int)(temp_pos+1)-temp_pos), right_coor-temp_pos);
      Area[i] += delta * bin[(int)temp_pos];
      temp_pos += delta;
    }

    if(new_len[i] == 0){
      Area[i] = 0;
    }

    if(i < exon_N-1){
      left_coor += new_len[i];
      right_coor += new_len[i+1];
    }
  }
}
