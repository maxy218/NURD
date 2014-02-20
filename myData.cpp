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

double get_log_likelihood(const gene& ps){
    int N = ps.N;
    int M = ps.M;

    double likeli = 0.0;
	
	for(int j = 0; j < N; j++){
		double tempSum1 = 0.0;
		double tempSum2 = 0.0;
		for(int i = 0; i < M; i++){
			double tmp_double = 0.0;
			tmp_double = ps.c[i*N+j]*ps.theta[i];
			tempSum1 += ps.l[j]*tmp_double;
			tempSum2 += tmp_double;
		}
		likeli = likeli-ps.w*tempSum1+ps.x[j]*log(ps.l[j]*ps.w*tempSum2+epsilon);
	}

    return likeli;
}

double get_gradient_of_log_likelihood(const gene& ps, int i){
    int N = ps.N;
    int M = ps.M;

    double gradient = 0.0;

    for(int j = 0; j < N; j++){
		gradient -= ps.w*ps.l[j]*ps.c[i*N+j];
		if(ps.x[j]*ps.c[i*N+j] != 0)
		{
			double tmp_sum = epsilon;
			for(int k = 0; k < M; k++)
			{
				tmp_sum += ps.c[k*N+j]*ps.theta[k];
			}
			
			gradient += ps.x[j]*ps.c[i*N+j]/tmp_sum;
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

void get_LBC_curve(gene& g, vector<double>& LBC){
    int M = g.M;
    int N = g.N;

    double tempSum = 0.0;
    double LBC_sum = 0.0;
    for(int j = 0; j < N; j++){
        tempSum = 0.0;
        for(int i = 0;i<M;i++){
            tempSum += g.theta[i]*g.a[i*N+j];
        }
        if(tempSum < epsilon){
            tempSum += epsilon;
        }
        LBC[j] = g.x[j]/(g.l[j]*tempSum);
        LBC_sum += LBC[j];
    }

	for(int j = 0; j < N; j++){
		LBC[j] = LBC[j]*N/LBC_sum;
	}
}

void get_LBC_matrix(gene& g){
    int M = g.M;
    int N = g.N;

    vector<double> LBC = vector<double>(N);
    get_LBC_curve(g, LBC);

    vector<int> newLength = vector<int>(N);
    int newTotalLength = 0;

    for(int j = 0; j < N; j++){
        newLength[j] = 0;
        for(int i = 0; i < M; i++){
            newLength[j] += (int)g.a[i*N+j];
        }
        newLength[j] *= g.l[j];
        newTotalLength += newLength[j];
    }

    vector<double> newLBC = vector<double>(newTotalLength);

    int curPos = 0;
    for(int j = 0; j < N; j++){
        for(int subj = 0; subj < newLength[j]; subj++){
            newLBC[curPos+subj]=LBC[j];
        }
        curPos += newLength[j];
    }

    double tempSum = 0.0;
    for(int j = 0; j < newTotalLength; j++){
        tempSum += newLBC[j];
    }
    for(int j = 0; j < newTotalLength; j++){
        newLBC[j] = newLBC[j]*newTotalLength/tempSum;
    }

    vector<double> temp = vector<double>(N);
    vector<int> temp_length = vector<int>(N);

    int gene_length = 0;
    for(int j = 0; j < N; j++){
        gene_length += g.l[j];
    }

    for(int i = 0; i < M; i++){
        for(int j = 0; j < N; j++){
            temp_length[j] = g.l[j]*(int)g.a[i*N+j];
        }

        get_Curve_from_bin(newLBC, newTotalLength, temp_length, N, temp);
		for(int j = 0; j < N; j++){
			g.LBC[i*N+j] = temp[j]*bin_N/newTotalLength;
		}
    }
}

void get_GBC_matrix(gene& g){
    int M = g.M;
    int N = g.N;
    vector<double> temp = vector<double>(N);
    vector<int> temp_length = vector<int>(N);

    for(int i = 0; i < M; i++){
		for(int j = 0; j < N; j++){
			temp_length[j] = g.l[j]*(int)g.a[i*N+j];
		}
		get_Curve_from_bin(bin, bin_N, temp_length, N, temp);
		for(int j = 0; j < N; j++){
			g.GBC[i*N+j] = temp[j];
		}
    }
}

//N:number of bin
//exon_N:number of exon
void get_Curve_from_bin(const vector<double>& bin, int N, const vector<int>& length, int exon_N, vector<double>& Area){
    vector<double> new_Length = vector<double>(exon_N);

    int total_Length = 0;
    for(int i = 0; i < exon_N; i++){
        total_Length += length[i];
    }
    for(int i = 0; i < exon_N; i++){
        new_Length[i] = (double)length[i]/total_Length*N;
    }

    double left_coor = 0.0;
    double right_coor = new_Length[0];

    double temp_pos = left_coor;
    double delta = 0;
    for(int i = 0; i < exon_N; i++){
        Area[i] = 0.0;

        while(right_coor > temp_pos){
            delta = min(((int)(temp_pos+1)-temp_pos),right_coor-temp_pos);
            Area[i] += delta*bin[(int)temp_pos];
            temp_pos += delta;
        }

        if(new_Length[i] == 0){
            Area[i] = 0;
        }

        if(i < exon_N-1){
            left_coor += new_Length[i];
            right_coor += new_Length[i+1];
        }
    }
}
