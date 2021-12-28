#include"HKPBC.h"
#include <iostream>
#include <cstdlib>
#include <cmath>
#include <assert.h>
#include <vector>
#include <random>
#include <fstream> 
#include <sstream>
#include<algorithm> 
#include<ctime>  
#include <stdio.h>
#include <stdlib.h>


#define Lx 20
#define Ly 20
#define N L*L //total nodes
#define steps 10*N //steps of the LRW
#define iteration 50 //total iterations
#define M_PI 3.14159265358979323846


using namespace std;


vector<vector<double>> levy(int steps, double C);
vector<vector<int>> split (const std::vector<int>& v); 
vector<int> Divisors(int n);


int main(void)
{
	srand(time(0));
	random_device rd;
	mt19937 gen(rd()); //Mersenne Twister RNG
	int ff = rand();

	ofstream output;//outputt;//outputt;
	output.open("test" + to_string(L) + "-" + to_string(ff) + ".txt");
//	outputt.open("delta" + to_string(L) + "-" + to_string(ff) + ".txt");
//	ofstream niegh;
//	niegh.open("niegh" + to_string(L) + ".txt");


	
	for (int k = 0; k < iteration; k++)
	{
		vector<vector<double>> indexes;


	    vector<int> divisions;
		divisions = Divisors(indexes.size());
		int num_iter = divisions[ceil(divisions.size()/2)];
		int num_proc = int(N/num_iter);

		printf("num_proc = %d \n", num_proc);
		printf("num_iter = %d \n", num_iter);
//		printf("iter = %d \n", k);


		std::future<vector<double>> futuer;
		std::vector<vector<vector<double>>> futures;

		for (int j = 0; j< num_proc; j++)//num_proc
		{
			std::future<vector<vector<double>>> futu = std::async(std::launch::async, water_clusters, j, num_iter, indexes, total_lat, plat, PBS_Pos);
			auto results = futu.get();
			std::copy(results[0].begin(), results[0].end(),  back_inserter(final_total_landmass));
			std::copy(results[1].begin(), results[1].end(),  back_inserter(final_big_cluster));
			printf("%d \n", j);
		}
		
		//saving the data in npy format
		string path1 = "/home/complex/c++/Rewrite_c++/HKPBC/CMB_test_data/"; //"/share/users/m_movahed/cmb_cluster/s22_c++/"
		string path = "/home/complex/c++/Rewrite_c++/HKPBC/CMB_test_data/"; //"/share/users/m_movahed/cmb_cluster/s22_c++/"
		string filename = "total_landmass_" + to_string(one)+ "_" + to_string(PBS_counter) + ".txt";
		string filename1 = "big_cluster_" + to_string(one) + "_" + to_string(PBS_counter) + ".txt";
		std::ofstream fout(path + filename);
		fout.precision(17);
		std::copy(final_total_landmass.begin(), final_total_landmass.end(),std::ostream_iterator<double>(fout, "\n"));

		std::ofstream fout1(path1 + filename1);
		fout1.precision(17);
		std::copy(final_big_cluster.begin(), final_big_cluster.end(),std::ostream_iterator<double>(fout1, "\n"));

		return 0;
	
}


vector<vector<double>> levy(int steps, double C)
{
	vector<vector<double>> final;
	uniform_int_distribution<> distr(0, 1);
	normal_distribution<> d{0, (1/(2*c)^(0.5))}; //the first entery is mean and the second is the std
	random_device rd;
	mt19937 gen(rd()); //Mersenne Twister RNG

	for ( int i = 0; i < steps; i++)
	{
		vector<double> result;
		duoble angle = 2 * M_PI * distr(gen);
		result.push_back(angle);
		cout<<"angle: "<<angle<<endl;
		double r = (1/(d(gen))^2);
		while (r > L)
		{
			r = (1/(d(gen))^2);
			cout<<"r: "<<r<<endl;
		}
		result.push_back(r);
		double x = r * cos(angle);
		result.push_back(x);
		double y = r * sin(angle);
		result.push_back(y);
		final.push_back(result);
	}

	return final;

}

std::vector<std::vector<int>> split (const std::vector<int>& v) 
{
    int n = v.size();
    int size_max = n / Nsplit + (n % Nsplit != 0);
    std::vector<std::vector<int>> split;
    for (int ibegin = 0; ibegin < n; ibegin += size_max) {
        int iend = ibegin + size_max;
        if (iend > n) iend = n;
        split.emplace_back (std::vector<int>(v.begin() + ibegin, v.begin() + iend));
    }
    return split;
}


//divisors function to caculate the number of processes and iterations
vector<int> Divisors(int n)
{
    // Vector to store half of the divisors
    vector<int> v;
    for (int i = 1; i <= sqrt(n); i++) {
        if (n % i == 0) {
            // check if divisors are equal
            if (n / i != i) 
            {
                v.push_back(n / i);
            }
        }
    }
    return v;
}

