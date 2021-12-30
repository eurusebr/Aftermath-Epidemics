#include"HKPBC.h"
#include <iostream>
#include <cstdlib>
#include <cmath>
#include <assert.h>
#include <vector>
#include <random>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <ctime>
#include <stdio.h>
#include <stdlib.h>
#include <cstdlib>
#include <iostream>
#include <thread>
#include <map>
#include <string>
#include <valarray>
#include <future>
#include <iterator>


#define Lx 20
#define Ly 20
#define N Lx*Ly //total nodes
#define steps 10*N //steps of the LRW
#define iteration 50 //total iterations
#define M_PI 3.14159265358979323846
#define Nsplit 10
#define Levy_Alpha 0.1


using namespace std;


vector<vector<double>> levy(int size, double C);
vector<vector<int>> split (const std::vector<int>& v); 
vector<int> Divisors(int n);
vector<int> path_maker(vector<vector<double>>& total);
vector<vector<double>> find_clusters(int ii, int const& num_iter, vector<int> indexes, int PBS_Pos);


int main(void)
{
	srand(time(0));
	random_device rd;
	mt19937 gen(rd()); //Mersenne Twister RNG
	int ff = rand();

	ofstream output;//outputt;//outputt;
	output.open("test" + to_string(Lx) + "-" + to_string(ff) + ".txt");
//	outputt.open("delta" + to_string(L) + "-" + to_string(ff) + ".txt");
//	ofstream niegh;
//	niegh.open("niegh" + to_string(L) + ".txt");


	
	for (int k = 0; k < iteration; k++)
	{
		vector<vector<double>> total;
		vector<int> indexes, divisions;
		vector<float> final_total_landmass, final_big_cluster;

		double C = Levy_Alpha;
		total = levy(steps, C);
		indexes = path_maker(total);
		divisions = Divisors(indexes.size());
		int num_iter = divisions[ceil(divisions.size()/2)];
		int num_proc = int(steps/num_iter);

		printf("num_proc = %d \n", num_proc);
		printf("num_iter = %d \n", num_iter);

		auto PBS_index = split(indexes);
    	int PBS_counter = two;
    	int PBS_Pos = PBS_counter*PBS_index[0].size();


		// future<vector<double>> futuer;
		// vector<vector<vector<double>>> futures;

		for (int j = 0; j< num_proc; j++)//num_proc
		{
			std::future<vector<vector<double>>> futu = std::async(std::launch::async, find_clusters, j, num_iter, indexes, PBS_Pos);
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

	}
	return 0;
}


vector<vector<double>> levy(int size, double C)
{
	vector<vector<double>> final;
	uniform_int_distribution<> distr(0, 1);
	// normal_distribution<> d{0, (1/pow(2*C, 0.5))}; //the first entery is mean and the second is the std
	random_device rd, rand;
	mt19937 gen(rd()); //Mersenne Twister RNG
	mt19937 genn(rand()); //Mersenne Twister RNG


	for ( int i = 0; i < steps; i++)
	{
		vector<double> result;
		double angle = 2 * M_PI * distr(gen);
		result.push_back(angle);
		cout<<"angle: "<<angle<<endl;
		double r = (1/pow(distr(genn), -1/C));

		while (r > Lx | r == 0)
		{
			r = (1/pow(distr(genn), -1/C));
		}
		cout<<"r: "<<r<<endl;
		result.push_back(r);
		double x = r * cos(angle);
		result.push_back(x);
		double y = r * sin(angle);
		result.push_back(y);
		final.push_back(result);
	}

	return final;

}

vector<vector<int>> split (const std::vector<int>& v) 
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

vector<int> path_maker(vector<vector<double>>& total )
{
	vector<int> matrix;
	random_device rd;
	mt19937 gen(rd()); //Mersenne Twister RNG
	uniform_int_distribution<int> ran_pos(0, N-1);
	int x = ran_pos(gen);
	int y = ran_pos(gen);
	matrix.push_back((Ly*x + y));
	for (int i = 0; i < total.size(); i++)
	{
		x = x + total[i][2] == Lx ? 0 : x + total[i][2];
		y = y + total[i][3] == Ly ? 0 : y + total[i][3];
		matrix.push_back((Ly*x + y));	
	}
	return matrix;
}

vector<vector<double>> find_clusters(int ii, int const& num_iter, vector<int> indexes, int PBS_Pos)
{
    int l = ii*num_iter;
    vector<double> total_landmass; //fraction of land
    vector<double> big_cluster; //bigger cluster
    vector<vector<double>> result; //to return both of the above
    double pp; // to calculate the fraction of land
    int s; //vector of the biggest labels in th descending order
    HKPBC object(Lx,Ly,N); //introducing the hk class


    for (int ii=l; ii< l + num_iter; ii++)
    {
        //making the array of labeling
        vector<int> myBoolArray(N);
        for (int m=0; m< (ii+1)+PBS_Pos; m++)
        {
            myBoolArray[indexes[m]] = 1;
        }
        //labeling the array and returning the cluster labels
        s = object.HK(myBoolArray);

        pp = ((ii+1)+PBS_Pos) / double(N); //normal fraction of land
        big_cluster.push_back(1.0*s); // biggest cluster
        total_landmass.push_back(pp);
    }

    //returning the result
    result.push_back(total_landmass);
    result.push_back(big_cluster);

    return result;
}
