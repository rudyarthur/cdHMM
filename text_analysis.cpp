#include <iostream>
#include <vector>
#include <algorithm>
#include <map>
#include <fstream>
#include <sstream>
#include <iterator>
#include <string>
#include <locale>         
#include <limits>
#include <ctype.h>
#include <math.h>
#include <getopt.h>
#include "discretehmm.hpp" 

using namespace std;




void usage(){
	cerr << "usage:" << endl;
	cerr << "./text_analysis filename.txt" << endl; 
	exit(1);
}

void read_txt(string filename, vector< char > &data){
	
	data.resize(0);
	ifstream in(filename.c_str());
	
	if(!in.is_open()){
		cout << "Failed to open file." << endl;
		exit(1);
	}
	
	char c;
	while(in.get(c)) // loop through the file character by character
	{	
		if( isalpha(c) || isspace(c) ){
			data.push_back(c);
		}
	}	
		
	in.close();
}

int main(int argc,char **argv){

	if(argc < 2){ usage(); }
	string filename = argv[1];
	
	vector<char> data;
	read_txt(filename, data);
	
	for(int i=0; i<data.size(); ++i){
		cout << i << " " << data[i] << endl;
	}
	/*int N=2;
	int M=2;
		
	multinomialHMM<int> hmm(N,M,0,hmm_max_iters); 
		
	//Print basic info
	hmm->info();

	//Generate observation sequence with these paramters
	hmm->init();

	hmm->A = A;
	hmm->pi = pi; 
	hmm->setB( B ); 

	vector<T> O(num_obs);
	hmm->generate_seq(O, num_obs, false);

	//Fit that observation sequence
	cerr << "Running";
	for(int i=0; i<num_restarts; ++i){
		cerr << ".";
		hmm->init();
		hmm->fit(O, eps);
		hmm->sortparams();
	}
	cerr << endl;
	
	output_check<T>(hmm, A, B, pi);*/
	
	return 0;
	
	
} 


