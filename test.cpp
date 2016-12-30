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
#include <math.h>
#include <getopt.h>
#include "cdhmm.hpp" 


using namespace std;
using namespace cdHMM;

template <typename T> void output_check(HMM<T> *hmm, vector<vector<double> > &A, vector<vector<double> > &B, vector<double> &pi){
	int N = hmm->N;
	int M = hmm->M;
	cout << hmm->type << " check" << fixed << endl;
	cout << "A fit :: A true" << endl;
	double max = 0;
	for(int i=0; i<N; ++i){
		for(int j=0; j<N; ++j){
			cout << hmm->maxA[i][j] << " ";
			double diff = fabs(hmm->maxA[i][j] - A[i][j]);
			if( diff > max ){ max = diff; }
		} 
		cout << "\t::\t";
		for(int j=0; j<N; ++j){
			cout << A[i][j] << " ";
		} 
		cout << "\n";
	}
	cout << "Max A diff = " << max << endl;
	cout << "B fit :: B true" << endl;
	max = 0;
	for(int i=0; i<N; ++i){
		for(int j=0; j<M; ++j){
			cout << hmm->maxB[i][j] << " ";
			double diff = fabs(hmm->maxB[i][j] - B[i][j]);
			if( diff > max ){ max = diff; }
		} 
		cout << "\t::\t";
		for(int j=0; j<M; ++j){
			cout << B[i][j] << " ";
		} 
		cout << "\n";
	}
	cout << "Max B diff = " << max << endl;
	cout << "pi fit :: pi true" << endl;
	max = 0;
	for(int i=0; i<N; ++i){
			cout << hmm->maxpi[i] << "\t::\t" << pi[i] << "\n";
			double diff = fabs(hmm->maxpi[i] - pi[i]);
			if( diff > max ){ max = diff; }
	} 
	cout << "Max pi diff = " << max << endl;
	
}


void usage(){
	cerr << "usage:\n";
	cerr << "./test -t type\n";
	cerr << 
"type = gaussian, exponential, gamma,\n\
       lognormal, pareto, laplace,\n\
       discrete, poisson,\n\
       gumbel, weibull, extremevalue,\n\
       multi.\n";
	cerr << "./test -t all\n";
	exit(1);
}

template <typename T> void run_hmm(HMM<T> *hmm, vector<vector<double> > &A, 
vector<vector<double> > &B, vector<double> &pi, 
int num_obs, double eps, int num_restarts){
	
	//Print basic info
	hmm->info();

	//Generate observation sequence with these paramters
	hmm->init();

	hmm->setA( A );
	hmm->setpi( pi ); 
	hmm->setB( B ); 

	vector<T> O(num_obs);
	hmm->generate_seq(O, num_obs, false);

	
	//Fit that observation sequence
	cerr << "Running";
	for(int i=0; i<num_restarts; ++i){
		cerr << ".";
		hmm->init();
		//hmm->print_iter = true;
		hmm->fit(O, eps);
		hmm->sortparams();
	}
	cerr << endl;
	
	output_check<T>(hmm, A, B, pi);

}

int main(int argc,char **argv){
	
	int c;
	static struct option loptions[] =    {
        {"test",1,0,'t'},
        {0,0,0,0}
    };
	
	string test_string = "none";
	
    while ((c = getopt_long(argc, argv, "t:",loptions,NULL)) >= 0) {  
        switch (c)
        {
        case 't': test_string = (optarg); break;
        case '?': usage();
        default: 
			usage();
        }
    }
    if( test_string == "none" ){ usage(); }

	transform(test_string.begin(), test_string.end(), test_string.begin(), ::tolower);    
    bool doall = (test_string == "all") ? true : false;
	
	int hmm_max_iters = 1000;
	//Gaussian emission
	if( doall || test_string == "gaussian") {
		
		int N=2;
		int M=2;
		
		vector<vector<double> > A(N, vector<double>(N));
		vector<vector<double> > B(N, vector<double>(M));
		vector<double> pi(N);
		
		A[0][0] = 0.9; A[0][1] = 0.1; 
		A[1][0] = 0.1; A[1][1] = 0.9;

		B[0][0] = 0.9; B[0][1] = 0.1;  
		B[1][0] = 0.1; B[1][1] = 0.9; 
		
		pi[0] = 0.6; pi[1] = 0.4;
	
		gaussianHMM<double> hmm(N,0,hmm_max_iters); 
		run_hmm(&hmm,  A, B, pi, 10000, 1e-5, 10);
		
	}
	//Exponential emission
	if( doall || test_string == "exponential") {
		
		int N=2;
		int M=1;
		
		vector<vector<double> > A(N, vector<double>(N));
		vector<vector<double> > B(N, vector<double>(M));
		vector<double> pi(N);
		
		A[0][0] = 0.9; A[0][1] = 0.1; 
		A[1][0] = 0.1; A[1][1] = 0.9;

		B[0][0] = 0.9;
		B[1][0] = 0.1;  
		
		pi[0] = 0.6; pi[1] = 0.4;
		
		exponentialHMM<double> hmm(N,0,hmm_max_iters);
		hmm.lb[0] = 5; hmm.lb[1] = 5; 
		run_hmm(&hmm,  A, B, pi, 10000, 1e-5, 10);

	}
	//Gamma emission
	if( doall || test_string == "gamma") {
		
		int N=2;
		int M=2;
		
		vector<vector<double> > A(N, vector<double>(N));
		vector<vector<double> > B(N, vector<double>(M));
		vector<double> pi(N);
		
		A[0][0] = 0.9; A[0][1] = 0.1; 
		A[1][0] = 0.1; A[1][1] = 0.9;
		
		B[0][0] = 4.5; B[0][1] = 1.5;  
		B[1][0] = 0.5; B[1][1] = 1.0; 
		
		pi[0] = 0.6; pi[1] = 0.4;
		

		gammaHMM<double> hmm(N,0,hmm_max_iters); 
		run_hmm(&hmm,  A, B, pi, 10000, 1e-5, 10);

	}
	//Log Normal
	if( doall || test_string == "lognormal") {
		
		int N=2;
		int M=2;
		
		vector<vector<double> > A(N, vector<double>(N));
		vector<vector<double> > B(N, vector<double>(M));
		vector<double> pi(N);
		
		A[0][0] = 0.9; A[0][1] = 0.1; 
		A[1][0] = 0.1; A[1][1] = 0.9;

		B[0][0] = 10; B[0][1] = 0.1;
		B[1][0] = 0; B[1][1] = 0.1; 
		
		pi[0] = 0.6; pi[1] = 0.4;
		
		lognormalHMM<double> hmm(N,0,hmm_max_iters); 
		run_hmm(&hmm,  A, B, pi, 10000, 1e-5, 10);

		
	}//Pareto
	if( doall || test_string == "pareto") {
		
		int N=2;
		int M=1;
		
		vector<vector<double> > A(N, vector<double>(N));
		vector<vector<double> > B(N, vector<double>(M));
		vector<double> pi(N);
		
		A[0][0] = 0.9; A[0][1] = 0.1; 
		A[1][0] = 0.1; A[1][1] = 0.9;

		B[0][0] = 2; 
		B[1][0] = 1; 
		
		pi[0] = 0.6; pi[1] = 0.4;
		
		paretoHMM<double> hmm(N,0,hmm_max_iters);
		hmm.lb[0] = 1; hmm.lb[1] = 1;
		cerr << "lower bounds = " << hmm.lb[0] << " " << hmm.lb[1] << endl; 
		
		run_hmm(&hmm,  A, B, pi, 10000, 1e-5, 10);

		
	}
	//Laplace
	if( doall || test_string == "laplace") {
		
		int N=2;
		int M=2;
		
		vector<vector<double> > A(N, vector<double>(N));
		vector<vector<double> > B(N, vector<double>(M));
		vector<double> pi(N);
		
		A[0][0] = 0.9; A[0][1] = 0.1; 
		A[1][0] = 0.1; A[1][1] = 0.9;

		B[0][0] = 15; B[0][1] = 0.5;
		B[1][0] = 0; B[1][1] = 1; 
		
		pi[0] = 0.6; pi[1] = 0.4;
		
		laplaceHMM<double> hmm(N,0,hmm_max_iters); 
		run_hmm(&hmm,  A, B, pi, 10000, 1e-5, 10);

		
	}
	//Categorical
	if( doall || test_string == "discrete") {
		
		int N=2;
		int M=2;
		
		vector<vector<double> > A(N, vector<double>(N));
		vector<vector<double> > B(N, vector<double>(M));
		vector<double> pi(N);
		
		A[0][0] = 0.9; A[0][1] = 0.1; 
		A[1][0] = 0.1; A[1][1] = 0.9;

		B[0][0] = 0.8; B[0][1] = 0.2; 
		B[1][0] = 0.2; B[1][1] = 0.8; 
		
		pi[0] = 0.6; pi[1] = 0.4;
		
		discreteHMM<int> hmm(N,M,0,hmm_max_iters); 
		run_hmm(&hmm,  A, B, pi, 10000, 1e-5, 10);

	}
	//Poisson emission
	if( doall || test_string == "poisson") {
		
		int N=2;
		int M=1;
		
		vector<vector<double> > A(N, vector<double>(N));
		vector<vector<double> > B(N, vector<double>(M));
		vector<double> pi(N);
		
		A[0][0] = 0.9; A[0][1] = 0.1; 
		A[1][0] = 0.1; A[1][1] = 0.9;

		B[0][0] = 10;
		B[1][0] = 1;  
		
		pi[0] = 0.6; pi[1] = 0.4;
		

		poissonHMM<int> hmm(N,0,hmm_max_iters); 
		run_hmm(&hmm,  A, B, pi, 10000, 1e-5, 10);

	} 
	//Gumbel
	if( doall || test_string == "gumbel") {
		
		int N=2;
		int M=2;
		
		vector<vector<double> > A(N, vector<double>(N));
		vector<vector<double> > B(N, vector<double>(M));
		vector<double> pi(N);
		
		A[0][0] = 0.9; A[0][1] = 0.1; 
		A[1][0] = 0.1; A[1][1] = 0.9;

		B[0][0] = 40; B[0][1] = 4;
		B[1][0] = 0.5; B[1][1] = 2; 
		
		pi[0] = 0.6; pi[1] = 0.4;
		

		gumbelHMM<double> hmm(N,0,hmm_max_iters); 
		run_hmm(&hmm,  A, B, pi, 10000, 1e-5, 10);

			
	}
	//Weibull
	if( doall || test_string == "weibull") {
		
		int N=2;
		int M=2;
		
		vector<vector<double> > A(N, vector<double>(N));
		vector<vector<double> > B(N, vector<double>(M));
		vector<double> pi(N);

		A[0][0] = 0.9; A[0][1] = 0.1; 
		A[1][0] = 0.1; A[1][1] = 0.9;

		B[0][0] = 5; B[0][1] = 1;
		B[1][0] = 0.9; B[1][1] = 2; 

		pi[0] = 0.6; pi[1] = 0.4;
		
		weibullHMM<double> hmm(N,0,hmm_max_iters); 
		run_hmm(&hmm,  A, B, pi, 10000, 1e-2, 10);

			
	}
	//Extreme value
	if( doall || test_string == "extremevalue") {
		
		int N=2;
		int M=3;
		
		vector<vector<double> > A(N, vector<double>(N));
		vector<vector<double> > B(N, vector<double>(M));
		vector<double> pi(N);

		//A[0][0] = 1;
		A[0][0] = 0.9; A[0][1] = 0.1; 
		A[1][0] = 0.1; A[1][1] = 0.9;


		B[0][0] = 500; B[0][1] = 10; B[0][2] = 0.125;
		B[1][0] = 100; B[1][1] = 10; B[1][2] = 0.125;

		//pi[0] = 1;
		pi[0] = 0.6; pi[1] = 0.4;

		extremevalueHMM<double> hmm(N,0,hmm_max_iters); 		
		run_hmm(&hmm,  A, B, pi, 10000, 1e-2, 10);

			
	}
	//multi emission
	if( doall || test_string == "multi") {
		
		int N=3;
		int M=2;
		
		vector<vector<double> > A(N, vector<double>(N));
		vector<vector<double> > B(N, vector<double>(M));
		vector<double> pi(N);
		
		A[0][0] = 0.9;  A[0][1] = 0.05;  A[0][2] = 0.05; 
		A[1][0] = 0.05; A[1][1] = 0.9;   A[1][2] = 0.05; 
		A[2][0] = 0.05; A[2][1] = 0.05;  A[2][2] = 0.9; 

		B[0][0] = 100;  B[0][1] = 0;
		B[1][0] = 5; B[1][1] = 1;  
		B[2][0] = 3; B[2][1] = 0.1;  


		pi[0] = 0.6; pi[1] = 0.2; pi[2] = 0.2;

	
		exponentialHMM<double> hmm_exp(1,0,hmm_max_iters); 
		lognormalHMM<double> hmm_lognorm(1,0,hmm_max_iters); 
		gaussianHMM<double> hmm_norm(1,0,hmm_max_iters);
		vector< HMM<double>* > hmms = { &hmm_exp, &hmm_lognorm, &hmm_norm };
		
		multiHMM<double> hmm(hmms,0,hmm_max_iters);
		run_hmm(&hmm,  A, B, pi, 30000, 1e-5, 10);


	}
	
	return 0;
	
	
} 


