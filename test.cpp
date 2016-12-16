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
#include "gaussianhmm.hpp" 
#include "exponentialhmm.hpp" 
#include "gammahmm.hpp" 
#include "lognormalhmm.hpp" 
#include "paretohmm.hpp"
#include "laplacehmm.hpp" 
#include "discretehmm.hpp" 
#include "poissonhmm.hpp" 
#include "gumbelhmm.hpp" 
#include "weibullhmm.hpp" 
#include "extremevaluehmm.hpp" 
#include "multihmm.hpp" 

using namespace std;

template <typename T> void output_check(HMM<T> *hmm, vector<vector<double> > &A, vector<vector<double> > &B, vector<double> &pi){
	int N = hmm->N;
	int M = hmm->M;
	cout << hmm->type << " check" << endl;
	cout << "A fit :: A true" << endl;
	double max = 0;
	for(int i=0; i<N; ++i){
		for(int j=0; j<N; ++j){
			cout << hmm->A[i][j] << " ";
			double diff = fabs(hmm->A[i][j] - A[i][j]);
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
			cout << hmm->B[i][j] << " ";
			double diff = fabs(hmm->B[i][j] - B[i][j]);
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
			cout << hmm->pi[i] << "\t::\t" << pi[i] << "\n";
			double diff = fabs(hmm->pi[i] - pi[i]);
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
       gumbel\n";
	cerr << "./test -t all\n";
	cerr << "./test" << endl;
	exit(1);
}

int main(int argc,char **argv){
	
	int c;
	static struct option loptions[] =    {
        {"test",1,0,'t'},
        {0,0,0,0}
    };
	
	string test_string = "all";
	
    while ((c = getopt_long(argc, argv, "t:",loptions,NULL)) >= 0) {  
        switch (c)
        {
        case 't': test_string = (optarg); break;
        case '?': usage();
        default: 
			test_string = "all";
        }
    }

	transform(test_string.begin(), test_string.end(), test_string.begin(), ::tolower);    
    bool doall = (test_string == "all") ? true : false;
	
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
	
		gaussianHMM<double> hmm(N,0,1000000); 
		hmm.info();
		
		//Generate observation sequence with these paramters
		hmm.init();
		hmm.A = A; hmm.B = B; hmm.pi = pi;
		
		int T = 10000;
		vector<double> O(T);
		hmm.generate_seq(O, T, false);
		
		//Fit that observation sequence
		hmm.init();
		hmm.fit(O, 1e-5);
		hmm.sortparams();
		output_check<double>(&hmm, A, B, pi);
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
		
		exponentialHMM<double> hmm(N,0,1000000); 
		hmm.info();

		hmm.lb[0] = 5;
		hmm.lb[1] = 5;
		hmm.init();

		hmm.A = A; hmm.B = B; hmm.pi = pi;
		
		int T = 10000;
		vector<double> O(T);
		hmm.generate_seq(O, T, false);

		hmm.init();
		hmm.fit(O, 1e-10);
		hmm.sortparams();
		output_check<double>(&hmm, A, B, pi);
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
		

		gammaHMM<double> hmm(N,0,1000000); 
		hmm.info();
		hmm.init();

		hmm.A = A; hmm.B = B; hmm.pi = pi;
		
		int T = 10000;
		vector<double> O(T);
		hmm.generate_seq(O, T, false);
		
		hmm.init();
		hmm.fit(O, 1e-5);
		hmm.sortparams();
		output_check<double>(&hmm, A, B, pi);

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
		
		lognormalHMM<double> hmm(N,0,1000000); 
		hmm.info();
		hmm.init();

		hmm.A = A; hmm.B = B; hmm.pi = pi;
		
		int T = 10000;
		vector<double> O(T);
		hmm.generate_seq(O, T, false);
		
		hmm.init();
		hmm.fit(O, 1e-5);
		hmm.sortparams();
		output_check<double>(&hmm, A, B, pi);
		
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
		
		paretoHMM<double> hmm(N,0,1000000);
		hmm.info(); 
		hmm.init();

		hmm.A = A; hmm.B = B; hmm.pi = pi;
		hmm.lb[0] = 1;
		hmm.lb[1] = 1;
			
		int T = 10000;
		vector<double> O(T);
		hmm.generate_seq(O, T, false);
		
		hmm.init();
		hmm.fit(O, 1e-5);
		hmm.sortparams();
		output_check<double>(&hmm, A, B, pi);
		
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
		
		laplaceHMM<double> hmm(N,0,1000000); 
		hmm.info(); 
		hmm.init();

		hmm.A = A; hmm.B = B; hmm.pi = pi;
		
		srand(1234567);
		
		int T = 10000;
		vector<double> O(T);
		hmm.generate_seq(O, T, false);
		
		hmm.init();
		hmm.fit(O, 1e-10);
		hmm.sortparams();
		output_check<double>(&hmm, A, B, pi);
		
	}
	//Categorical
	if( doall || test_string == "discrete") {
		
		int N=2;
		int M=3;
		
		vector<vector<double> > A(N, vector<double>(N));
		vector<vector<double> > B(N, vector<double>(M));
		vector<double> pi(N);
		
		A[0][0] = 0.9; A[0][1] = 0.1; 
		A[1][0] = 0.1; A[1][1] = 0.9;

		B[0][0] = 0.8; B[0][1] = 0.0; B[0][2] = 0.2;
		B[1][0] = 0.1; B[1][1] = 0.4; B[1][2] = 0.5;
		
		pi[0] = 0.6; pi[1] = 0.4;
		
		multinomialHMM<int> hmm(N,M,0,1000000); 
		hmm.info(); 
		hmm.init();

		hmm.A = A; hmm.B = B; hmm.pi = pi;
		
		int T = 20000;
		vector<int> O(T);
		hmm.generate_seq(O, T, false);
		
		hmm.init();
		hmm.fit(O, 1e-5);
		hmm.sortparams();
		output_check(&hmm, A, B, pi);
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
		

		poissonHMM<int> hmm(N,0,1000000); 
		hmm.info(); 
		hmm.init();

		hmm.A = A; hmm.B = B; hmm.pi = pi;
		
		int T = 10000;
		vector<int> O(T);
		hmm.generate_seq(O, T, false);
		
		hmm.init();
		hmm.fit(O, 1e-10);
		hmm.sortparams();
		output_check(&hmm, A, B, pi);
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
		

		gumbelHMM<double> hmm(N,0,1000000); 
		hmm.info(); 
		hmm.init();

		hmm.A = A; hmm.B = B; hmm.pi = pi;
		
		srand(123456789);
		
		int T = 10000;
		vector<double> O(T);
		hmm.generate_seq(O, T, false);
		
		hmm.init();
		hmm.fit(O, 1e-10);
		hmm.sortparams();
		output_check(&hmm, A, B, pi);
			
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
		
		weibullHMM<double> hmm(N,0,1000000); 
		hmm.info(); 
		hmm.init();

		hmm.A = A; hmm.B = B; hmm.pi = pi;
				
		srand(123456789);
		
		int T = 10000;
		vector<double> O(T);
		hmm.generate_seq(O, T, false);
		
		hmm.init();
		hmm.fit(O, 1e-5);
		hmm.sortparams();
		output_check(&hmm, A, B, pi);
			
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


		extremevalueHMM<double> hmm(N,0,1000000); 		
		hmm.info(); 
		hmm.init();

		hmm.A = A; hmm.B = B; hmm.pi = pi;
					
		srand(123456789);
		
		int T = 10000;
		vector<double> O(T);
		hmm.generate_seq(O, T, false);
          
		hmm.init();
		hmm.print_iter = true;
		hmm.fit(O, 1e-5);
		hmm.sortparams();
		output_check(&hmm, A, B, pi);
			
	}
	//multi emission
	if( doall || test_string == "multi") {
		
		int N=3;
		int M=2;
		
		vector<vector<double> > A(N, vector<double>(N));
		vector<vector<double> > B(N, vector<double>(M));
		vector<double> pi(N);
		
		//A[0][0] = 0.9; A[0][1] = 0.1; 
		//A[1][0] = 0.1; A[1][1] = 0.9;
		
		A[0][0] = 0.9;  A[0][1] = 0.05;  A[0][2] = 0.05; 
		A[1][0] = 0.05; A[1][1] = 0.9;   A[1][2] = 0.05; 
		A[2][0] = 0.05; A[2][1] = 0.05;  A[2][2] = 0.9; 

		B[0][0] = 100;  B[0][1] = 0;
		B[1][0] = 5; B[1][1] = 1;  
		B[2][0] = 3; B[2][1] = 0.1;  

		//pi[0] = 0.6; pi[1] = 0.4;
		pi[0] = 0.6; pi[1] = 0.2; pi[2] = 0.2;

		
		vector<vector<double> > Bexp(1, vector<double>(1));
		Bexp[0][0] = 100;  
		vector<vector<double> > Blognorm(1, vector<double>(2));
		Blognorm[0][0] = 5; Blognorm[0][1] = 1; 
		vector<vector<double> > Bnorm(1, vector<double>(2));
		Bnorm[0][0] = 3; Bnorm[0][1] = 0.1; 

			
	
		exponentialHMM<double> hmm_exp(1,0,1000000); hmm_exp.B = Bexp;
		gaussianHMM<double> hmm_norm(1,0,1000000); hmm_norm.B = Bnorm;
		lognormalHMM<double> hmm_lognorm(1,0,1000000); hmm_lognorm.B = Blognorm;
		//vector< HMM<double>* > hmms = { &hmm_exp };
		//vector< HMM<double>* > hmms = { &hmm_exp, &hmm_norm };
		vector< HMM<double>* > hmms = { &hmm_exp, &hmm_lognorm, &hmm_norm };
		
		multiHMM<double> hmm(hmms,0,1000000); 
		hmm.init();
		hmm.info();
		hmm.A = A; hmm.pi = pi; hmm.setB(B);


		int T = 30000;
		vector<double> O(T);
		hmm.generate_seq(O, T, false);
		
		hmm.init();
		hmm.print_iter = true;
		hmm.fit(O, 1e-5);
		hmm.sortparams();
		output_check(&hmm, A, B, pi);
	}
	
	return 0;
	
	
} 


