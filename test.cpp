#include <iostream>
#include <vector>
#include <map>
#include <fstream>
#include <sstream>
#include <iterator>
#include <string>
#include <limits>
#include <math.h>
#include "gaussianhmm.hpp" 
#include "exponentialhmm.hpp" 
#include "gammahmm.hpp" 
/*#include "poissonhmm.hpp" 
#include "lognormalhmm.hpp" 
#include "discretehmm.hpp" 
#include "laplacehmm.hpp" 
#include "jointhmm.hpp" 
#include "gumbelhmm.hpp" 
#include "weibullhmm.hpp" 
#include "extremevaluehmm.hpp" 
#include "paretohmm.hpp" 
#include "multihmm.hpp" */

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

int main(int argc,char **argv){
	
	bool test_gaussian = true;
	bool test_exponential = false;
	bool test_gamma = true;
	bool test_poisson = false;
	bool test_multinomial = false;
	bool test_lognormal = false;
	bool test_joint = false;
	bool test_laplace = false;
	bool test_gumbel = false;
	bool test_weibull = false;
	bool test_extremevalue = true;
	bool test_pareto = false;
	bool test_multi = false;
	
	//Gaussian emission
	if( test_gaussian ) {
		
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
	if( test_exponential ) {
		
		int N=2;
		int M=1;
		
		vector<vector<double> > A(N, vector<double>(N));
		vector<vector<double> > B(N, vector<double>(M));
		vector<double> pi(N);
		
		A[0][0] = 0.9; A[0][1] = 0.1; 
		A[1][0] = 0.1; A[1][1] = 0.9;

		B[0][0] = 0.1;
		B[1][0] = 0.9;  
		
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
	if( test_gamma ) {
		
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

	}//Poisson emission
	/*if( test_poisson ) {
		
		int N=2;
		int M=1;
		
		vector<vector<double> > A(N, vector<double>(N));
		vector<vector<double> > B(N, vector<double>(M));
		vector<double> pi(N);
		
		A[0][0] = 0.9; A[0][1] = 0.1; 
		A[1][0] = 0.1; A[1][1] = 0.9;

		B[0][0] = 1;
		B[1][0] = 10;  
		
		pi[0] = 0.6; pi[1] = 0.4;
		
		cout << "Poisson Emission" << endl;

		poissonHMM<int> hmm(N,0,1000000); 
		hmm.init();

		hmm.A = A; hmm.B = B; hmm.pi = pi;
		
		int T = 10000;
		vector<int> O(T);
		hmm.generate_seq(O, T, false);
		
		hmm.init();
		hmm.fit(O, false, 1e-10);
		cout << "A = "; hmm.printA();
		cout << "B = "; hmm.printB();
		cout << "pi = "; hmm.printpi();
		
		hmm.init();
		hmm.fit(O, true, 1e-10);
		cout << "A = "; hmm.printA();
		cout << "B = "; hmm.printB();
		cout << "pi = "; hmm.printpi();
	} 
	//Categorical
	if( test_multinomial ) {
		
		int N=2;
		int M=3;
		
		vector<vector<double> > A(N, vector<double>(N));
		vector<vector<double> > B(N, vector<double>(M));
		vector<double> pi(N);
		
		A[0][0] = 0.9; A[0][1] = 0.1; 
		A[1][0] = 0.1; A[1][1] = 0.9;

		B[0][0] = 0.1; B[0][1] = 0.4; B[0][2] = 0.5;
		B[1][0] = 0.8; B[1][1] = 0.0; B[1][2] = 0.2;
		
		pi[0] = 0.6; pi[1] = 0.4;
		
		cout << "Multinomial Emission" << endl;

		multinomialHMM<int> hmm(N,M,0,1000000); 
		hmm.init();

		hmm.A = A; hmm.B = B; hmm.pi = pi;
		
		int T = 10000;
		vector<int> O(T);
		hmm.generate_seq(O, T, false);
		
		hmm.init();
		hmm.fit(O, false, 1e-10);
		cout << "A = "; hmm.printA();
		cout << "B = "; hmm.printB();
		cout << "pi = "; hmm.printpi();
		
		hmm.init();
		hmm.fit(O, true, 1e-10);
		cout << "A = "; hmm.printA();
		cout << "B = "; hmm.printB();
		cout << "pi = "; hmm.printpi();
	}
	//Log Normal
	if( test_lognormal ) {
		
		int N=2;
		int M=2;
		
		vector<vector<double> > A(N, vector<double>(N));
		vector<vector<double> > B(N, vector<double>(M));
		vector<double> pi(N);
		
		A[0][0] = 0.9; A[0][1] = 0.1; 
		A[1][0] = 0.1; A[1][1] = 0.9;

		B[0][0] = 0; B[0][1] = 0.1;
		B[1][0] = 10; B[1][1] = 0.1; 
		
		pi[0] = 0.6; pi[1] = 0.4;
		
		cout << "Lognormal Emission" << endl;

		lognormalHMM<double> hmm(N,0,1000000); 
		hmm.init();

		hmm.A = A; hmm.B = B; hmm.pi = pi;
		
		int T = 10000;
		vector<double> O(T);
		hmm.generate_seq(O, T, false);
		
		hmm.init();
		hmm.fit(O, false, 1e-10);
		cout << "A = "; hmm.printA();
		cout << "B = "; hmm.printB();
		cout << "pi = "; hmm.printpi();
		
		hmm.init();
		hmm.fit(O, true, 1e-10);
		cout << "A = "; hmm.printA();
		cout << "B = "; hmm.printB();
		cout << "pi = "; hmm.printpi();
	}
	//joint emission
	if( test_joint ) {
		
		int N=2;
		int M=2;
		
		vector<vector<double> > A(N, vector<double>(N));
		vector<double> pi(N);
		
		A[0][0] = 0.9; A[0][1] = 0.1; 
		A[1][0] = 0.1; A[1][1] = 0.9;

		vector<vector<double> > Bexp(1, vector<double>(1));
		Bexp[0][0] = 100;  
		vector<vector<double> > Bnorm(1, vector<double>(2));
		Bnorm[0][0] = 3; Bnorm[0][1] = 0.1; 
		
		pi[0] = 0.6; pi[1] = 0.4;
	
		exponentialHMM<double> hmm_exp(1,0,1000000); 
		gaussianHMM<double> hmm_norm(1,0,1000000); 
		jointHMM<double> hmm(&hmm_exp,&hmm_norm,0,1000000); 
		hmm.init();
		hmm.A = A; hmm.pi = pi;
		hmm_exp.B = Bexp; hmm_norm.B = Bnorm;

		cout << "Joint Emission" << endl;

		int T = 10000;
		vector<double> O(T);
		hmm.generate_seq(O, T, false);
		
		
		hmm.init();
		hmm.fit(O, false);
		cout << "A = "; hmm.printA();
		cout << "B = "; hmm.printB();
		cout << "pi = "; hmm.printpi();
		
		hmm.init();
		hmm.fit(O, true);
		cout << "A = "; hmm.printA();
		cout << "B = "; hmm.printB();
		cout << "pi = "; hmm.printpi();
	}
	//Laplace
	if( test_laplace ) {
		
		int N=2;
		int M=2;
		
		vector<vector<double> > A(N, vector<double>(N));
		vector<vector<double> > B(N, vector<double>(M));
		vector<double> pi(N);
		
		A[0][0] = 0.9; A[0][1] = 0.1; 
		A[1][0] = 0.1; A[1][1] = 0.9;

		B[0][0] = 0; B[0][1] = 1;
		B[1][0] = 15; B[1][1] = 0.5; 
		
		pi[0] = 0.6; pi[1] = 0.4;
		
		cout << "Laplace Emission" << endl;

		laplaceHMM<double> hmm(N,0,1000000); 
		hmm.init();

		hmm.A = A; hmm.B = B; hmm.pi = pi;
		
		srand(1234567);
		
		int T = 10000;
		vector<double> O(T);
		hmm.generate_seq(O, T, false);
		
		hmm.init(); 
		hmm.fit(O, false, 1e-10);
		cout << "A = "; hmm.printA();
		cout << "B = "; hmm.printB();
		cout << "pi = "; hmm.printpi();
		
		hmm.init(); 
		hmm.fit(O, true, 1e-10);
		cout << "A = "; hmm.printA();
		cout << "B = "; hmm.printB();
		cout << "pi = "; hmm.printpi();
		
	}
	//Gumbel
	if( test_gumbel ) {
		
		int N=2;
		int M=2;
		
		vector<vector<double> > A(N, vector<double>(N));
		vector<vector<double> > B(N, vector<double>(M));
		vector<double> pi(N);
		
		A[0][0] = 0.9; A[0][1] = 0.1; 
		A[1][0] = 0.1; A[1][1] = 0.9;

		B[0][0] = 0.5; B[0][1] = 2;
		B[1][0] = 40; B[1][1] = 4; 
		
		pi[0] = 0.6; pi[1] = 0.4;
		
		cout << "Gumbel Emission" << endl;

		gumbelHMM<double> hmm(N,0,1000000); 
		hmm.init();

		hmm.A = A; hmm.B = B; hmm.pi = pi;
		
		srand(123456789);
		
		int T = 10000;
		vector<double> O(T);
		hmm.generate_seq(O, T, false);
		
		hmm.init(); 
		hmm.fit(O, false, 1e-5);
		cout << "A = "; hmm.printA();
		cout << "B = "; hmm.printB();
		cout << "pi = "; hmm.printpi();

		hmm.init(); 
		hmm.fit(O, true, 1e-5);
		cout << "A = "; hmm.printA();
		cout << "B = "; hmm.printB();
		cout << "pi = "; hmm.printpi();
			
	}
	//Weibull
	if( test_weibull ) {
		
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
		
		cout << "Weibull Emission" << endl;

		weibullHMM<double> hmm(N,0,1000000); 
		hmm.init();

		hmm.A = A; hmm.B = B; hmm.pi = pi;
				
		srand(123456789);
		
		int T = 10000;
		vector<double> O(T);
		hmm.generate_seq(O, T, false);
		
		hmm.init(); 
		
		hmm.fit(O, false, 1e-5);
		cout << "A = "; hmm.printA();
		cout << "B = "; hmm.printB();
		cout << "pi = "; hmm.printpi();

		hmm.init(); 
		hmm.fit(O, true, 1e-5);
		cout << "A = "; hmm.printA();
		cout << "B = "; hmm.printB();
		cout << "pi = "; hmm.printpi();
			
	}
	//Extreme value
	if( test_extremevalue ) {
		
		int N=2;
		int M=3;
		
		vector<vector<double> > A(N, vector<double>(N));
		vector<vector<double> > B(N, vector<double>(M));
		vector<double> pi(N);

		//A[0][0] = 1;
		A[0][0] = 0.9; A[0][1] = 0.1; 
		A[1][0] = 0.1; A[1][1] = 0.9;

		B[0][0] = 10000; B[0][1] = 3000; B[0][2] = 0.125;
		B[1][0] = -10000; B[1][1] = 3000; B[1][2] = 0.125 ;

		//pi[0] = 1;
		pi[0] = 0.6; pi[1] = 0.4;
		
		//test GEV functions
		/*  gsl_vector *theta;
		  theta = gsl_vector_alloc (3);
		  gsl_vector_set (theta, 0, B[0][0]);
		  gsl_vector_set (theta, 1, B[0][1]);
		  gsl_vector_set (theta, 2, B[0][2]);
		  
		  gsl_vector *dl = gsl_vector_alloc (3);
		  gsl_matrix *J = gsl_matrix_alloc (3,3);
   
		  vector<double> O = {9000, 9500, 10000};
          double nlog = GEV_Lhood (theta, O, dl, J);
     
          cout << nlog << "\n" << endl;
          cout << gsl_vector_get (dl, 0) << " " << gsl_vector_get (dl, 1) << " " << gsl_vector_get (dl, 2) << "\n" << endl;
          cout << gsl_matrix_get (J, 0, 0) << " " << gsl_matrix_get (J, 0, 1) << " " << gsl_matrix_get (J, 0, 2) << endl;
          cout << gsl_matrix_get (J, 1, 0) << " " << gsl_matrix_get (J, 1, 1) << " " << gsl_matrix_get (J, 1, 2) << endl;
          cout << gsl_matrix_get (J, 2, 0) << " " << gsl_matrix_get (J, 2, 1) << " " << gsl_matrix_get (J, 2, 2) << endl;
        */
        
      
		//cout << "Extreme Value Emission" << endl;
/*
		extremevalueHMM<double> hmm(N,0,1000000); 
		hmm.init();

		hmm.A = A; hmm.B = B; hmm.pi = pi;
			
					
		srand(123456789);
		
		int T = 10000;
		vector<double> O(T);
		hmm.generate_seq(O, T, false);
          
		hmm.init(); 
		hmm.print_iter = true;
		hmm.solver_eps = 1e-6;
		hmm.solver_max_iter = 1000;
		hmm.B = B;
		hmm.fit(O, true, 1e-5);
		cout << "A = "; hmm.printA();
		cout << "B = "; hmm.printB();
		cout << "pi = "; hmm.printpi();
			
	}//Pareto
	if( test_pareto ) {
		
		int N=2;
		int M=1;
		
		vector<vector<double> > A(N, vector<double>(N));
		vector<vector<double> > B(N, vector<double>(M));
		vector<double> pi(N);
		
		//A[0][0] = 1;
		A[0][0] = 0.9; A[0][1] = 0.1; 
		A[1][0] = 0.1; A[1][1] = 0.9;

		B[0][0] = 2; 
		B[1][0] = 1; 
		
		//pi[0] = 1;
		pi[0] = 0.6; pi[1] = 0.4;
		
		//cout << "Pareto Emission" << endl;

		paretoHMM<double> hmm(N,0,1000000); 
		hmm.init();

		hmm.A = A; hmm.B = B; hmm.pi = pi;
		hmm.lb[0] = 1;
		hmm.lb[1] = 1;
			
		int T = 10000;
		vector<double> O(T);
		hmm.generate_seq(O, T, false);
		
		hmm.init(); 
		hmm.fit(O, false, 1e-5);
		cout << "A = "; hmm.printA();
		cout << "B = "; hmm.printB();
		cout << "pi = "; hmm.printpi();
		
		hmm.init();
		hmm.fit(O, true, 1e-5);
		cout << "A = "; hmm.printA();
		cout << "B = "; hmm.printB();
		cout << "pi = "; hmm.printpi();
	}//multi emission
	if( test_multi ) {
		
		int N=3;
		int M=2;
		
		vector<vector<double> > A(N, vector<double>(N));
		vector<double> pi(N);
		
		//A[0][0] = 1;

		//A[0][0] = 0.9; A[0][1] = 0.1; 
		//A[1][0] = 0.1; A[1][1] = 0.9;
		
		A[0][0] = 0.9;  A[0][1] = 0.05;  A[0][2] = 0.05; 
		A[1][0] = 0.05; A[1][1] = 0.9;   A[1][2] = 0.05; 
		A[2][0] = 0.05; A[2][1] = 0.05;  A[2][2] = 0.9; 

		vector<vector<double> > Bexp(1, vector<double>(1));
		Bexp[0][0] = 100;  
		vector<vector<double> > Bnorm(1, vector<double>(2));
		Bnorm[0][0] = 3; Bnorm[0][1] = 0.1; 
		vector<vector<double> > Blognorm(1, vector<double>(2));
		Blognorm[0][0] = 5; Blognorm[0][1] = 1; 
			
		//pi[0]=1;
		//pi[0] = 0.6; pi[1] = 0.4;
		pi[0] = 0.6; pi[1] = 0.2; pi[2] = 0.2;
	
		exponentialHMM<double> hmm_exp(1,0,1000000); 
		gaussianHMM<double> hmm_norm(1,0,1000000); 
		lognormalHMM<double> hmm_lognorm(1,0,1000000);
		//vector< HMM<double>* > hmms = { &hmm_exp };
		//vector< HMM<double>* > hmms = { &hmm_exp, &hmm_norm };
		vector< HMM<double>* > hmms = { &hmm_exp, &hmm_norm, &hmm_lognorm };
		
		multiHMM<double> hmm(hmms,0,1000000); 
		hmm.init();
		hmm.A = A; hmm.pi = pi;
		hmm_exp.B = Bexp; hmm_norm.B = Bnorm; hmm_lognorm.B = Blognorm;

		//cout << "Joint Emission" << endl;

		int T = 30000;
		vector<double> O(T);
		hmm.generate_seq(O, T, false);
		
		
		hmm.init();
		hmm.fit(O, false, 1e-10);
		cout << "A = "; hmm.printA();
		cout << "B = "; hmm.printB();
		cout << "pi = "; hmm.printpi();
		
		hmm.init();
		hmm.fit(O, true, 1e-10);
		cout << "A = "; hmm.printA();
		cout << "B = "; hmm.printB();
		cout << "pi = "; hmm.printpi();
	}
	
	return 0;
	*/
	
} 


