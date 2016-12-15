#pragma once

#include <iostream>
#include <vector>
#include <map>
#include <fstream>
#include <sstream>
#include <iterator>
#include <string>
#include <algorithm>
#include <iomanip>      
#include "hmm.hpp"

using namespace std;

class boot {
public:
	int Nboot;
	vector<double> data;
	
	double av;
	double err;
	double uerr;
	double lerr;
	bool calcdone;
	bool shuffle;
	
	double confidence;
	
	boot(){
		Nboot = 100;
		data.resize(Nboot);
		calcdone = false;
		shuffle = false;
		confidence = 0.6827;
	}
	
	boot(int Nboot_){
		Nboot = Nboot_;
		data.resize(Nboot);
		calcdone = false;
		shuffle = false;
		confidence = 0.6827;
	}
	
	boot(vector<double> &data_){
		Nboot = data_.size();
		data.resize(Nboot);
		for(int i=0; i<Nboot; ++i){
			data[i] = data_[i];
		}
		calcdone = false;
		shuffle = false;
		confidence = 0.6827;
	}
	
	friend void swap(boot& first, boot& second)
    {
        swap(first.Nboot, second.Nboot);
        swap(first.data, second.data);
    }	
    boot& operator=(boot other)
	{
		swap(*this, other); 
		return *this;
	}
	
	boot& operator+=(const boot &rhs) {
		Nboot = rhs.Nboot;
		data.resize(Nboot);
		if( shuffle ){
			vector<double> shuf = rhs.data;
			random_shuffle ( shuf.begin(), shuf.end() );
			for(int i=0; i<Nboot; ++i){
				data[i] += shuf[i];
			}
		} else {
			for(int i=0; i<Nboot; ++i){
				data[i] += rhs.data[i];
			}
		}
		
		return *this;
    }
	const boot operator+(const boot &other) const {
		return boot(*this) += other;
    }
    
    boot& operator*=(const boot &rhs) {
		Nboot = rhs.Nboot;
		data.resize(Nboot);
		if( shuffle ){
			vector<double> shuf = rhs.data;
			random_shuffle ( shuf.begin(), shuf.end() );
			for(int i=0; i<Nboot; ++i){
				data[i] *= shuf[i];
			}
		} else {
			for(int i=0; i<Nboot; ++i){
				data[i] *= rhs.data[i];
			}
		}
		return *this;
    }
	const boot operator*(const boot &other) const {
		return boot(*this) *= other;
    }
   
	void calc(double av_){
		
		vector<double> sorted_data = data;
		sort( sorted_data.begin(), sorted_data.end() );
		
		av = av_;
		/*av = 0;
		for(int b=0; b<sorted_data.size(); ++b){ 
			av += sorted_data[b]; 
		} av /= sorted_data.size();*/
		/*err = 0;
		for(int b=0; b<sorted_data.size(); ++b){ 
			err += ( sorted_data[b] - av )*( sorted_data[b] - av ); 
		} err /= (sorted_data.size()-1);
		err = sqrt(err);*/
		
		int lt=0, gt=0;
		for(int b=0; b<sorted_data.size(); ++b){ 
			if( sorted_data[b] < av ){ ++lt; }
			else{ ++gt; } 
			//cout<<sorted_data[b]<<endl;
		}
		//if symm lt ~ gt
		//toss = 0.5*(1.0 - confidence);
		//left_boundary = toss * (2lt)
		//right_boundary = B - toss * (2gt)
		
		double toss = 0.5*(1.0 - confidence);
		double lb = sorted_data[ (int)(toss*2*lt) ];
		double ub = sorted_data[ (int)(sorted_data.size() - 1 - toss*2*gt) ];

		lerr = (av-lb);
		uerr = (ub-av);
		
		//cerr << av << " " << lt << " " << gt << " " << lb << " " << ub << endl; exit(1);
		
		calcdone = true;
	}
	
	void pstream(stringstream& os) const {
		os << fixed << setw(4) << setprecision(4) << av << " \\pm _{" << lerr << "}^{" << uerr << "}";
	}

	friend ostream& operator<<(std::ostream& os, const boot& b){
		std::stringstream ss;
		b.pstream(ss);
		os << ss.str();
	  	return os;
	};
	
	
};

void contract(vector<int> &in, vector<int> &out, int w){
	out.resize(0);
	int t=0;

	for(t=0; t<in.size(); t+=w){
		int obs = 0;
		for(int i=0; i<w; ++i){
			obs += in[t + i];
		}
		/*if(obs==0){ out.push_back(0); }
		else{ out.push_back(1); }*/
		out.push_back(obs);
	}
	/*if( (in.size()%w)!=0 ){
		t -= w;
		int obs = 0;
		for(int i=t; i<in.size(); ++i){
			obs += in[i];
		}
		if(obs==0){ out.push_back(0); }
		else{ out.push_back(1); }
	}*/
}

/*void integrateB( Ref<MatrixXd> B, Ref<MatrixXd> IB, vector<vector< boot> > &bootB, vector<vector< boot> > &IbootB){

  int N = B.rows();
  for(int i=0; i<N; ++i){
	IB(i,0) = B(i,0)*B(i,0);
	IB(i,1) = 2*B(i,0)*B(i,1) + B(i,1)*B(i,1);
  }

  IbootB.resize(N);
  for(int i=0; i<N; ++i){
	boot tmp = bootB[i][0]*bootB[i][0]; tmp.calc(bootB[i][0].av*bootB[i][0].av);
	IbootB[i].push_back(tmp);
	boot tmp2 = bootB[i][1]*bootB[i][0]+bootB[i][0]*bootB[i][1]+bootB[i][1]*bootB[i][1]; 
	tmp2.calc(2*bootB[i][0].av*bootB[i][1].av + bootB[i][1].av*bootB[i][1].av);
	IbootB[i].push_back(tmp2);
  }
  
}*/

void print( vector<vector< boot> > &B, int p_=4 ){
	for(int i=0; i<B.size(); ++i){ 
		for(int j=0; j<B[i].size()-1; ++j){ 
			cout << fixed << setw(p_) << setprecision(p_) << B[i][j] << " &"; 
		} cout << fixed << setw(p_) << setprecision(p_) << B[i][B[i].size()-1] << " \\\\\n"; 
	} 
}

void graphBoot( vector<vector< boot> > &A, vector<vector< boot> > &B, int p_=4 ){
	cout << "digraph G {\n";
	for(int i=0; i<A.size(); ++i){ 
		cout << "\tState" << i << ";\n";
	}
	
	for(int j=0; j<B[0].size(); ++j){ 
		//cout << "\tObs" << j << " [ shape=box ]; \n";
		
		cout << "{ rank=sink; ";
		if( j==0 ){ cout << "\tNormal [ shape=box ];"; }
		else if ( j==1 ){ cout << "\tSpike [ shape=box ];"; }
		else if ( j==2 ){ cout << "\tSeizure [ shape=box ];"; }
		cout << " }\n";
	}
	
	for(int i=0; i<A.size(); ++i){ 
		for(int j=0; j<A[i].size(); ++j){ 
			if( A[i][j].av - A[i][j].lerr > 1e-5 ){
			cout << "\tState" << i << " -> State" << j;
			cout << " [ label= \"" << A[i][j].av << "\" ];\n"; }
		} 
	} 
	
	for(int i=0; i<B.size(); ++i){ 
		for(int j=0; j<B[i].size(); ++j){ 
			if( B[i][j].av - B[i][j].lerr > 1e-5 ){
			//cout << "\tState" << i << " -> Obs" << j;
			cout << "\tState" << i;
			if( j==0 ){ cout << " -> Normal"; }
			else if ( j==1 ){ cout << " -> Spike"; }
			else if ( j==2 ){ cout << " -> Seizure"; }
			cout << " [ style=dotted, label= \"" << B[i][j].av << "\" ];\n"; }
		} 
	} 
	
	cout << "}\n";
}

template<typename number> void generate_bootstrap( HMM<number> &hmm_eval, vector<vector< boot> > &bA, vector<vector< boot> > &bB , int T, int Nboot, double eps){

	int N = hmm_eval.N;
	int M = hmm_eval.M;
		
	boot empty(Nboot);
	bA.resize(N); for(int i=0; i<N; ++i){ bA[i].resize(N); for(int j=0; j<N; ++j){ bA[i][j] = empty; } }
	bB.resize(N); for(int i=0; i<N; ++i){ bB[i].resize(M); for(int j=0; j<M; ++j){ bB[i][j] = empty; } }
	
	vector<vector<number> > At;
	vector<vector<number> > Bt;
	vector<number> pit;
  
	for(int b=0; b<Nboot; ++b){
		cerr << b << endl;
		vector<int> newO;
		hmm_eval.generate_seq(newO, T);

		vector<number> res_best(2, -numeric_limits<number>::infinity() );
		int repeat = 1;
		for(int r=0;r<repeat;++r){
			HMM<number> hmm_boot(N,M,0,1000000);
			hmm_boot.init();
			
			hmm_boot.A  = hmm_eval.A;
			hmm_boot.B  = hmm_eval.B;
			hmm_boot.pi = hmm_eval.pi;  
		
			//hmm_boot.fuzz(0.1);
			
			vector<number> res = hmm_boot.fit( newO , eps );	
			if( res[1] > 0 ){
				cout << "big err " << endl;
				cout << res[0] << " " << res[1] << endl;
				exit(1);
			}
			hmm_boot.sortparams();
			if( res[1] > res_best[1] ){
				res_best[0] = res[0]; res_best[1] = res[1];
				At = hmm_boot.A; Bt = hmm_boot.B; pit = hmm_boot.pi; 
			}
		}	
	
		for(int i=0; i<N; ++i){ for(int j=0; j<N; ++j){ bA[i][j].data[b] = At[i][j]; } }
		for(int i=0; i<N; ++i){ for(int j=0; j<M; ++j){ bB[i][j].data[b] = Bt[i][j]; } }

	}

	for(int i=0; i<N; ++i){ for(int j=0; j<N; ++j){ bA[i][j].calc( hmm_eval.A[i][j] ); } } 
    for(int i=0; i<N; ++i){ for(int j=0; j<M; ++j){ bB[i][j].calc( hmm_eval.B[i][j] ); } } 
	
}
