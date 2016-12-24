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
       
void usage(){
	cerr << "usage:\n";
	cerr << "./hmm -in type -model 1 ... input.txt\n";
	cerr << "Reading options\n";
	cerr << "\t--in char :: input is text\n";
	cerr << "\t--in int :: input is column of integers\n";
	cerr << "\t--in float :: input is column of floats\n";
	cerr << "\t--in double :: input is column of doubles\n";
	cerr << "\t--col c :: read column c\n";
	cerr << "\t--tolower :: skips all punctuation except space, converts to lower case\n";
	cerr << "Possibilities for -model\n";
	cerr << "\t--gaussian N :: N states with gaussian emission\n";
	cerr << "\t--exponential N :: N states with exponential emission\n";
	cerr << "\t--gamma N :: N states with gamma emission\n";
	cerr << "\t--lognormal N :: N states with lognormal emission\n";
	cerr << "\t--pareto N :: N states with pareto emission\n";
	cerr << "\t--laplace N :: N states with laplace emission\n";
	cerr << "\t--discrete N :: N states with discrete emission\n";
	cerr << "\t--discreteobs M :: M emission types in each discrete state (coded 0, 1, 2, ...)\n";
	cerr << "\t--poisson N :: N states with poisson emission\n";
	cerr << "\t--gumbel N :: N states with gumbel emission\n";
	cerr << "\t--weibull N :: N states with weibull emission\n";
	cerr << "\t--extremevalue N :: N states with extremevalue emission\n";
	cerr << "Fitting options\n";
	cerr << "\t--restarts :: number of HMM runs from random (default 1)\n";
	cerr << "\t--precision:: HMM stopping precision (default 1e-5)\n";
	cerr << "\t--min_hmm_iters :: maximum number of hmm iterations (default 0)\n";
	cerr << "\t--max_hmm_iters :: minimum number of hmm iterations (default 1000)\n";
	cerr << "\t--generate :: generate this many observations from the HMM(default 0)\n";
	exit(1);
}

template <typename T> void run(map<string, int> &states, vector<string> &model,
double eps, int restarts, int min_hmm_iters, int max_hmm_iters, int generate,
vector<T> &data, map<char, int> &charmap, string read_type, int M=1
){

		vector< HMM<T>* > hmms;
		for(int i=0; i<model.size(); ++i){
			HMM<T>* modelhmm = HMMFactory<T>::generateHMM(model[i], states[model[i]], 0, 1000, M);
			hmms.push_back( modelhmm );
		}

		multiHMM<T> hmm(hmms,min_hmm_iters,max_hmm_iters);
		hmm.info();
		//Fit that observation sequence
		cerr << "Running";
		for(int i=0; i<restarts; ++i){
			cerr << ".";
			hmm.init();
			hmm.fit(data, eps);
			hmm.sortparams();
		}
		cerr << endl;

		cout << fixed << "Log Likelihood: " << hmm.maxlhood << endl;
		cout << fixed << "Transition Matrix" << endl;
		for(int i=0; i<hmm.N; ++i){
			for(int j=0; j<hmm.N; ++j){
				cout << hmm.maxA[i][j] << "\t";
			} cout << endl;
		}
		cout << "Emission Matrix" << endl;
		if(read_type == "char"){
			for(int j=0; j<hmm.M; ++j){
				//bad
				map<char,int>::iterator it;
				for(it = charmap.begin(); it!=charmap.end(); ++it){
					if( it->second == j ) break;
				}
				if( it != charmap.end() ){ 
					cout << it->first << ":\t"; 
				}
				for(int i=0; i<hmm.N-1; ++i){
					cout << hmm.maxB[i][j] << "\t";
				} cout << hmm.maxB[hmm.N-1][j] << endl;
			}
			cout << endl; 
		} else {
			for(int i=0; i<hmm.N; ++i){
				for(int j=0; j<hmm.M; ++j){
					cout << hmm.maxA[i][j] << "\t";
				} cout << endl;
			}
		}
		if( generate > 0){
			cout << "Synthethic data" << endl;
			vector<T> O;
			hmm.generate_seq(O, generate, true);
		}
		
		for(int i=0; i<hmms.size(); ++i){ delete hmms[i]; }
		
		
}

int main(int argc,char **argv){
	
	if(argc<3) usage();
	 
	int c;
	static struct option loptions[] =    {
        {"in",1,0,'i'},
        {"col",1,0,'c'},
        {"tolower",0,0,'j'},
        {"gaussian",1,0,'g'},
        {"exponential",1,0,'e'},
        {"gamma",1,0,'m'},
        {"lognormal",1,0,'l'},
        {"pareto",1,0,'p'},
        {"laplace",1,0,'a'},
        {"discrete",1,0,'d'},
        {"discreteobs",1,0,'1'},
        {"poisson",1,0,'o'},
        {"gumbel",1,0,'b'},
        {"weibull",1,0,'w'},
        {"extremevalue",1,0,'v'},
        {"restarts",1,0,'r'},
        {"precision",1,0,'q'},
        {"min_hmm_iters",1,0,'n'},
        {"max_hmm_iters",1,0,'x'},
        {"generate",1,0,'y'},
        {0,0,0,0}
    };
	
	string read_type = "none";
	int col = 0;
	
	map<string, int> states;
	vector<string> model;
	 
	double eps = 1e-5;
	int restarts = 1;
	int min_hmm_iters = 0;
	int max_hmm_iters = 1000;
	int generate = 0;
	bool converttolower = false;
	int M=-1;
    while ((c = getopt_long(argc, argv, "i:c:g:e:m:l:p:a:d:1:o:b:w:v:r:q:n:x:y:j",loptions,NULL)) >= 0) {  
        switch (c)
        {
        case 'i': read_type = (string)optarg; break;
        case 'c': col = atoi(optarg); break;
        case 'j': converttolower = true; break;
        case 'g': states[ "gaussian" ] = atoi(optarg); model.push_back("gaussian"); break;
        case 'e': states[ "exponential" ] = atoi(optarg); model.push_back("exponential"); break;
        case 'm': states[ "gamma" ] = atoi(optarg); model.push_back("gamma"); break;
        case 'l': states[ "lognormal" ] = atoi(optarg); model.push_back("lognormal"); break;
        case 'p': states[ "pareto" ] = atoi(optarg); model.push_back("pareto"); break;
        case 'a': states[ "laplace" ] = atoi(optarg); model.push_back("laplace"); break;
        case 'd': states[ "discrete" ] = atoi(optarg); model.push_back("discrete"); break;
        case '1': M = atoi(optarg); break;
        case 'o': states[ "poisson" ] = atoi(optarg); model.push_back("poisson"); break;
        case 'b': states[ "gumbel" ] = atoi(optarg); model.push_back("gumbel"); break;
        case 'w': states[ "weibull" ] = atoi(optarg); model.push_back("weibull"); break;
        case 'v': states[ "extremevalue" ] = atoi(optarg); model.push_back("extremevalue"); break;
        case 'r': restarts = atoi(optarg); break;
        case 'q': eps = atof(optarg); break;
        case 'n': min_hmm_iters = atoi(optarg); break;
        case 'x': max_hmm_iters = atoi(optarg); break;
        case 'y': generate = atoi(optarg); break;
        case '?': usage();
        default: 
			usage();
        }
    }
    if( read_type == "none" ){ cerr << "--in argument is mandatory!" << endl; usage(); }
    if( model.size() == 0 ){ cerr << "specify at least one emission model!" << endl; usage(); }
    
    string input = argv[optind];
	cerr <<"Parsing " << input << " as " << read_type << endl; 
	
	map<char, int> charmap;
	if( read_type == "double" ){
		vector<double> data;
		read_txt(input, data, col, read_type);
		cerr << "Read " << data.size() << " values." << endl;
		
		run(states, model, eps, restarts, min_hmm_iters, max_hmm_iters, generate, data, charmap, read_type, M);
	}
	if( read_type == "float" ){
		vector<float> data;
		read_txt(input, data, col, read_type);
		cerr << "Read " << data.size() << " values." << endl;
		
		run(states, model, eps, restarts, min_hmm_iters, max_hmm_iters, generate, data, charmap, read_type, M);
	}
	if( read_type == "int" ){
		vector<int> data;
		read_txt(input, data, col, read_type);
		cerr << "Read " << data.size() << " values." << endl;
		
		run(states, model, eps, restarts, min_hmm_iters, max_hmm_iters, generate, data, charmap, read_type, M);
	}
	if( read_type == "char" ){
		vector<int> data;
		M = read_words(input, data, charmap, converttolower);
		cerr << "Read " << data.size() << " values." << endl;
		cerr << "Read " << M << " unique chars." << endl;
		
		if(model.size() > 1 || model[0] != "discrete"){
			cerr << "char option only works with discrete HMM" << endl;
		}
		run(states, model, eps, restarts, min_hmm_iters, max_hmm_iters, generate, data, charmap, read_type, M);
	}
	
	
	
	return 0;
	
	
} 


