#include <iostream>
#include <vector>
#include <map>
#include <fstream>
#include <sstream>
#include <iterator>
#include <string>
#include <algorithm>

using namespace std;

bool pair_compare_first ( const pair<double, int>& l, const pair<double, int>& r){ return l.first < r.first; }


class patient {
public:
	string syn;
	int num;
	int id;
	vector<int> type;
	vector<double> time;
	vector<double> duration;
	double wake;
	double sleep;
	
	patient(){
			id = 0;
			syn = "";
			num = 0;
			type.resize(0);
			time.resize(0);
			duration.resize(0);
			wake = 0;
			sleep = 0;
	}
	
	void set_id( int id_ ){ id = id_; };
	void set_syn( string syn_ ){ syn = syn_; };
	void add_type( int t ){ type.push_back(t); };
	void add_time( double t ){ time.push_back(t); };
	void add_duration( double t ){ duration.push_back(t); };
	void add_wake( double t ){ wake = (t); };
	void add_sleep( double t ){ sleep = (t); };
    void set_num(){ num = type.size(); }
    int get_num(){ set_num(); return num; }
    
    void print(){
		set_num();
		cout << "ID = " << id << "\n"
		<< "Syndrome = " << syn << "\n"
		<< num << " records" << endl;
	}
	
	void full_print(){
		print();
		cout << "Wake: " << wake << "\n";
		cout << "Sleep: " << sleep << "\n";
		cout << "Times: "; for(int i=0; i<time.size(); ++i){ cout << time[i] << " ";} cout << "\n";
		cout << "Types: "; for(int i=0; i<type.size(); ++i){ cout << type[i] << " ";} cout << "\n";
		cout << "Durations: "; for(int i=0; i<duration.size(); ++i){ cout << duration[i] << " ";} cout << endl;
	}
	
	
	void gen_obs(int T, vector<int> &obs, double long_dur = numeric_limits<double>::infinity() ){

		if( time.size() != duration.size() ){ cerr << "#time != #duration" << endl; exit(1); }
		
		obs.resize(T); for(int t=0; t<T; ++t){ obs[t] = 0; }

		//convert time
		for(int i=0; i<time.size(); ++i){ 
			
			for( int j=(int)( time[i] * (double)T ); 
			j<=(int)( (time[i] * (double)T) + (duration[i]* (double)T/(86400.0)) ); ++j ){

				if( duration[i] < long_dur ){
					obs[ j ] = 1;
				} else {
					obs[ j ] = 2;
				}
				
			}
		}

	}	

	void gen_obs_bin(int T, vector<int> &obs, double long_dur = numeric_limits<double>::infinity() ){

		if( time.size() != duration.size() ){ cerr << "#time != #duration" << endl; exit(1); }
		
		obs.resize(T); for(int t=0; t<T; ++t){ obs[t] = 0; }

		//convert time
		for(int i=0; i<time.size(); ++i){ 
			
			for( int j=(int)( time[i] * (double)T ); 
			j<=(int)( (time[i] * (double)T) + (duration[i]* (double)T/(86400.0)) ); ++j ){

					++obs[ j ];
				
			}
		}

	}	
	
	void sort_times(){
		vector< pair<double, int> > slist(time.size()); 
		for(int i=0; i<time.size(); ++i){ slist[i].first = time[i]; slist[i].second = i; } 
		sort(slist.begin(), slist.end(), pair_compare_first);
		
		vector<double> tmp( time.size() );
		for(int i=0; i<time.size(); ++i){ tmp[ i ] = time[ slist[i].second ]; } time = tmp;
		for(int i=0; i<time.size(); ++i){ tmp[ i ] = duration[ slist[i].second ];  } duration = tmp;
		vector<int> tmpi(time.size() );
		for(int i=0; i<time.size(); ++i){ tmpi[ i ] = type[ slist[i].second ];  } type = tmpi;
	}
	
	int interval_count(double start, double end){

		int st=-1;
		for(int i=0; i<time.size(); ++i){ if(time[i] >= start){st = i; break;} }
		int nd = -1;
		for(int i=0; i<time.size(); ++i){ if(time[i] >= end){nd = i; break;} }

		//periodic count
		return (nd - st + time.size())%time.size();
	}
	
	int wake_shift(int T){ return (int)( wake * (double) T ); }
	void wake_shift(int T, vector<int> &obs){ 
		vector<int> tmp;
		for(int t=wake_shift(T); tmp.size() < obs.size(); ++t){ tmp.push_back( obs[ t%T ] ); }
		obs.resize(tmp.size()); for(int t=0; t<tmp.size(); ++t){ obs[t] = tmp[t]; }
	}
	int sleep_shift(int T){ return (int)( sleep * (double) T ); }
	void sleep_shift(int T, vector<int> &obs){ 
		vector<int> tmp;
		for(int t=sleep_shift(T); tmp.size() < obs.size(); ++t){ tmp.push_back( obs[ t%T ] ); }
		obs.resize(tmp.size()); for(int t=0; t<tmp.size(); ++t){ obs[t] = tmp[t]; }
	}
	
	void awake(int T, vector<int> &obs){ //periodic!!!
		vector<int> tmp;
		for(int t=wake_shift(T); t%T != sleep_shift(T); ++t){ tmp.push_back( obs[ t%T ] ); }
		obs.resize(tmp.size()); for(int t=0; t<tmp.size(); ++t){ obs[t] = tmp[t]; }
	}
	
	void asleep(int T, vector<int> &obs){ //periodic!!!
		vector<int> tmp;
		for(int t=sleep_shift(T); t%T!=wake_shift(T); ++t){ tmp.push_back( obs[ t%T ] ); }
		obs.resize(tmp.size()); for(int t=0; t<tmp.size(); ++t){ obs[t] = tmp[t]; }
	}

};

void read_csv(string filename, map<int, patient> &data, char delim = ','){
	
	ifstream in(filename.c_str());
	
	if(!in.is_open()){
		cout << "Failed to open file." << endl;
		exit(1);
	}
	
	string line = ""; 
	string item;
	map< string, int > key;

	//first line -> key
	getline(in,line);
	stringstream ss(line); 
	int id = 0; while (getline(ss, item, delim)) { key[ item ] = id++; }

	//read rest of data
	while(getline(in,line)) // loop through the file
	{	
		
		stringstream parse(line); 
		id = 0;
		int pid = 0;
		while (getline(parse, item, delim)) { 
			
			if(id == key["UniqueID"]){
				pid = stoi(item);
				if( data.find( pid ) == data.end() ){	//new ID	
					patient p; p.set_id( pid );
					data[ pid ] = p;
				} 
			}
			
			if(id == key[ "Syndrome" ]){ data[ pid ].set_syn( item ); }
			if(id == key[ "discharge_time11" ]){ item.pop_back(); data[ pid ].add_time( stod(item)/100.0 ); } //% version
			if(id == key[ "discharge_duration11" ]){ data[ pid ].add_duration( stod(item) ); }
			if(id == key[ "discharge_type11" ]){ data[ pid ].add_type( stoi(item) ); }
			if(id == key[ "Time_Sleep_End"] ){ item.pop_back(); data[ pid ].add_wake( stod(item)/100.0 ); }
			if(id == key[ "Time_Sleep_Start"] ){ item.pop_back(); data[ pid ].add_sleep( stod(item)/100.0 ); }
			
			id++; 
			
		}
	}	
	
	/*int ct = 0;
	for(map<int, patient>::iterator iter = data.begin(); iter != data.end(); ++iter){
		iter->second.full_print();
		++ct;
	}
	cout << "tot = " << ct << endl;*/
		
	in.close();
}



void read_txt(string filename, vector< vector<double> > &data, char delim = ' '){
	
	data.resize(0);
	ifstream in(filename.c_str());
	
	if(!in.is_open()){
		cout << "Failed to open file." << endl;
		exit(1);
	}
	
	string line = ""; 
	string item;
	int id;
	
	//read rest of data
	while(getline(in,line)) // loop through the file
	{	
		
		stringstream parse(line); 
		id = 0;
		vector<double> tmp;
		while (getline(parse, item, delim)) { 
			if( tmp.size() < 2){
				tmp.push_back( stod(item) );
			} else {
				tmp.push_back( stod(item) );
			}
		}
		data.push_back(tmp);
	}	
		
	in.close();
}
