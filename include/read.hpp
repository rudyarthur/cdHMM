#include <iostream>
#include <vector>
#include <algorithm>
#include <fstream>
#include <sstream>
#include <iterator>
#include <string>
#include <locale>         
#include <limits>
#include <ctype.h>
#include <math.h>
#include <set>
#include <map>

namespace cdHMM {
/*
map<char, int> charmap = {
{'a' , 0},
{'b' , 1},
{'c' , 2},
{'d' , 3},
{'e' , 4},
{'f' , 5},
{'g' , 6},
{'h' , 7},
{'i' , 8},
{'j' , 9},
{'k' , 10},
{'l' , 11},
{'m' , 12},
{'n' , 13},
{'o' , 14},
{'p' , 15},
{'q' , 16},
{'r' , 17},
{'s' , 18},
{'t' , 19},
{'u' , 20},
{'v' , 21},
{'w' , 22},
{'x' , 23},
{'y' , 24},
{'z' , 25},
{' ' , 26}
};

map<int, char> intmap = {
{0,'a'},
{1,'b'},
{2,'c'},
{3,'d'},
{4,'e'},
{5,'f'},
{6,'g'},
{7,'h'},
{8,'i'},
{9,'j'},
{10,'k'},
{11,'l'},
{12,'m'},
{13,'n'},
{14,'o'},
{15,'p'},
{16,'q'},
{17,'r'},
{18,'s'},
{19,'t'},
{20,'u'},
{21,'v'},
{22,'w'},
{23,'x'},
{24,'y'},
{25,'z'},
{26,' '}
};
*/
template <typename T> void read_txt(string filename, vector< T > &data, int col, string read_type, char delim = ' '){
	
	data.resize(0);
	ifstream in(filename.c_str());
	
	if(!in.is_open()){
		cout << "Failed to open file." << endl;
		exit(1);
	}
	
	string line = ""; 
	string item;
	
	//read rest of data
	while(getline(in,line)) // loop through the file
	{	
		
		stringstream parse(line); 
		vector<T> tmp;
		if( read_type == "double" ){
			while (getline(parse, item, delim)) { 
				if( item != "" ) tmp.push_back( stod(item) );
			} 
		}
		else if( read_type == "float" ){
			while (getline(parse, item, delim)) { 
				if( item != "" ) tmp.push_back( stof(item) );
			} 
		}
		else if( read_type == "int" ){
			while (getline(parse, item, delim)) { 
				if( item != "" ) tmp.push_back( stoi(item) );
			} 
		} else {
			cerr << "unknown numeric read type" << endl;  exit(1);
		}
		data.push_back(tmp[col]);
	}	
		
	in.close();
}

int read_words(string filename, vector< int > &data, map<char, int> &charmap, bool converttolower = false){

	set<char> charset;
	data.resize(0);
	ifstream in(filename.c_str());
	
	if(!in.is_open()){
		cout << "Failed to open file." << endl;
		exit(1);
	}
	
	char c;
	while(in.get(c)) // loop through the file character by character
	{	
		
		bool use=true;
		if(converttolower){
			use = false;
			if( isalpha(c) ){
				c = tolower(c); use = true;
			} else if(isspace(c) ){
				c = ' ';
				use = true;
			}
		}
		if( use ){
			int oldnum = charset.size();
			charset.insert( c );
			int idx;
			if( oldnum == charset.size() ){	//already added this char
				idx = charmap[c];
			} else { //new char
				charmap[c] = oldnum;
				idx = oldnum;
			}
			data.push_back( idx );
		}
		
	}	

	in.close();
	return charset.size();
}

}
