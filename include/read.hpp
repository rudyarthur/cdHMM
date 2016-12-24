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
 /*!
  * Read a column of numbers
  * @param[in] filename Name of the input file
  * @param[out] data Vector containing the data read
  * @param[in] col Column of the file that is to be read
  * @param[in] read_type Type of data to be read
  * @param[in] delim Separator of the columns. Default whitespace.
  * */
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
 /*!
  * Read a test file
  * @param[in] filename Name of the input file
  * @param[out] data Vector containing the data read (as int)
  * @param[in] charmap map of characters to integers
  * @param[in] converttolower ignore pnctuation and capitalization
  * */
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
