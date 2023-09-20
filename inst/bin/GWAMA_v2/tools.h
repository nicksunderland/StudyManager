#ifndef _TOOLS_H_
#define _TOOLS_H_

#include <iostream>
#include <map>
#include <vector>
#include <stdlib.h>
#include <fstream>
#include <math.h>
#include <cctype> // std::toupper
#include <string>
#include <algorithm>
using namespace std;

void sortVec(vector <double>& x, int size);

int Tokenize(const  string& str1,
                      vector<string>& tokens,
			const string& delimiters);
string uc(string s);	//uppercase
bool checkAlleles(string & s1, string & s2);	//check if alleles are ok and change numbers to letters if necessary
string flip(string s);	//flip the alleles if 





#endif
