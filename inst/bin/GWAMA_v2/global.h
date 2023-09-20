/*************************************************************************
GWAMA software:  May, 2010

Contributors:
    * Andrew P Morris amorris@well.ox.ac.uk
    * Reedik Magi reedik@well.ox.ac.uk

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
"AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*************************************************************************/
#ifndef _GLOBAL_H_
#define _GLOBAL_H_

#include <iostream>
#include <map>
#include <vector>
#include <stdlib.h>
#include <fstream>

using namespace std;
class filter
{
public:
    string _columnName;  // Column name for the filter check
    string _equation;    // >, <, >=, <=, !=, ==
    double _value;       // numerical value for filtering
    filter(string columnName, string equation, double value);            // constructor
    bool valid();             //consistency check
};
class global
{
public:
	int _errNR;		//count of errors and warnings found
	double _thresholdPValDir; //p-value threshold for directions
	bool _genomicControl;			//genomic control for input files
	bool _genomicControlOutput;			//genomic control for output
	bool _binaryTrait;		//by default binary trait
	bool _randomEffect;    //use random effect
	bool _genderHet;	//use gender-specific analysis
	bool _noAlleles;	//use gender-specific analysis
    bool _indel;        //use indel alleles
	vector <string> _alternativeMarkerNames;
	vector <string> _alternativeEffectAlleles;
	vector <string> _alternativeOtherAlleles;
	vector <string> _alternativeEffectAlleleFreqs;
	vector <string> _alternativeStrands;
	vector <string> _alternativeBetas;
	vector <string> _alternativeSEs;
	vector <string> _alternativeORs;
	vector <string> _alternativeOR_95Ls;
	vector <string> _alternativeOR_95Us;
	vector <string> _alternativeNs;
	vector <string> _alternativeImputeds;
	string _fileList;	//gwama.in file
	string _markerMap;   //map file name
	string _outputRoot;  //output file root (w/o file extension)
	string _version;

	int _markerCount;
	int _studyCount;
	vector <string> studyList;
	vector <string> studyGenders;
    vector <filter> filters;

	double outputLambda;
	double outputLambda_male;;
	double outputLambda_female;


	global();		//constructor
};



#endif

