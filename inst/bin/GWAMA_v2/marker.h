/*************************************************************************
GWAMA software:  May, 2009

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


#pragma once
#include <iostream>
#include <map>
#include <vector>
#include <stdlib.h>
#include <fstream>
#include "study.h"
#include "global.h"
#include "problem.h"

using namespace std;

class marker
{

private:
	std::string _name;

	std::string _direction;
	string _effectAllele;
	string _otherAllele;
	int _chromosome;
	int _position;
	int _N;
	bool _Nok; //if any population is missing N, then it turns false
	double _EAF;

	double _B;
	double _W;
	double _W2;
	double _Q;
	double _Br;
	double _Wr;
	double _Qr;

	//gender specific 
	int _Nmale;
	double _EAFmale;

	double _Bmale;
	double _Wmale;
	double _W2male;
	double _Qmale;
	double _Brmale;
	double _Wrmale;
	double _Qrmale;

	int _Nfemale;
	double _EAFfemale;

	double _Bfemale;
	double _Wfemale;
	double _W2female;
	double _Qfemale;
	double _Brfemale;
	double _Wrfemale;
	double _Qrfemale;

	vector <double> all_se;
	vector <double> male_se;
	vector <double> female_se;
	vector <double> all_beta;
	vector <double> male_beta;
	vector <double> female_beta;

/*
private:
	std::vector <bool> _imputed;
	std::vector <double> _effectAlleleFreq;
	std::vector <char> _effectAlleleVec;
	std::vector <char> _nonEffectAlleleVec;
	std::vector <double> _beta;
	std::vector <double> _se;
	std::vector <int> _studyNr;
	std::vector <int> _n;
*/

public:
	int _studyCount;
	int male_studyCount;
	int female_studyCount;

	marker(std::string name, int popCount);
//	int addStudy(bool imp, char strand, char effectAllele, char nonEffectAllele, double effectAlleleFreq, 
//		double beta, double se, int studyNr,vector <study> & studies,int n, ofstream & ERR, int & errNR, ofstream & LOG);
	int addCohort(int fileNr, global & g, double directLambda, double imputedLambda, 
		string myMarker,string myEffectAllele, string myOtherAllele, double myEffectAlleleFreq,
		bool myStrand, double myBeta, double mySE, bool myImputed, int myN, ofstream & ERR, ofstream & LOG, problem & markerProblems);

	bool addMap(int chromosome, int position);
	bool calculateBWQ(vector <study> &, bool useLambda,bool useRandom, double threshold ,ofstream & ERR, int & errNR, ofstream & LOG);
	bool doRandom(global & g);

	double sumBeta(bool useRandom);
	double sumSE(bool useRandom);
	double Lower95Beta(bool useRandom);
	double Upper95Beta(bool useRandom);
	double OR(bool useRandom);
	double Lower95OR(bool useRandom);
	double Upper95OR(bool useRandom);
	double pValue(double outputLambda,bool useRandom);
	double genderdiff_pValue(global & g,bool useRandom);
	double genderhet_pValue(global & g,bool useRandom);
	double male_sumBeta(bool useRandom);
	double male_sumSE(bool useRandom);
	double male_Lower95Beta(bool useRandom);
	double male_Upper95Beta(bool useRandom);
	double male_OR(bool useRandom);
	double male_Lower95OR(bool useRandom);
	double male_Upper95OR(bool useRandom);
	double male_pValue(global & g,bool useRandom);
double male_z(bool useRandom);
	double female_sumBeta(bool useRandom);
	double female_sumSE(bool useRandom);
	double female_Lower95Beta(bool useRandom);
	double female_Upper95Beta(bool useRandom);
	double female_OR(bool useRandom);
	double female_Lower95OR(bool useRandom);
	double female_Upper95OR(bool useRandom);
	double female_pValue(global & g,bool useRandom);
double female_z(bool useRandom);

	double z(bool useRandom);
	double qStatistic();
	double qPValue();
	double i2();
	std::string printMarker(global & g);

	~marker(void);
};
