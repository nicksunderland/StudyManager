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


#include <iostream>
#include <map>
#include <vector>
#include <stdlib.h>
#include <fstream>
#include "marker.h"
#include "study.h"
#include <math.h>
#include <cctype> // std::toupper
#include <zlib.h>
#include "global.h"
#include "commandLine.h"
#include "readFile.h"
#include "tools.h"

using namespace std;





int 
main(int argc, char *argv[]) 
{
	string gwamaversion = "2.1";
	global _g;

	char sb [11];
	vector <string> fileGenders;
	vector <string> snpNames;
	vector <study> studies;
	bool proc;

	proc = readCommandline( argc, argv, _g);
	if (!proc){cerr << "Error reading command line options! Exit program." << endl; exit(1);}
	string gwamagc = _g._outputRoot + ".gc.out";
	string gwamaout = _g._outputRoot + + ".out";
	string gwamaerr = _g._outputRoot + + ".err.out";
	string gwamalog = _g._outputRoot + + ".log.out";
	ofstream LAMBDA_FILE (gwamagc.c_str());
	ofstream ERR (gwamaerr.c_str());
	ofstream LOG (gwamalog.c_str());

	LAMBDA_FILE << "Cohort\tdirectly_genotyped_markers_lambda\t\tdirectly_genotyped_markers_count\timputed_markers_lambda\timputed_markers_count\n";
	LOG << "Running GWAMA " << gwamaversion << endl;
	proc = readFilelist(_g, ERR,  LOG);
	if (!proc){cerr << "Error reading file list! Exit program." << endl; exit(1);}

	map <string, int> markerPosition;
	vector <marker> markerList;
	for (int i = 0; i < _g._studyCount; i++)
	{
		cout << "------------------------------------\nReading file: " << _g.studyList[i] << endl;
		if (_g._genderHet && _g.studyGenders.size()>i) cout << "Gender status: " << _g.studyGenders[i] << endl;
		else if (_g._genderHet){ cerr << "Error reading gender status from filelist. Exit program<< "; exit(1);}
		proc = readCohort(i,_g,markerPosition,markerList,ERR,LOG, LAMBDA_FILE);
	}
	ofstream OUT (gwamaout.c_str());

	if (_g._markerCount==1) {cerr << "No markers found! Exit program!"<<endl; exit(1);}

	if (_g._markerMap != "N")
	{
		cout << "------------------------------------\nReading map file"<<endl;
		proc = readMapFile(_g._markerMap,markerList,markerPosition, ERR, _g._errNR, LOG);
	}
	cout << "------------------------------------\nPreparing output..."<<endl;
	if (_g._genomicControlOutput)
	{
		vector <double> _chiStat;
		vector <double> _chiStat_male;
		vector <double> _chiStat_female;
		for (int i = 0; i < _g._markerCount-1; i++)
		{
			 if (markerList[i]._studyCount>0)
			 {_chiStat.push_back((markerList[i].sumBeta(0)*markerList[i].sumBeta(0))/(markerList[i].sumSE(0)*markerList[i].sumSE(0)));}
			if (_g._genderHet)
			{
				if (markerList[i].male_studyCount>0)
				{_chiStat_male.push_back((markerList[i].male_sumBeta(0)*markerList[i].male_sumBeta(0))/(markerList[i].male_sumSE(0)*markerList[i].male_sumSE(0)));}
				if (markerList[i].female_studyCount>0)
				{_chiStat_female.push_back((markerList[i].female_sumBeta(0)*markerList[i].female_sumBeta(0))/(markerList[i].female_sumSE(0)*markerList[i].female_sumSE(0)));}

			}
		}
	

		double _median = 0;
		sortVec( _chiStat, _chiStat.size());
		if ((_chiStat.size())%2!=0) _median = _chiStat[((_chiStat.size())/2)-1];
		else _median = (_chiStat[(_chiStat.size())/2 -1] + _chiStat[(_chiStat.size())/2])/2;
		_g.outputLambda = _median/0.4549364;		//median of chi-sq from R ... qchisq(0.5, df= 1)
		if (_g.outputLambda>1)LOG << "Output p-value is corrected with lambda=" << _g.outputLambda << endl;
		else LOG << "Output p-value is not corrected. Lambda=" << _g.outputLambda << endl;
		if (_g.outputLambda>1)cout <<"------------------------------------\nOutput p-value is corrected with lambda=" << _g.outputLambda << endl;
		else cout <<"------------------------------------\nOutput p-value is not corrected. Lambda=" << _g.outputLambda << endl;
			if (_g._genderHet)
			{
				_median = 0;
				sortVec( _chiStat_male, _chiStat_male.size());
				if ((_chiStat_male.size())%2!=0) _median = _chiStat_male[((_chiStat_male.size())/2)-1];
				else _median = (_chiStat_male[(_chiStat_male.size()-1)/2 -1] + _chiStat_male[(_chiStat_male.size()-1)/2])/2;
				_g.outputLambda_male = _median/0.4549364;		//median of chi-sq from R ... qchisq(0.5, df= 1)
				if (_g.outputLambda_male>1)LOG << "Output p-value is corrected with lambda=" << _g.outputLambda_male << endl;
				else LOG << "Output p-value is not corrected. Lambda=" << _g.outputLambda_male << endl;
				if (_g.outputLambda_male>1)cout <<"------------------------------------\nOutput p-value is corrected with lambda=" << _g.outputLambda_male << endl;
				else cout <<"------------------------------------\nOutput p-value is not corrected. Lambda=" << _g.outputLambda_male << endl;

				_median = 0;
				sortVec( _chiStat_female, _chiStat_female.size());
				if ((_chiStat_female.size())%2!=0) _median = _chiStat_female[((_chiStat_female.size())/2)-1];
				else _median = (_chiStat_female[(_chiStat_female.size()-1)/2 -1] + _chiStat_female[(_chiStat_female.size()-1)/2])/2;
				_g.outputLambda_female = _median/0.4549364;		//median of chi-sq from R ... qchisq(0.5, df= 1)
				if (_g.outputLambda_female>1)LOG << "Output p-value is corrected with lambda=" << _g.outputLambda_female << endl;
				else LOG << "Output p-value is not corrected. Lambda=" << _g.outputLambda_female << endl;
				if (_g.outputLambda_female>1)cout <<"------------------------------------\nOutput p-value is corrected with lambda=" << _g.outputLambda_female << endl;
				else cout <<"------------------------------------\nOutput p-value is not corrected. Lambda=" << _g.outputLambda_female << endl;


			}

	}
	if (_g._markerMap != "N")
		{OUT <<"chromosome\tposition\t";}
	if (_g._binaryTrait)
		{OUT << "rs_number\treference_allele\tother_allele\teaf\tOR\tOR_se\tOR_95L\tOR_95U\tz\tp-value\t_-log10_p-value\tq_statistic\tq_p-value\ti2\tn_studies\tn_samples\teffects";}
	else
		{OUT << "rs_number\treference_allele\tother_allele\teaf\tbeta\tse\tbeta_95L\tbeta_95U\tz\tp-value\t_-log10_p-value\tq_statistic\tq_p-value\ti2\tn_studies\tn_samples\teffects";}

	if (_g._genderHet)
	{
		if (_g._binaryTrait)
		{OUT << "\tmale_eaf\tmale_OR\tmale_OR_se\tmale_OR_95L\tmale_OR_95U\tmale_z\tmale_p-value\tmale_n_studies\tmale_n_samples\tfemale_eaf\tfemale_OR\tfemale_OR_se\tfemale_OR_95L\tfemale_OR_95U\tfemale_z\tfemale_p-value\tfemale_n_studies\tfemale_n_samples\tgender_differentiated_p-value\tgender_heterogeneity_p-value";}
		else
		{OUT << "\tmale_eaf\tmale_beta\tmale_se\tmale_beta_95L\tmale_beta_95U\tmale_z\tmale_p-value\tmale_n_studies\tmale_n_samples\tfemale_eaf\tfemale_beta\tfemale_se\tfemale_beta_95L\tfemale_beta_95U\tfemale_z\tfemale_p-value\tfemale_n_studies\tfemale_n_samples\tgender_differentiated_p-value\tgender_heterogeneity_p-value";}
	}
		OUT << endl;

	
	for (int i = 0; i < _g._markerCount-1; i++)
	{
		if (_g._randomEffect)markerList[i].doRandom(_g);
		OUT << markerList[i].printMarker( _g) << endl;
	}
	
	OUT.close();
    LOG << "Analysis finished." << endl;

	ERR.close();
	LOG.close();
	cout << "------------------------------------\nGWAMA program finished current job successfully!" << endl;
	cout << "Please check " << _g._outputRoot << ".log.out for full information." << endl;
	cout << "Analysis finished." << endl;
	return 0;
}





