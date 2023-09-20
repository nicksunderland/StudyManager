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


#include "marker.h"
#include <iostream>
#include <map>
#include <vector>
#include <stdlib.h>
#include <fstream>
#include <math.h>
#include "statistics.h"
#include "study.h"
#include "chisquaredistr.h"
#include "global.h"
#include "tools.h"

using namespace std;


double
abs_d(double x)
{
if (x < 0) return x * -1;
return x;

}

marker::marker(std::string name, int popCount)
{
		_name = name;
		for (int i = 0; i < popCount; i++)
		{
			_direction += "?";
		}
		_effectAllele="N";
		_otherAllele="N";
		_studyCount = 0;
		male_studyCount = 0;
		female_studyCount = 0;
		_position = -9;
		_chromosome = -9;
		_N=0;
		_Nok=true;
		_EAF=0;
		_B=0;
		_W=0;
		_W2=0;
		_Q=0;

		_Nmale=0;
		_EAFmale=0;
		_Bmale=0;
		_Wmale=0;
		_W2male=0;
		_Qmale=0;
		_Brmale=0;
		_Wrmale=0;
		_Qrmale=0;

		_Nfemale=0;
		_EAFfemale=0;
		_Bfemale=0;
		_Wfemale=0;
		_W2female=0;
		_Qfemale=0;
		_Brfemale=0;
		_Wrfemale=0;
		_Qrfemale=0;

}

marker::~marker(void)
{
}

bool 
marker::addMap(int chromosome, int position)
{
	_chromosome=chromosome;
	_position=position;
	return true;
}


int 
marker::addCohort(int fileNr, global & g, double directLambda, double imputedLambda, 
		string myMarker,string myEffectAllele, string myOtherAllele, double myEffectAlleleFreq,
		bool myStrand, double myBeta, double mySE, bool myImputed, int myN, ofstream & ERR, 
		ofstream & LOG, problem & markerProblems)
{
	//READING DATA
	char sb [11];
	double _se = 0;
	if (g._genomicControl)
	{
		if (imputedLambda<1)imputedLambda = 1;
		if (directLambda<1)directLambda = 1;
		if (myImputed) _se = mySE * sqrt(imputedLambda);
		else _se = mySE * sqrt(directLambda);
	}
	else _se = mySE;

	if (_effectAllele.compare("N")==0)_effectAllele = myEffectAllele;
	if (_otherAllele.compare("N")==0)_otherAllele = myOtherAllele;

	double _beta = 0;
	double _eaf = 0;
	if (_effectAllele == myEffectAllele && _otherAllele == myOtherAllele)
	{
		_beta = myBeta;
		if (myEffectAlleleFreq>0) _eaf = myEffectAlleleFreq;
	}
	else if (_effectAllele == myOtherAllele && _otherAllele == myEffectAllele)//we need to align according to other allele
	{
		_beta = myBeta * -1;
		if (myEffectAlleleFreq>0)  _eaf = 1 - myEffectAlleleFreq;
	}
	else	//something is wrong with alleles - i'll switch strand and check if that helps
	{
		if (_effectAllele == flip(myEffectAllele) && _otherAllele == flip(myOtherAllele))
		{
			_beta = myBeta;
			if (myEffectAlleleFreq>0)  _eaf = myEffectAlleleFreq;
			sprintf(sb, "W%.9d" , g._errNR); g._errNR++;
			LOG << "Warning: Marker " << myMarker << " has wrong strand. (" << sb << ")" <<endl;
			ERR << sb << " " << g.studyList[fileNr] << " warning: Marker " << _name << " has wrong strand. Strand flipped." <<endl;
			markerProblems.problemStrand++;
		}
		else if (_effectAllele == flip(myOtherAllele) && _otherAllele == flip(myEffectAllele))
		{
			_beta = myBeta * -1;
			if (myEffectAlleleFreq>0) _eaf = 1 - myEffectAlleleFreq;
			sprintf(sb, "W%.9d" , g._errNR); g._errNR++;
			LOG << "Warning: Marker " << myMarker << " has wrong strand. (" << sb << ")" <<endl;
			ERR << sb << " " << g.studyList[fileNr] << " warning: Marker " << _name << " has wrong strand. Strand flipped." <<endl;
			markerProblems.problemStrand++;
		}
		else //cannot resolve the problem - marker dropped 
		{
			sprintf(sb, "E%.9d" , g._errNR); g._errNR++;
			LOG << "Error: Marker " << myMarker << " has wrong alleles. (" << sb << ")" <<endl;
			ERR << sb << " " << g.studyList[fileNr] << " error: Marker " << _name << " has wrong alleles. Marker dropped." <<endl;
			markerProblems.wrongAlleles++;
			return 1;
		}
	}
	markerProblems.markersOK++;
	//ADD STATS
	double _wstat = 1/ (_se*_se);
	_B+=_beta*_wstat;
	_W+=_wstat;
	_W2+=_wstat*_wstat;
	_Q+=(_beta*_beta)*_wstat;


	if (g._randomEffect)
	{
		all_beta.push_back(_beta);
		all_se.push_back(_se);
	}

	if (g._genderHet)
	{
		if (g.studyGenders[fileNr]=="M")
		{
			_Bmale+=_beta*_wstat;		
			_Wmale+=_wstat;
			_W2male+=_wstat*_wstat;
			_Qmale+=(_beta*_beta)*_wstat;
			if (g._randomEffect)
			{
				male_beta.push_back(_beta);
				male_se.push_back(_se);
			}

		}
		else if
		(g.studyGenders[fileNr]=="F")
		{
			_Bfemale+=_beta*_wstat;		
			_Wfemale+=_wstat;
			_W2female+=_wstat*_wstat;
			_Qfemale+=(_beta*_beta)*_wstat;
			if (g._randomEffect)
			{
				female_beta.push_back(_beta);
				female_se.push_back(_se);
			}

		}
	}



	//ADD N, POPCOUNT & WEIGHTED EAF
	_studyCount++;

			
        if (abs_d(_EAF - _eaf)>0.3 && _EAF >0 && _eaf >0) //large diversity of allele frequencies - printing warning to log file
        {
                sprintf(sb, "W%.9d" , g._errNR); g._errNR++;
                LOG << "Warning: Marker " << _name << " has large effect allele frequency discrepancy. (" << sb << ")" <<endl;

                ERR << sb << " Warning: Marker " << myMarker << " has large effect allele frequency discrepancy." <<endl;
                ERR << sb << " Weighted EAF: " << _EAF << ", but file " <<  g.studyList[fileNr] << " has EAF: " << _eaf <<endl;
		markerProblems.problemStrand++;
        }




	if (myN>0 && _eaf>0)
	{
		if (_N>=0) _EAF = ((_EAF*_N)+(_eaf*myN))/(_N+myN);
		if (_N>=0)_N+=myN;
		if (g._genderHet)
		{

			if (g.studyGenders[fileNr]=="M")
			{
				male_studyCount++;
				if (_Nmale>=0) _EAFmale = ((_EAFmale*_Nmale)+(_eaf*myN))/(_Nmale+myN);
				if (_Nmale>=0)_Nmale+=myN;
			}
			else if
			(g.studyGenders[fileNr]=="F")
			{
				female_studyCount++;
				if (_Nfemale>=0) _EAFfemale = ((_EAFfemale*_Nfemale)+(_eaf*myN))/(_Nfemale+myN);
				if (_Nfemale>=0)_Nfemale+=myN;
			}
		}
	}
	else if (myN>0)
	{
		if (_N>=0) _N+=myN;
		if (g._genderHet)
		{

			if (g.studyGenders[fileNr]=="M")
			{
				male_studyCount++;
				if (_Nmale>=0) _Nmale+=myN;
			}
			else if
			(g.studyGenders[fileNr]=="F")
			{
				female_studyCount++;
				if (_Nfemale>=0) _Nfemale+=myN;
			}
		}
	}
	else	//stop using N as one of populations didnt had number
	{
		if (_eaf>0)
		{
			if (_EAF>0)_EAF = (_EAF + _eaf)/2;
			else _EAF = _eaf;
		}
		if (g._genderHet)
		{

			if (g.studyGenders[fileNr]=="M")
			{
				male_studyCount++;
			}
			else if
				(g.studyGenders[fileNr]=="F")
			{
				female_studyCount++;
			}
		}
		_Nok=false;
	}

	//direction
	double x = abs_d(_beta * _beta * (1/(_se * _se)));
	if (bdm_chi2(x, 1)<=g._thresholdPValDir)
	{
		if (_beta>0)_direction[fileNr]='+';
		else if (_beta<0)_direction[fileNr]='-';
        else _direction[fileNr]='0';
	}

	return 0;
}


std::string marker::printMarker(global & g)
{
	char sb [1024];		//temporary char buffer
	std::string x="";	//starting to collect all information into this string
	double sbd;
	bool useRandom = g._randomEffect;
	bool map = (g._markerMap!="N");
	bool useMetaLambda = g._genomicControlOutput;
	bool _bt = g._binaryTrait;
	double outputLambda = g.outputLambda;
	double outputLambda_male = g.outputLambda_male;
	double outputLambda_female = g.outputLambda_female;

if (map)
{
	sprintf(sb, "%d\t" , _chromosome); 
	x.append(sb);
	sprintf(sb, "%d\t" , _position); 
	x.append(sb);
}


	x.append(_name);	// printing SNP name
	
	x.append("\t"); //printing reference allele
	x+=_effectAllele;

	x.append("\t"); //printing other allele
	x+=_otherAllele;

	x.append("\t"); // EAF
	if (_EAF>0)sprintf(sb, "%.6f" , _EAF); 
	else sprintf(sb, "%d" , -9); 
	x.append(sb);

	if (_bt)
	{
		x.append("\t");//printing OR
		sbd = OR(useRandom);  
		sprintf(sb, "%.6f" , sbd);
		x.append(sb);

		x.append("\t");//printing se
		sbd = (OR(useRandom)-Lower95OR(useRandom))/1.96;  
		sprintf(sb, "%.6f" , sbd);
		x.append(sb);


		x.append("\t");//printing L95 OR
		sbd = Lower95OR(useRandom);  
		sprintf(sb, "%.6f" , sbd);
		x.append(sb);

		x.append("\t");//printing U95 OR
		sbd = Upper95OR(useRandom);  
		sprintf(sb, "%.6f" , sbd);
		x.append(sb);	
	}
	else 
	{
		x.append("\t");//printing average beta
		sbd = sumBeta(useRandom);  
		sprintf(sb, "%.6f" , sbd);
		x.append(sb);

		x.append("\t");//printing se
		sbd = sumSE(useRandom);  
		sprintf(sb, "%.6f" , sbd);
		x.append(sb);


		x.append("\t");//printing L95 beta
		sbd = Lower95Beta(useRandom);  
		sprintf(sb, "%.6f" , sbd);
		x.append(sb);

		x.append("\t");//printing U95 beta
		sbd = Upper95Beta(useRandom);  
		sprintf(sb, "%.6f" , sbd);
		x.append(sb);
	}

	x.append("\t");//printing z
	sbd = z(useRandom);  
	sprintf(sb, "%.6f" , sbd);
	x.append(sb);

	x.append("\t");//p-value
	if (useMetaLambda) sbd = pValue(outputLambda, useRandom);
	else sbd = pValue(1, useRandom);
	if (sbd>=0.001) {sprintf(sb, "%.6f" , sbd);}
	else {sprintf(sb, "%.2E" , sbd);}
	x.append(sb);

	x.append("\t");//-10 log p-value
	sprintf(sb, "%.6f" , -log10(sbd));
	x.append(sb);

	x.append("\t");//q-statistic
	sbd = qStatistic();  
	sprintf(sb, "%.6f" , sbd);
	x.append(sb);

	x.append("\t");//qPValue
	sbd = qPValue();  
	if (sbd>=0.001) {sprintf(sb, "%.6f" , sbd);}
	else if (sbd<0.001){sprintf(sb, "%.2E" , sbd);}
	else {sprintf(sb, "1");}
	x.append(sb);

	x.append("\t");//i2
	sbd = i2();  
	sprintf(sb, "%.6f" , sbd);
	x.append(sb);

	x.append("\t"); // printing study nr
	sprintf(sb, "%d" , _studyCount); 
	x.append(sb);

	x.append("\t"); // printing sample nr
	if (_Nok)sprintf(sb, "%d" , _N); 
	else sprintf(sb, "%d" , -9); 
	x.append(sb);

	x.append("\t");				//printing effect directions
	x.append(_direction);

/////////////
	if (g._genderHet)
	{

/////////////MALE
		if (male_studyCount>0)
		{
			x.append("\t"); // EAF
			if (_EAFmale>0)sprintf(sb, "%.6f" , _EAFmale); 
			else sprintf(sb, "%d" , -9); 
			x.append(sb);

			if (_bt)
			{
				x.append("\t");//printing OR
				sbd = male_OR(useRandom);  
				sprintf(sb, "%.6f" , sbd);
				x.append(sb);

				x.append("\t");//printing se
				sbd = (male_OR(useRandom)-male_Lower95OR(useRandom))/1.96;  
				sprintf(sb, "%.6f" , sbd);
				x.append(sb);


				x.append("\t");//printing L95 OR
				sbd = male_Lower95OR(useRandom);  
				sprintf(sb, "%.6f" , sbd);
				x.append(sb);

				x.append("\t");//printing U95 OR
				sbd = male_Upper95OR(useRandom);  
				sprintf(sb, "%.6f" , sbd);
				x.append(sb);	
			}
			else 
			{
				x.append("\t");//printing average beta
				sbd = male_sumBeta(useRandom);  
				sprintf(sb, "%.6f" , sbd);
				x.append(sb);

				x.append("\t");//printing se
				sbd = male_sumSE(useRandom);  
				sprintf(sb, "%.6f" , sbd);
				x.append(sb);


				x.append("\t");//printing L95 beta
				sbd = male_Lower95Beta(useRandom);  
				sprintf(sb, "%.6f" , sbd);
				x.append(sb);

				x.append("\t");//printing U95 beta
				sbd = male_Upper95Beta(useRandom);  
				sprintf(sb, "%.6f" , sbd);
				x.append(sb);
			}

			x.append("\t");//printing z
			sbd = male_z(useRandom);  
			sprintf(sb, "%.6f" , sbd);
			x.append(sb);

			x.append("\t");//p-value
			if (useMetaLambda) sbd = male_pValue(g, useRandom);
			else sbd = male_pValue(g, useRandom);
			if (sbd>=0.001) {sprintf(sb, "%.6f" , sbd);}
			else {sprintf(sb, "%.2E" , sbd);}
			x.append(sb);

			x.append("\t"); // printing study nr
			sprintf(sb, "%d" , male_studyCount); 
			x.append(sb);

			x.append("\t"); // printing sample nr
			if (_Nok)sprintf(sb, "%d" , _Nmale); 
			else sprintf(sb, "%d" , -9); 
			x.append(sb);
		}
		else
		{	
			x.append("\t-9\t-9\t-9\t-9\t-9\t-9\t-9\t-9\t-9");
		}

///////////FEMALE
		if (female_studyCount>0)
		{
			x.append("\t"); // EAF
			if (_EAFfemale>0)sprintf(sb, "%.6f" , _EAFfemale); 
			else sprintf(sb, "%d" , -9); 
			x.append(sb);

			if (_bt)
			{
				x.append("\t");//printing OR
				sbd = female_OR(useRandom);  
				sprintf(sb, "%.6f" , sbd);
				x.append(sb);

				x.append("\t");//printing se
				sbd = (female_OR(useRandom) - female_Lower95OR(useRandom))/1.96;
				sprintf(sb, "%.6f" , sbd);
				x.append(sb);


				x.append("\t");//printing L95 OR
				sbd = female_Lower95OR(useRandom);  
				sprintf(sb, "%.6f" , sbd);
				x.append(sb);

				x.append("\t");//printing U95 OR
				sbd = female_Upper95OR(useRandom);  
				sprintf(sb, "%.6f" , sbd);
				x.append(sb);	
			}
			else 
			{
				x.append("\t");//printing average beta
				sbd = female_sumBeta(useRandom);  
				sprintf(sb, "%.6f" , sbd);
				x.append(sb);

				x.append("\t");//printing se
				sbd = female_sumSE(useRandom);  
				sprintf(sb, "%.6f" , sbd);
				x.append(sb);


				x.append("\t");//printing L95 beta
				sbd = female_Lower95Beta(useRandom);  
				sprintf(sb, "%.6f" , sbd);
				x.append(sb);

				x.append("\t");//printing U95 beta
				sbd = female_Upper95Beta(useRandom);  
				sprintf(sb, "%.6f" , sbd);
				x.append(sb);
			}

			x.append("\t");//printing z
			sbd = female_z(useRandom);  
			sprintf(sb, "%.6f" , sbd);
			x.append(sb);

			x.append("\t");//p-value
			if (useMetaLambda) sbd = female_pValue(g, useRandom);
			else sbd = female_pValue(g, useRandom);
			if (sbd>=0.001) {sprintf(sb, "%.6f" , sbd);}
			else {sprintf(sb, "%.2E" , sbd);}
			x.append(sb);

			x.append("\t"); // printing study nr
			sprintf(sb, "%d" , female_studyCount); 
			x.append(sb);

			x.append("\t"); // printing sample nr
			if (_Nok)sprintf(sb, "%d" , _Nfemale); 
			else sprintf(sb, "%d" , -9); 
			x.append(sb);
		}
		else
		{	
			x.append("\t-9\t-9\t-9\t-9\t-9\t-9\t-9\t-9\t-9");
		}

		if (male_studyCount>0 && female_studyCount>0)
		{
			x.append("\t");//p-value
			if (useMetaLambda) sbd = genderdiff_pValue(g, useRandom);
			else sbd = genderdiff_pValue(g, useRandom);
			if (sbd>=0.001) {sprintf(sb, "%.6f" , sbd);}
			else {sprintf(sb, "%.2E" , sbd);}
			x.append(sb);

			x.append("\t");//p-value
			if (useMetaLambda) sbd = genderhet_pValue(g, useRandom);
			else sbd = genderhet_pValue(g, useRandom);
			if (sbd>=0.001) {sprintf(sb, "%.6f" , sbd);}
			else {sprintf(sb, "%.2E" , sbd);}
			x.append(sb);

		}
		else
		{	
			x.append("\t-9\t-9");
		}
	}
/////////////////////






	return x;

}

bool
marker::doRandom(global & g)
{
		double statQ = _Q - _W*((_B / _W)*(_B / _W));
		double tau2 =0;
		if ((_W-(_W2/_W))!=0)tau2=(statQ -(_studyCount-1))/(_W-(_W2/_W));
		if (tau2<0)tau2=0;

		_Wr = 0;
		_Br = 0;
		_Qr = 0;
		for (int i = 0; i < all_se.size(); i++)
		{
			double se2;
			se2 = all_se[i]*all_se[i];
			_Wr += (1 / (se2 + tau2));
			_Br += all_beta[i] * (1/(se2 + tau2));
			_Qr += all_beta[i] * all_beta[i] * (1/se2);
		}

		if (g._genderHet)
		{
/////////////////
				statQ = _Qmale - _Wmale*((_Bmale / _Wmale)*(_Bmale / _Wmale));
				tau2 =0;
				if ((_Wmale-(_W2male/_Wmale))!=0)tau2=(statQ -(male_studyCount-1))/(_Wmale-(_W2male/_Wmale));
				if (tau2<0)tau2=0;

				_Wrmale = 0;
				_Brmale = 0;
				_Qrmale = 0;
				for (int i = 0; i < male_se.size(); i++)
				{
					double se2;
					se2 = male_se[i]*male_se[i];
					_Wrmale += (1 / (se2 + tau2));
					_Brmale += male_beta[i] * (1/(se2 + tau2));
					_Qrmale += male_beta[i] * male_beta[i] * (1/se2);
				}



//////////////////
				statQ = _Qfemale - _Wfemale*((_Bfemale / _Wfemale)*(_Bfemale / _Wfemale));
				tau2 =0;
				if ((_Wfemale-(_W2female/_Wfemale))!=0)tau2=(statQ -(female_studyCount-1))/(_Wfemale-(_W2female/_Wfemale));
				if (tau2<0)tau2=0;

				_Wrfemale = 0;
				_Brfemale = 0;
				_Qrfemale = 0;
				for (int i = 0; i < female_se.size(); i++)
				{
					double se2;
					se2 = female_se[i]*female_se[i];
					_Wrfemale += (1 / (se2 + tau2));
					_Brfemale += female_beta[i] * (1/(se2 + tau2));
					_Qrfemale += female_beta[i] * female_beta[i] * (1/se2);
				}


/////////////////



		}
return true;

}



double marker::sumBeta(bool useRandom)
{
	if (useRandom)return _Br / _Wr;
	return _B / _W;
}
double marker::sumSE(bool useRandom)
{
	if (useRandom)return (1/(sqrt(_Wr)));
	return (1/(sqrt(_W)));
}

double marker::Lower95Beta(bool useRandom)
{
	if (useRandom)return (_Br / _Wr) - (1.96/(sqrt(_Wr)));
	return (_B / _W) - (1.96/(sqrt(_W)));
}
double marker::Upper95Beta(bool useRandom)
{
	if (useRandom)return (_Br / _Wr) + (1.96/(sqrt(_Wr)));
	return (_B / _W) + (1.96/(sqrt(_W)));
}
double marker::OR(bool useRandom)
{
	if (useRandom)return exp(_Br / _Wr);
	return exp(_B / _W);
}
double marker::Upper95OR(bool useRandom)
{
	if (useRandom)return exp((_Br / _Wr) + (1.96/(sqrt(_Wr))));
	return exp((_B / _W) + (1.96/(sqrt(_W))));
}
double marker::Lower95OR(bool useRandom)
{
	if (useRandom)	return exp((_Br / _Wr) - (1.96/(sqrt(_Wr)))) ;
	return exp((_B / _W) - (1.96/(sqrt(_W)))) ;
}
double marker::z(bool useRandom)
{
	if (useRandom)	return (_Br / _Wr) * sqrt(_Wr);
	return (_B / _W) * sqrt(_W);
}
double marker::pValue(double outputLambda, bool useRandom)
{
	double x;
	double beta;
	double se;
	if (outputLambda>1)
	{	
		if (useRandom)
		{
			 beta = _Br / _Wr;
			se = (1/(sqrt(_Wr)));
		}
		else
		{
			beta = _B / _W;
			 se = (1/(sqrt(_W)));
		}
		se = se * sqrt(outputLambda);
		x = (beta/se)*(beta/se);
	}
	else 
	{	
		if (useRandom)x = ((_Br / _Wr) * sqrt(_Wr))*((_Br / _Wr) * sqrt(_Wr));
		else x = ((_B / _W) * sqrt(_W))*((_B / _W) * sqrt(_W));
	}
	double z = abs_d(x);
	return bdm_chi2(z, 1);
	return 1;
}

double 
marker::genderdiff_pValue(global & g, bool useRandom)
{
	double x;
	double beta;
	double se;
	if (g.outputLambda_male>1)
	{	
		if (useRandom)
		{
			 beta = _Brmale / _Wrmale;
			se = (1/(sqrt(_Wrmale)));
		}
		else
		{
			beta = _Bmale / _Wmale;
			 se = (1/(sqrt(_Wmale)));
		}
		se = se * sqrt(g.outputLambda_male);
		x = (beta/se)*(beta/se);
	}
	else 
	{	
		if (useRandom)x = ((_Brmale / _Wrmale) * sqrt(_Wrmale))*((_Brmale / _Wrmale) * sqrt(_Wrmale));
		else x = ((_Bmale / _Wmale) * sqrt(_Wmale))*((_Bmale / _Wmale) * sqrt(_Wmale));
	}
	double z = abs_d(x);

//FEMALE
	if (g.outputLambda_female>1)
	{	
		if (useRandom)
		{
			 beta = _Brfemale / _Wrfemale;
			se = (1/(sqrt(_Wrfemale)));
		}
		else
		{
			beta = _Bfemale / _Wfemale;
			 se = (1/(sqrt(_Wfemale)));
		}
		se = se * sqrt(g.outputLambda_female);
		x = (beta/se)*(beta/se);
	}
	else 
	{	
		if (useRandom)x = ((_Brfemale / _Wrfemale) * sqrt(_Wrfemale))*((_Brfemale / _Wrfemale) * sqrt(_Wrfemale));
		else x = ((_Bfemale / _Wfemale) * sqrt(_Wfemale))*((_Bfemale / _Wfemale) * sqrt(_Wfemale));
	}

	z += abs_d(x);
	return bdm_chi2(z, 2);

	return 1;
}
double 
marker::genderhet_pValue(global & g, bool useRandom)
{
	double x;
	double beta;
	double se;
	if (g.outputLambda>1)
	{	
		if (useRandom)
		{
			 beta = _Br / _Wr;
			se = (1/(sqrt(_Wr)));
		}
		else
		{
			beta = _B / _W;
			 se = (1/(sqrt(_W)));
		}
		se = se * sqrt(g.outputLambda);
		x = (beta/se)*(beta/se);
	}
	else 
	{	
		if (useRandom)x = ((_Br / _Wr) * sqrt(_Wr))*((_Br / _Wr) * sqrt(_Wr));
		else x = ((_B / _W) * sqrt(_W))*((_B / _W) * sqrt(_W));
	}
	double z = abs_d(x);


	if (g.outputLambda_male>1)
	{	
		if (useRandom)
		{
			 beta = _Brmale / _Wrmale;
			se = (1/(sqrt(_Wrmale)));
		}
		else
		{
			beta = _Bmale / _Wmale;
			 se = (1/(sqrt(_Wmale)));
		}
		se = se * sqrt(g.outputLambda_male);
		x = (beta/se)*(beta/se);
	}
	else 
	{	
		if (useRandom)x = ((_Brmale / _Wrmale) * sqrt(_Wrmale))*((_Brmale / _Wrmale) * sqrt(_Wrmale));
		else x = ((_Bmale / _Wmale) * sqrt(_Wmale))*((_Bmale / _Wmale) * sqrt(_Wmale));
	}
	double z2 = abs_d(x);

//FEMALE
	if (g.outputLambda_female>1)
	{	
		if (useRandom)
		{
			 beta = _Brfemale / _Wrfemale;
			se = (1/(sqrt(_Wrfemale)));
		}
		else
		{
			beta = _Bfemale / _Wfemale;
			 se = (1/(sqrt(_Wfemale)));
		}
		se = se * sqrt(g.outputLambda_female);
		x = (beta/se)*(beta/se);
	}
	else 
	{	
		if (useRandom)x = ((_Brfemale / _Wrfemale) * sqrt(_Wrfemale))*((_Brfemale / _Wrfemale) * sqrt(_Wrfemale));
		else x = ((_Bfemale / _Wfemale) * sqrt(_Wfemale))*((_Bfemale / _Wfemale) * sqrt(_Wfemale));
	}

	z2 += abs_d(x);

	return bdm_chi2(abs_d(z2-z), 1);

	return 1;
}



double marker::qStatistic()
{
	return _Q - _W*((_B / _W)*(_B / _W));
}
double marker::qPValue()
{
	double x = _Q - _W*((_B / _W)*(_B / _W));
	if (x>0 && _studyCount > 1) return 1-chisquaredistribution(_studyCount-1, x);	
	return 1;
}
double marker::i2()
{
	double h2 = (_Q-_W*((_B / _W)*(_B / _W))) /(_studyCount-1);
	if (h2<1) return 0;
	return (h2-1)/h2;
}


////////////////////////////
//MALE
double marker::male_sumBeta(bool useRandom)
{
	if (useRandom)return _Brmale / _Wrmale;
	return _Bmale / _Wmale;
}
double marker::male_sumSE(bool useRandom)
{
	if (useRandom)return (1/(sqrt(_Wrmale)));
	return (1/(sqrt(_Wmale)));
}

double marker::male_Lower95Beta(bool useRandom)
{
	if (useRandom)return (_Brmale / _Wrmale) - (1.96/(sqrt(_Wrmale)));
	return (_Bmale / _Wmale) - (1.96/(sqrt(_Wmale)));
}
double marker::male_Upper95Beta(bool useRandom)
{
	if (useRandom)return (_Brmale / _Wrmale) + (1.96/(sqrt(_Wrmale)));
	return (_Bmale / _Wmale) + (1.96/(sqrt(_Wmale)));
}
double marker::male_OR(bool useRandom)
{
	if (useRandom)return exp(_Brmale / _Wrmale);
	return exp(_Bmale / _Wmale);
}
double marker::male_Upper95OR(bool useRandom)
{
	if (useRandom)return exp((_Brmale / _Wrmale) + (1.96/(sqrt(_Wrmale))));
	return exp((_Bmale / _Wmale) + (1.96/(sqrt(_Wmale))));
}
double marker::male_Lower95OR(bool useRandom)
{
	if (useRandom)	return exp((_Brmale / _Wrmale) - (1.96/(sqrt(_Wrmale)))) ;
	return exp((_Bmale / _Wmale) - (1.96/(sqrt(_Wmale)))) ;
}
double marker::male_z(bool useRandom)
{
	if (useRandom)	return (_Brmale / _Wrmale) * sqrt(_Wrmale);
	return (_Bmale / _Wmale) * sqrt(_Wmale);
}
double marker::male_pValue(global & g, bool useRandom)
{
	double x;
	double beta;
	double se;
	double outputLambda = g.outputLambda_male;
	if (outputLambda>1)
	{	
		if (useRandom)
		{
			 beta = _Brmale / _Wrmale;
			se = (1/(sqrt(_Wrmale)));
		}
		else
		{
			beta = _Bmale / _Wmale;
			 se = (1/(sqrt(_Wmale)));
		}
		se = se * sqrt(outputLambda);
		x = (beta/se)*(beta/se);
	}
	else 
	{	
		if (useRandom)x = ((_Brmale / _Wrmale) * sqrt(_Wrmale))*((_Brmale / _Wrmale) * sqrt(_Wrmale));
		else x = ((_Bmale / _Wmale) * sqrt(_Wmale))*((_Bmale / _Wmale) * sqrt(_Wmale));
	}
	double z = abs_d(x);
	return bdm_chi2(z, 1);
	return 1;
}

////////////////////////////////
//FEMALE

double marker::female_sumBeta(bool useRandom)
{
	if (useRandom)return _Brfemale / _Wrfemale;
	return _Bfemale / _Wfemale;
}
double marker::female_sumSE(bool useRandom)
{
	if (useRandom)return (1/(sqrt(_Wrfemale)));
	return (1/(sqrt(_Wfemale)));
}

double marker::female_Lower95Beta(bool useRandom)
{
	if (useRandom)return (_Brfemale / _Wrfemale) - (1.96/(sqrt(_Wrfemale)));
	return (_Bfemale / _Wfemale) - (1.96/(sqrt(_Wfemale)));
}
double marker::female_Upper95Beta(bool useRandom)
{
	if (useRandom)return (_Brfemale / _Wrfemale) + (1.96/(sqrt(_Wrfemale)));
	return (_Bfemale / _Wfemale) + (1.96/(sqrt(_Wfemale)));
}
double marker::female_OR(bool useRandom)
{
	if (useRandom)return exp(_Brfemale / _Wrfemale);
	return exp(_Bfemale / _Wfemale);
}
double marker::female_Upper95OR(bool useRandom)
{
	if (useRandom)return exp((_Brfemale / _Wrfemale) + (1.96/(sqrt(_Wrfemale))));
	return exp((_Bfemale / _Wfemale) + (1.96/(sqrt(_Wfemale))));
}
double marker::female_Lower95OR(bool useRandom)
{
	if (useRandom)	return exp((_Brfemale / _Wrfemale) - (1.96/(sqrt(_Wrfemale)))) ;
	return exp((_Bfemale / _Wfemale) - (1.96/(sqrt(_Wfemale)))) ;
}
double marker::female_z(bool useRandom)
{
	if (useRandom)	return (_Brfemale / _Wrfemale) * sqrt(_Wrfemale);
	return (_Bfemale / _Wfemale) * sqrt(_Wfemale);
}
double marker::female_pValue(global & g, bool useRandom)
{
	double x;
	double beta;
	double se;
	double outputLambda = g.outputLambda_female;
	if (outputLambda>1)
	{	
		if (useRandom)
		{
			 beta = _Brfemale / _Wrfemale;
			se = (1/(sqrt(_Wrfemale)));
		}
		else
		{
			beta = _Bfemale / _Wfemale;
			 se = (1/(sqrt(_Wfemale)));
		}
		se = se * sqrt(outputLambda);
		x = (beta/se)*(beta/se);
	}
	else 
	{	
		if (useRandom)x = ((_Brfemale / _Wrfemale) * sqrt(_Wrfemale))*((_Brfemale / _Wrfemale) * sqrt(_Wrfemale));
		else x = ((_Bfemale / _Wfemale) * sqrt(_Wfemale))*((_Bfemale / _Wfemale) * sqrt(_Wfemale));
	}
	double z = abs_d(x);
	return bdm_chi2(z, 1);
	return 1;
}
