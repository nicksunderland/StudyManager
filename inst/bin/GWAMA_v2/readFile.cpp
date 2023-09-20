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

#include <iostream>
#include <map>
#include <vector>
#include <stdlib.h>
#include <fstream>
#include "marker.h"
#include "study.h"
#include <math.h>
#include <zlib.h>
#include <cctype> // std::toupper
#include "global.h"
#include "tools.h"
#include "cohort.h"
#include "problem.h"

#define LENS 1000000
using namespace std;

bool isFiltered(double value, string equation, double threshold)
{
    if (equation=="==" && value==threshold) return true;
    if (equation=="!=" && value!=threshold) return true;
    if (equation==">=" && value>=threshold) return true;
    if (equation=="<=" && value<=threshold) return true;
    if (equation==">" && value>threshold) return true;
    if (equation=="<" && value<threshold) return true;
    return false;
}

bool 
readFilelist(global & _g, ofstream & ERR,  ofstream & LOG)
{
   	ifstream IF (_g._fileList.c_str());
	if (IF.is_open())
	{
		LOG << "Reading input file list:" << endl;
		while (! IF.eof() )
    	{
            	string line;
            	vector<string> tokens;
            	getline (IF,line);
            	int n = Tokenize(string(line), tokens, " ");
            	if (n)
            	{
					LOG << "\t" << string(tokens[0]) << endl;
					_g.studyList.push_back(string(tokens[0]));
					if (_g._genderHet && n>1)
					{
						if (uc(string(tokens[1]))=="M" 
							|| uc(string(tokens[1]))=="F") 
							_g.studyGenders.push_back(uc(string(tokens[1])));
						else
						{
							cerr <<"Cannot read gender status for file : " << string(tokens[0]) << ". Exit program!" << endl;
							exit (1);
						}
					}
					else if
						(_g._genderHet)
					{
							cerr <<"Cannot read gender status for file : " << string(tokens[0]) << ". Exit program!" << endl;
							exit (1);

					}
				}
		}
		LOG << "END-OF-FILE" << endl;
	}
	IF.close();
	_g._studyCount = _g.studyList.size();
	LOG << "Study count: " << _g._studyCount <<  endl;
	cout << "File " <<_g._fileList << " contained "<< _g._studyCount << " studies." << endl;
	return true;
}

bool
getLambda(string name, double & directLambda, double & imputedLambda, 
		  int & directCount, int & imputedCount, int columnBeta, 
		  int columnSE, int columnImputed, int columnCount, vector <int> filterColumns, global & g)
{
	int currentLine =0;
	int _imputedCount = 0;
	int _directCount = 0;
	vector <double> _imputedChiStat;
	vector <double> _directChiStat;
    
     if (name.substr(name.length()-2)=="gz")
    {
        gzFile F =gzopen(name.c_str(),"r");  
       
        
        char *buffer = new char[LENS];
        while(0!=gzgets(F,buffer,LENS))
            {
                string line;
                vector<string> tokens;

                int n = Tokenize(buffer, tokens, " ");
                if (n)
                {
                    if (currentLine>0) // skip headers
                    {
                        bool filtered=0;
                        double beta=0;
                        double se=0;
                        if (g._binaryTrait)
                        {
                            double oratio = atof(tokens[columnBeta].c_str());
                            double oratio_95l = atof(tokens[columnSE].c_str());
                            if (oratio> 0 && oratio > oratio_95l && oratio_95l>0)
                            {
                                beta=log(oratio);
                                se = ((log(oratio)-log(oratio_95l))/1.96);
                            }
                        }
                        else
                        {
                            beta = atof(tokens[columnBeta].c_str());
                            se = atof(tokens[columnSE].c_str());
                        }
                        for (int j; j<filterColumns.size();j++)
                        {
                            if (filterColumns[j] != -1)
                            {
                                double valueX = atof(tokens[filterColumns[j]].c_str());
                                if (isFiltered(valueX, g.filters[j]._equation, g.filters[j]._value))
                                    filtered=1;
                            }
                        }
                        int imputed = 0;
                        if (columnImputed!=-9)
                        {
                            if (string(tokens[columnImputed])=="1" || string(tokens[columnImputed])=="0")
                            { 
                                if (atoi(tokens[columnImputed].c_str())==1)imputed=1;
                            }
                            
                        }
                        if (se>0 && filtered==0)
                        {
                            if (imputed)
                            {
                                _imputedCount++;
                                _imputedChiStat.push_back((beta*beta)/(se*se));
                            }
                            else
                            {
                                _directCount++;
                                _directChiStat.push_back((beta*beta)/(se*se));
                            }
                        }
                    }
                }
                if (int(currentLine/10000)==currentLine/10000.0)	cout << "Finding lambda. Line: " << currentLine << "\r";
                
                currentLine++;
            }
        
        
        if (_imputedCount>0)
        {
            double _medianI = 0;
            sortVec( _imputedChiStat, _imputedCount);
            if (_imputedCount%2!=0) _medianI = _imputedChiStat[((_imputedCount+1)/2)-1];
            else _medianI = (_imputedChiStat[_imputedCount/2 -1] + _imputedChiStat[_imputedCount/2])/2;
            imputedLambda = _medianI/0.4549364;		//median of chi-sq from R ... qchisq(0.5, df= 1)
            
        }
        
        if (_directCount>0)
        {
            double _medianD = 0;
            sortVec( _directChiStat, _directCount);
            if (_directCount%2!=0) _medianD = _directChiStat[((_directCount+1)/2)-1];
            else _medianD = (_directChiStat[_directCount/2 -1] + _directChiStat[_directCount/2])/2;
            directLambda = _medianD/0.4549364;        //median of chi-sq from R ... qchisq(0.5, df= 1)
        }
        directCount = _directCount;
        imputedCount = _imputedCount;
        
    }
 
    else
    {
        ifstream F (name.c_str());
        if (F.is_open())
        {	
            while (! F.eof() )
            {
                string line;
                vector<string> tokens;
                getline (F,line);
                int n = Tokenize(string(line), tokens, " ");
                if (n)
                {
                    if (currentLine>0) // skip headers
                    {
                        double beta=0;
                        double se=0;
                        bool filtered = 0;
                        if (g._binaryTrait)
                        {
                            double oratio = atof(tokens[columnBeta].c_str());
                            double oratio_95l = atof(tokens[columnSE].c_str());
                            if (oratio> 0 && oratio > oratio_95l && oratio_95l>0)
                            {
                            beta=log(oratio);
                            se = ((log(oratio)-log(oratio_95l))/1.96);
                            }
                        }
                        else
                        {
                        beta = atof(tokens[columnBeta].c_str());
                        se = atof(tokens[columnSE].c_str());
                        }
                        int imputed = 0;
                        if (columnImputed!=-9)
                        {
                            if (string(tokens[columnImputed])=="1" || string(tokens[columnImputed])=="0")
                            { 
                                if (atoi(tokens[columnImputed].c_str())==1)imputed=1;
                            }

                        }
                        for (int j; j<filterColumns.size();j++)
                        {
                            if (filterColumns[j] != -1)
                            {
                                double valueX = atof(tokens[filterColumns[j]].c_str());
                                if (isFiltered(valueX, g.filters[j]._equation, g.filters[j]._value))
                                    filtered=1;
                            }
                        }

                        if (se>0 && filtered==0)
                        {
                            if (imputed)
                            {
                                _imputedCount++;
                                _imputedChiStat.push_back((beta*beta)/(se*se));
                            }
                            else
                            {
                                _directCount++;
                                _directChiStat.push_back((beta*beta)/(se*se));
                            }
                        }
                    }
                }
                if (int(currentLine/10000)==currentLine/10000.0)	cout << "Finding lambda. Line: " << currentLine << "\r";

                currentLine++;
            }
        }

        if (_imputedCount>0)
        {
            double _medianI = 0;
            sortVec( _imputedChiStat, _imputedCount);
            if (_imputedCount%2!=0) _medianI = _imputedChiStat[((_imputedCount+1)/2)-1];
            else _medianI = (_imputedChiStat[_imputedCount/2 -1] + _imputedChiStat[_imputedCount/2])/2;
            imputedLambda = _medianI/0.4549364;		//median of chi-sq from R ... qchisq(0.5, df= 1)

        }

        if (_directCount>0)
        {
            double _medianD = 0;
            sortVec( _directChiStat, _directCount);
            if (_directCount%2!=0) _medianD = _directChiStat[((_directCount+1)/2)-1];
            else _medianD = (_directChiStat[_directCount/2 -1] + _directChiStat[_directCount/2])/2;
            directLambda = _medianD/0.4549364;        //median of chi-sq from R ... qchisq(0.5, df= 1)
        }
        directCount = _directCount;
        imputedCount = _imputedCount;
    
    }
    
    

	return true;
}

bool 
readMapFile(string markermap, vector <marker> & markerlist,
             map<string, int> & markerPosition,
				 ofstream & ERR, int & errNR, ofstream & LOG)
{
	char sb [11];
	cout << "Reading map: " << markermap << endl;
	int currentLine = 0;
	ifstream F (markermap.c_str());
	if (F.is_open())
    {
       	while (! F.eof() )
        {
	        string line;
        	vector<string> tokens;
        	getline (F,line);
			int chr = 0;
			int pos = 0;
			string currentmarker = "";
			int n = Tokenize(string(line), tokens, " ");
            if (n>2)
            {
				if (currentLine>0)
				{
					chr = atoi(tokens[0].c_str());
					currentmarker = string(tokens[1]);
					if (n==3) pos = atoi(tokens[2].c_str());
					else pos = atoi(tokens[3].c_str());
					if (markerPosition[currentmarker]!=0)
					{
						markerlist[markerPosition[currentmarker]-1].addMap(chr,pos);
					}
				}
			}
			currentLine ++;
		}
	}
	else {
		cerr << "Error: Cannot find or access " << markermap << ". Positions not added." << endl;
		sprintf(sb, "E%.9d" , errNR); errNR++;
		LOG << "Error: Cannot find or access " << markermap << ". Positions not added! (" << sb << ")" <<endl;
		ERR << sb << " Error: Cannot find or access " << markermap << ". Positions not added!" <<endl;
		ERR << sb << " Error: Make sure that this file is in given path and is readable." <<endl;
		ERR << sb << " Error: Default 'marker.map' file can be downloaded from GWAMA web page." <<endl;
		return false;
	}
	return true;
}


bool
readCohort(int fileNr, global & g,map<string, int> & markerPosition,vector <marker> & markerlist, ofstream & ERR, ofstream & LOG,  ofstream & LAMBDA_FILE)
{
	int currentLine = 0;
	int _countGoodMarkers=0;
	int _countBadMarkers=0;
    
    vector <int> _filterColumnNr;
    vector <unsigned int> _markerFilteredBy;
    for (int i = 0; i < g.filters.size(); i++)
    {
        _filterColumnNr.push_back(-1);
        _markerFilteredBy.push_back(0);
    }
	char sb [11];
	cohort thisCohort;
	problem _markerProblems;

	map<string, int>  _markernames; //collect all OK marker names for checking for duplicates
	int _markerNamesCount=1;

	thisCohort._name=g.studyList[fileNr];
//	cout << "Reading file: " << g.studyList[fileNr] << endl;
	LOG << " [" << fileNr+1 << "]" << " Reading file: " << g.studyList[fileNr]  << endl;
    
     if (g.studyList[fileNr].substr(g.studyList[fileNr].length()-2)=="gz")
    {
        //      	ifstream F (G.inputGenFile.c_str());
        gzFile F =gzopen(g.studyList[fileNr].c_str(),"r");   

        char *buffer = new char[LENS];
        while(0!=gzgets(F,buffer,LENS))
        {
                string line;
                vector<string> tokens;
                int n = Tokenize(buffer, tokens, " ");
                if (n)
                {
                    if (currentLine==0) // get headers
                    {
                        for (int i=0; i<n; i++)
                        {
                            //markername
                            for (unsigned int j=0; j<g._alternativeMarkerNames.size(); j++)
                            {
                                if (uc(string(tokens[i])).compare(g._alternativeMarkerNames[j])==0){thisCohort._columnMarkerName=i;}
                            }
                            //effect allele
                            for (unsigned int j=0; j<g._alternativeEffectAlleles.size(); j++)
                            {
                                if (uc(string(tokens[i])).compare(g._alternativeEffectAlleles[j])==0)thisCohort._columnEffectAllele=i;
                            }
                            //other allele
                            for (unsigned int j=0; j<g._alternativeOtherAlleles.size(); j++)
                            {
                                if (uc(string(tokens[i])).compare(g._alternativeOtherAlleles[j])==0)thisCohort._columnOtherAllele=i;
                            }
                            //effect allele freq
                            for (unsigned int j=0; j<g._alternativeEffectAlleleFreqs.size(); j++)
                            {
                                if (uc(string(tokens[i])).compare(g._alternativeEffectAlleleFreqs[j])==0)thisCohort._columnEffectAlleleFreq=i;
                            }
                            //strand
                            for (unsigned int j=0; j<g._alternativeStrands.size(); j++)
                            {
                                if (uc(string(tokens[i])).compare(g._alternativeStrands[j])==0)thisCohort._columnStrand=i;
                            }
                            //n
                            for (unsigned int j=0; j<g._alternativeNs.size(); j++)
                            {
                                if (uc(string(tokens[i])).compare(g._alternativeNs[j])==0)thisCohort._columnN=i;
                            }
                            //beta
                            for (unsigned int j=0; j<g._alternativeBetas.size(); j++)
                            {
                                if (uc(string(tokens[i])).compare(g._alternativeBetas[j])==0)thisCohort._columnBeta=i;
                            }
                            //se
                            for (unsigned int j=0; j<g._alternativeSEs.size(); j++)
                            {
                                if (uc(string(tokens[i])).compare(g._alternativeSEs[j])==0)thisCohort._columnSE=i;
                            }
                            //or
                            for (unsigned int j=0; j<g._alternativeORs.size(); j++)
                            {
                                if (uc(string(tokens[i])).compare(g._alternativeORs[j])==0)thisCohort._columnOR=i;
                            }
                            //or 95L
                            for (unsigned int j=0; j<g._alternativeOR_95Ls.size(); j++)
                            {
                                if (uc(string(tokens[i])).compare(g._alternativeOR_95Ls[j])==0)thisCohort._columnOR_95L=i;
                            }
                            //or 95U
                            for (unsigned int j=0; j<g._alternativeOR_95Us.size(); j++)
                            {
                                if (uc(string(tokens[i])).compare(g._alternativeOR_95Us[j])==0)thisCohort._columnOR_95U=i;
                            }
                            //imputed
                            for (unsigned int j=0; j<g._alternativeImputeds.size(); j++)
                            {
                                if (uc(string(tokens[i])).compare(g._alternativeImputeds[j])==0)thisCohort._columnImputed=i;
                            }
                            //filters
                            for (unsigned int j=0; j<g.filters.size(); j++)
                            {
                                if (uc(string(tokens[i])).compare(uc(g.filters[j]._columnName))==0)_filterColumnNr[j]=i;
                            }
                            
                        }
                        thisCohort._columnCount = n;
                        if (g._noAlleles)
                        {
                            thisCohort._columnEffectAllele=-9;
                            thisCohort._columnOtherAllele=-9;
                        }
                    }//header lines
                    else  //data lines
                    {
                        bool currentMarkerIsOK=true;
                        if (thisCohort._columnCount!=n)	// lets check if current row has right number of tokens
                        {
                            sprintf(sb, "E%.9d" , g._errNR); g._errNR++;
                            LOG <<  g.studyList[fileNr] << " row " << currentLine << " has different number of columns compared to header line ( " << n << ", " <<  thisCohort._columnCount << ") Skipping line!(" << sb << ")" << endl;
                            ERR << sb << " Error: "<< g.studyList[fileNr] << " row " << currentLine << " has different number of columns compared to header line ( " << n << ", " <<  thisCohort._columnCount << ") Skipping line!" << endl;
                            
                        }
                        if (currentLine==1)	//lets check if all necessary columns are present
                        {
                            if (thisCohort._columnMarkerName==-9 || 
                                (thisCohort._columnBeta==-9 && g._binaryTrait==false) ||
                                (thisCohort._columnSE==-9 && g._binaryTrait==false) ||
                                (thisCohort._columnOR==-9 && g._binaryTrait==true) ||
                                (thisCohort._columnOR_95L==-9 && g._binaryTrait==true) ||
                                (thisCohort._columnEffectAllele==-9 && g._noAlleles==false) ||
                                (thisCohort._columnOtherAllele==-9 && g._noAlleles==false))
                            {
                                cout << "\tmissing mandatory column! Skipping file." << endl;
                                sprintf(sb, "E%.9d" , g._errNR); g._errNR++;
                                LOG <<  g.studyList[fileNr] << " misses some of the mandatory columns. Skipping file!(" << sb << ")" << endl;
                                ERR << sb << " " << g.studyList[fileNr] << " misses some of the mandatory columns. Skipping file!" << endl;					
                                if (thisCohort._columnMarkerName==-9)ERR << sb << " Markername column is missing" << endl;
                                if (thisCohort._columnBeta==-9 && g._binaryTrait==false)ERR << sb << " Beta column is missing" << endl;
                                if (thisCohort._columnSE==-9 && g._binaryTrait==false)ERR << sb << " SE column is missing" << endl;
                                if (thisCohort._columnOR==-9 && g._binaryTrait==true)ERR << sb << " OR column is missing" << endl;
                                if (thisCohort._columnOR_95L==-9 && g._binaryTrait==true)ERR << sb << " OR _95L column is missing" << endl;
                                if (thisCohort._columnOR_95U==-9 && g._binaryTrait==true)ERR << sb << " OR _95U column is missing" << endl;
                                if (thisCohort._columnEffectAllele==-9 && g._noAlleles==false)ERR << sb << " Effect allele column is missing. If all effects are according to the same allele, then please use --no_alleles option" << endl;
                                if (thisCohort._columnOtherAllele==-9 && g._noAlleles==false)ERR << sb << " Other allele column is missing. If all effects are according to the same allele, then please use --no_alleles option" << endl;
                                return 0;
                            }
                            if (thisCohort._columnStrand==-9)
                            {
                                cout << "Strand column missing! Expecting always positive strand."<< endl;
                                LOG << "Strand column missing! Expecting always positive strand."<< endl;
                            }
                            for (int j=0;j<_filterColumnNr.size();j++)
                            {
                                if (_filterColumnNr[j]!=-1) cout << "Filtering by " << g.filters[j]._columnName << " column" << endl;
                            }
                            if (g._genomicControl)	//calculating lambdas for the file
                            {
                                bool lambdaSuccess;
                                if (g._binaryTrait==false)lambdaSuccess = getLambda(thisCohort._name, thisCohort._directLambda, 
                                    thisCohort._imputedLambda, thisCohort._directCount, thisCohort._imputedCount, thisCohort._columnBeta, thisCohort._columnSE, thisCohort._columnImputed, 
                                        thisCohort._columnCount, _filterColumnNr, g);
                                else lambdaSuccess = getLambda(thisCohort._name, thisCohort._directLambda, thisCohort._imputedLambda,
                                    thisCohort._directCount, thisCohort._imputedCount,thisCohort._columnOR, thisCohort._columnOR_95L, thisCohort._columnImputed, 
                                    thisCohort._columnCount, _filterColumnNr, g);
                                
                                if (lambdaSuccess)
                                {	
                                    cout << "GC lambda genotyped: " << thisCohort._directLambda << " (" << thisCohort._directCount<<  ") imputed: " << thisCohort._imputedLambda << " (" << thisCohort._imputedCount << ")"<< endl;
                                    
                                    LAMBDA_FILE << thisCohort._name;
                                    char sb [1024];		//temporary char buffer
                                    std::string x="";	//starting to collect all information into this string
                                    
                                    x.append("\t"); 
                                    sprintf(sb, "%.4f" , thisCohort._directLambda); 
                                    x.append(sb);
                                    
                                    x.append("\t"); 
                                    sprintf(sb, "%d" , thisCohort._directCount); 
                                    x.append(sb);
                                    
                                    x.append("\t"); 
                                    sprintf(sb, "%.4f" , thisCohort._imputedLambda); 
                                    x.append(sb);
                                    
                                    x.append("\t"); 
                                    sprintf(sb, "%d" , thisCohort._imputedCount); 
                                    x.append(sb);
                                    LAMBDA_FILE << x <<endl;
                                    
                                }
                            }
                            
                        }
                        //everything seems to be OK..lets read marker data
                        //marker name
                        string myMarker = string(tokens[thisCohort._columnMarkerName]);
                        _markerProblems.markersAll++;
                        //alleles
                        string myEffectAllele, myOtherAllele;
                        if (thisCohort._columnEffectAllele!=-9 && thisCohort._columnOtherAllele!=-9)
                        {
                            myEffectAllele = uc(string(tokens[thisCohort._columnEffectAllele]));
                            myOtherAllele = uc(string(tokens[thisCohort._columnOtherAllele]));
                            if (!checkAlleles(myEffectAllele, myOtherAllele) && !g._indel)	//problem with alleles - reporting
                            {
								sprintf(sb, "E%.9d" , g._errNR); g._errNR++;
								LOG <<  g.studyList[fileNr] << " has problem with alleles for marker " << myMarker << "!(" << sb << ")" << endl;
								ERR << sb << " " << g.studyList[fileNr] << " has problem with alleles for marker " << myMarker << "!"<< endl;
								ERR << sb << " " << g.studyList[fileNr] << " Given alleles: "<< string(tokens[thisCohort._columnEffectAllele]) << "/" << string(tokens[thisCohort._columnOtherAllele]) << ". Skipping marker!" << endl;
								currentMarkerIsOK=false;
								_markerProblems.wrongAlleles++;
                            }
                        }
                        else {myEffectAllele = "N";myOtherAllele = "N";}
                        //strand
                        bool myStrand = true;
                        if (thisCohort._columnStrand!=-9){if (string(tokens[thisCohort._columnStrand]).compare("-")==0){myStrand = false;}}
                        if (!myStrand && !g._indel)
                        {
                            myEffectAllele = flip(myEffectAllele);
                            myOtherAllele = flip(myOtherAllele);
                        }
                        else if (!myStrand & checkAlleles(myEffectAllele, myOtherAllele))
                        {
                            myEffectAllele = flip(myEffectAllele);
                            myOtherAllele = flip(myOtherAllele);
                        }
                        //eaf
                        double myEaf=-9;
                        if (thisCohort._columnEffectAlleleFreq!=-9)
                        {
                            if (atof(tokens[thisCohort._columnEffectAlleleFreq].c_str())>0 &&
                                atof(tokens[thisCohort._columnEffectAlleleFreq].c_str())<1)
                            {
                                myEaf = atof(tokens[thisCohort._columnEffectAlleleFreq].c_str());
                            }
                            else	//eaf out of range - report it
                            {
								sprintf(sb, "E%.9d" , g._errNR); g._errNR++;
								LOG <<  g.studyList[fileNr] << " has problem with effect allele frequency for marker " << myMarker << "!(" << sb << ")" << endl;
								ERR << sb << " " << g.studyList[fileNr] << " has problem with effect allele frequency " << myMarker << "!"<< endl;
								ERR << sb << " " << g.studyList[fileNr] << " Given value: "<< string(tokens[thisCohort._columnEffectAlleleFreq]) << ". Value not used!" << endl;					
                            }
                        }
                        //effect+se (or+ci)
                        double myBeta=-9;
                        double mySE=-9;
                        if (g._binaryTrait)
                        {
                            double oratio = atof(tokens[thisCohort._columnOR].c_str());
                            double oratio_95l = atof(tokens[thisCohort._columnOR_95L].c_str());
                            if (oratio>0 && oratio_95l< oratio && oratio_95l>0)
                            {
								myBeta=log(oratio);
								mySE = ((log(oratio)-log(oratio_95l))/1.96);
                            }
                            else // problem with or value
                            {
								sprintf(sb, "E%.9d" , g._errNR); g._errNR++;
								LOG <<  g.studyList[fileNr] << " has problem with odds ratio and its 95_L for marker " << myMarker << "!(" << sb << ")" << endl;
								ERR << sb << " " << g.studyList[fileNr] << " has problem with odds ratio or its CI of " << myMarker << "!"<< endl;
								ERR << sb << " " << g.studyList[fileNr] << " Given values: OR="<< string(tokens[thisCohort._columnOR]) << " CI_95L="<<string(tokens[thisCohort._columnOR_95L]) << ". Marker not used!" << endl;					
								_markerProblems.problemEffect++;
								currentMarkerIsOK=false;
                                
                            }
                        }
                        else
                        {
                            myBeta = atof(tokens[thisCohort._columnBeta].c_str());
                            mySE = atof(tokens[thisCohort._columnSE].c_str());
                            if (mySE<=0)
                            {
								sprintf(sb, "E%.9d" , g._errNR); g._errNR++;
								LOG <<  g.studyList[fileNr] << " has problem with beta and se for marker " << myMarker << "!(" << sb << ")" << endl;
								ERR << sb << " " << g.studyList[fileNr] << " has problem with odds ratio or its CI of " << myMarker << "!"<< endl;
								ERR << sb << " " << g.studyList[fileNr] << " Given values: BETA="<< string(tokens[thisCohort._columnBeta]) << " SE="<<string(tokens[thisCohort._columnSE]) << ". Marker not used!" << endl;					
								_markerProblems.problemEffect++;
								currentMarkerIsOK=false;
                            }
                        }
                        //imputed
                        bool myImputed = false;
                        if (thisCohort._columnImputed!=-9)
                        {
                            if (string(tokens[thisCohort._columnImputed])=="1" || string(tokens[thisCohort._columnImputed])=="0")
                            { 
                                if (atoi(tokens[thisCohort._columnImputed].c_str())==1){myImputed = true;}
                            }
                            else
                            {
								sprintf(sb, "E%.9d" , g._errNR); g._errNR++;
								LOG <<  g.studyList[fileNr] << " has problem with imputation status for marker " << myMarker << "!(" << sb << ")" << endl;
								ERR << sb << " " << g.studyList[fileNr] << " has problem with imputation status for " << myMarker << "!"<< endl;
								ERR << sb << " " << g.studyList[fileNr] << " Given values: IMPUTED="<< string(tokens[thisCohort._columnBeta]) << "Value must be 1-imputed, 0-directly genotyped." << endl;					
                            }
                        }
                        //n
                        int myN = -9;
                        if (thisCohort._columnN!=-9)
                        {
                            if (atoi(tokens[thisCohort._columnN].c_str())>0){myN = atoi(tokens[thisCohort._columnN].c_str());}
                            else 
                            {
								sprintf(sb, "E%.9d" , g._errNR); g._errNR++;
								LOG <<  g.studyList[fileNr] << " has problem with sample size for marker " << myMarker << "!(" << sb << ")" << endl;
								ERR << sb << " " << g.studyList[fileNr] << " has problem with sample size for " << myMarker << "!"<< endl;
								ERR << sb << " " << g.studyList[fileNr] << " Given values: N="<< string(tokens[thisCohort._columnN]) << "Value must be larger than zero." << endl;					
                            }
                        }
                        for (int j=0; j<_filterColumnNr.size();j++)
                        {
                            if (_filterColumnNr[j] != -1)
                            {
                                double valueX = atof(tokens[_filterColumnNr[j]].c_str());
                                if (isFiltered(valueX, g.filters[j]._equation, g.filters[j]._value))
                                {
                                    _markerFilteredBy[j]++;
                                    currentMarkerIsOK=0;
                                }
                            }
                        }
                        //all data read - will add marker to marker list
                        if (currentMarkerIsOK)
                        {
                            if (_markernames[myMarker]==0)
                            {
                                _markernames[myMarker]=_markerNamesCount;
                                _markerNamesCount++;
                                
                                if (markerPosition[myMarker]==0)
                                {
                                    markerlist.push_back(marker(myMarker, g._studyCount));
                                    if (markerlist[g._markerCount-1].addCohort(fileNr, g,
                                                                               thisCohort._directLambda, thisCohort._imputedLambda, 
                                                                               myMarker,myEffectAllele, myOtherAllele, myEaf,
                                                                               myStrand, myBeta, mySE, myImputed, myN, ERR, LOG, _markerProblems)==0)
                                        //										markerlist[markerCount-1].addStudy(currentimputed,
                                        //											currentstrand, currenteffectallele, currentnoneffectallele, currenteffectallelefreq,
                                        //											currentbeta, currentse, i, studies,n, ERR, errNR, LOG)	;
                                    {
                                        markerPosition[myMarker]=g._markerCount;
                                        g._markerCount++;
                                    }
                                }
                                else		//marker is already existing in database
                                {
                                    int x = markerPosition[myMarker];	//this is the line number of current marker
                                    markerlist[x-1].addCohort(fileNr, g,
                                                              thisCohort._directLambda, thisCohort._imputedLambda, 
                                                              myMarker,myEffectAllele, myOtherAllele, myEaf,
                                                              myStrand, myBeta, mySE, myImputed, myN, ERR, LOG, _markerProblems);
                                    
                                    //										markerlist[x-1].addStudy(currentimputed,
                                    //											currentstrand, currenteffectallele, currentnoneffectallele, currenteffectallelefreq,
                                    //											currentbeta, currentse, i, studies,n, ERR, errNR, LOG);
                                    
                                }
                            }
                            else
                            {
                                _markerProblems.problemMulti++;
                            }
                            _countGoodMarkers++;
                        }
                        else
                        {
                            _countBadMarkers++;
                        }
                        
                    }   //data lines end here
                    
                }
                if (int(currentLine/10000)==currentLine/10000.0)	cout << "Line: " << currentLine << "\r";
				currentLine++;
        }
        gzclose(F);
        
    }
    else
    {
    
    
		ifstream F (g.studyList[fileNr].c_str());
		if (F.is_open())
		{	
			while (! F.eof() )
        	{
			string line;
            vector<string> tokens;
            getline (F,line);
        	int n = Tokenize(string(line), tokens, " ");
        	if (n)
        	{
				if (currentLine==0) // get headers
				{
					for (int i=0; i<n; i++)
					{
						//markername
						for (unsigned int j=0; j<g._alternativeMarkerNames.size(); j++)
						{
							if (uc(string(tokens[i])).compare(g._alternativeMarkerNames[j])==0){thisCohort._columnMarkerName=i;}
						}
						//effect allele
						for (unsigned int j=0; j<g._alternativeEffectAlleles.size(); j++)
						{
							if (uc(string(tokens[i])).compare(g._alternativeEffectAlleles[j])==0)thisCohort._columnEffectAllele=i;
						}
						//other allele
						for (unsigned int j=0; j<g._alternativeOtherAlleles.size(); j++)
						{
							if (uc(string(tokens[i])).compare(g._alternativeOtherAlleles[j])==0)thisCohort._columnOtherAllele=i;
						}
						//effect allele freq
						for (unsigned int j=0; j<g._alternativeEffectAlleleFreqs.size(); j++)
						{
							if (uc(string(tokens[i])).compare(g._alternativeEffectAlleleFreqs[j])==0)thisCohort._columnEffectAlleleFreq=i;
						}
						//strand
						for (unsigned int j=0; j<g._alternativeStrands.size(); j++)
						{
							if (uc(string(tokens[i])).compare(g._alternativeStrands[j])==0)thisCohort._columnStrand=i;
						}
						//n
						for (unsigned int j=0; j<g._alternativeNs.size(); j++)
						{
							if (uc(string(tokens[i])).compare(g._alternativeNs[j])==0)thisCohort._columnN=i;
						}
						//beta
						for (unsigned int j=0; j<g._alternativeBetas.size(); j++)
						{
							if (uc(string(tokens[i])).compare(g._alternativeBetas[j])==0)thisCohort._columnBeta=i;
						}
						//se
						for (unsigned int j=0; j<g._alternativeSEs.size(); j++)
						{
							if (uc(string(tokens[i])).compare(g._alternativeSEs[j])==0)thisCohort._columnSE=i;
						}
						//or
						for (unsigned int j=0; j<g._alternativeORs.size(); j++)
						{
							if (uc(string(tokens[i])).compare(g._alternativeORs[j])==0)thisCohort._columnOR=i;
						}
						//or 95L
						for (unsigned int j=0; j<g._alternativeOR_95Ls.size(); j++)
						{
							if (uc(string(tokens[i])).compare(g._alternativeOR_95Ls[j])==0)thisCohort._columnOR_95L=i;
						}
						//or 95U
						for (unsigned int j=0; j<g._alternativeOR_95Us.size(); j++)
						{
							if (uc(string(tokens[i])).compare(g._alternativeOR_95Us[j])==0)thisCohort._columnOR_95U=i;
						}
						//imputed
						for (unsigned int j=0; j<g._alternativeImputeds.size(); j++)
						{
							if (uc(string(tokens[i])).compare(g._alternativeImputeds[j])==0)thisCohort._columnImputed=i;
						}
                        //filters
                        for (unsigned int j=0; j<g.filters.size(); j++)
                        {
                            if (uc(string(tokens[i])).compare(uc(g.filters[j]._columnName))==0)_filterColumnNr[j]=i;
                        }

					}
					thisCohort._columnCount = n;
					if (g._noAlleles)
					{
						thisCohort._columnEffectAllele=-9;
						thisCohort._columnOtherAllele=-9;
					}
				}//header lines
				else  //data lines
				{
					bool currentMarkerIsOK=true;
					if (thisCohort._columnCount!=n)	// lets check if current row has right number of tokens
					{
						sprintf(sb, "E%.9d" , g._errNR); g._errNR++;
						LOG <<  g.studyList[fileNr] << " row " << currentLine << " has different number of columns compared to header line ( " << n << ", " <<  thisCohort._columnCount << ") Skipping line!(" << sb << ")" << endl;
						ERR << sb << " Error: "<< g.studyList[fileNr] << " row " << currentLine << " has different number of columns compared to header line ( " << n << ", " <<  thisCohort._columnCount << ") Skipping line!" << endl;
						
					}
					if (currentLine==1)	//lets check if all necessary columns are present
					{
						if (thisCohort._columnMarkerName==-9 || 
							(thisCohort._columnBeta==-9 && g._binaryTrait==false) ||
							(thisCohort._columnSE==-9 && g._binaryTrait==false) ||
							(thisCohort._columnOR==-9 && g._binaryTrait==true) ||
							(thisCohort._columnOR_95L==-9 && g._binaryTrait==true) ||
							(thisCohort._columnEffectAllele==-9 && g._noAlleles==false) ||
							(thisCohort._columnOtherAllele==-9 && g._noAlleles==false))
						{
							cout << "\tmissing mandatory column! Skipping file." << endl;
							sprintf(sb, "E%.9d" , g._errNR); g._errNR++;
							LOG <<  g.studyList[fileNr] << " misses some of the mandatory columns. Skipping file!(" << sb << ")" << endl;
							ERR << sb << " " << g.studyList[fileNr] << " misses some of the mandatory columns. Skipping file!" << endl;					
							if (thisCohort._columnMarkerName==-9)ERR << sb << " Markername column is missing" << endl;
							if (thisCohort._columnBeta==-9 && g._binaryTrait==false)ERR << sb << " Beta column is missing" << endl;
							if (thisCohort._columnSE==-9 && g._binaryTrait==false)ERR << sb << " SE column is missing" << endl;
							if (thisCohort._columnOR==-9 && g._binaryTrait==true)ERR << sb << " OR column is missing" << endl;
							if (thisCohort._columnOR_95L==-9 && g._binaryTrait==true)ERR << sb << " OR _95L column is missing" << endl;
							if (thisCohort._columnOR_95U==-9 && g._binaryTrait==true)ERR << sb << " OR _95U column is missing" << endl;
							if (thisCohort._columnEffectAllele==-9 && g._noAlleles==false)ERR << sb << " Effect allele column is missing. If all effects are according to the same allele, then please use --no_alleles option" << endl;
							if (thisCohort._columnOtherAllele==-9 && g._noAlleles==false)ERR << sb << " Other allele column is missing. If all effects are according to the same allele, then please use --no_alleles option" << endl;
							return 0;
						}
						if (thisCohort._columnStrand==-9)
						{
							cout << "Strand column missing! Expecting always positive strand."<< endl;
							LOG << "Strand column missing! Expecting always positive strand."<< endl;
						}
                        for (int j=0;j<_filterColumnNr.size();j++)
                        {
                            if (_filterColumnNr[j]!=-1) cout << "Fltering by " << g.filters[j]._columnName << " column" << endl;
                        }
						if (g._genomicControl)	//calculating lambdas for the file
						{
							bool lambdaSuccess;
							if (g._binaryTrait==false)lambdaSuccess = getLambda(thisCohort._name, thisCohort._directLambda, 
								thisCohort._imputedLambda, thisCohort._directCount, thisCohort._imputedCount,
								thisCohort._columnBeta, thisCohort._columnSE, thisCohort._columnImputed, 
								thisCohort._columnCount,_filterColumnNr, g);
							else lambdaSuccess = getLambda(thisCohort._name, thisCohort._directLambda, thisCohort._imputedLambda,
								thisCohort._directCount, thisCohort._imputedCount,thisCohort._columnOR, thisCohort._columnOR_95L, thisCohort._columnImputed, 
								thisCohort._columnCount,_filterColumnNr, g);
							
							if (lambdaSuccess)
							{	
								cout << "GC lambda genotyped: " << thisCohort._directLambda << " (" << thisCohort._directCount<<  ") imputed: " << thisCohort._imputedLambda << " (" << thisCohort._imputedCount << ")"<< endl;
							
								LAMBDA_FILE << thisCohort._name;
								char sb [1024];		//temporary char buffer
								std::string x="";	//starting to collect all information into this string

								x.append("\t"); 
								sprintf(sb, "%.4f" , thisCohort._directLambda); 
								x.append(sb);

								x.append("\t"); 
								sprintf(sb, "%d" , thisCohort._directCount); 
								x.append(sb);

								x.append("\t"); 
								sprintf(sb, "%.4f" , thisCohort._imputedLambda); 
								x.append(sb);

								x.append("\t"); 
								sprintf(sb, "%d" , thisCohort._imputedCount); 
								x.append(sb);
								LAMBDA_FILE << x <<endl;

							}
						}
					
					}
					//everything seems to be OK..lets read marker data
					//marker name
					string myMarker = string(tokens[thisCohort._columnMarkerName]);
					_markerProblems.markersAll++;
					//alleles
					string myEffectAllele, myOtherAllele;
					if (thisCohort._columnEffectAllele!=-9 && thisCohort._columnOtherAllele!=-9)
					{
						myEffectAllele = uc(string(tokens[thisCohort._columnEffectAllele]));
						myOtherAllele = uc(string(tokens[thisCohort._columnOtherAllele]));
						if (!checkAlleles(myEffectAllele, myOtherAllele) && !g._indel)	//problem with alleles - reporting
						{
								sprintf(sb, "E%.9d" , g._errNR); g._errNR++;
								LOG <<  g.studyList[fileNr] << " has problem with alleles for marker " << myMarker << "!(" << sb << ")" << endl;
								ERR << sb << " " << g.studyList[fileNr] << " has problem with alleles for marker " << myMarker << "!"<< endl;
								ERR << sb << " " << g.studyList[fileNr] << " Given alleles: "<< string(tokens[thisCohort._columnEffectAllele]) << "/" << string(tokens[thisCohort._columnOtherAllele]) << ". Skipping marker!" << endl;
								currentMarkerIsOK=false;
								_markerProblems.wrongAlleles++;
						}
					}
					else {myEffectAllele = "N";myOtherAllele = "N";}
					//strand
					bool myStrand = true;
					if (thisCohort._columnStrand!=-9){if (string(tokens[thisCohort._columnStrand]).compare("-")==0){myStrand = false;}}
					if (!myStrand && !g._indel)
					{
						myEffectAllele = flip(myEffectAllele);
						myOtherAllele = flip(myOtherAllele);
					}
					//eaf
					double myEaf=-9;
					if (thisCohort._columnEffectAlleleFreq!=-9)
					{
						if (atof(tokens[thisCohort._columnEffectAlleleFreq].c_str())>0 &&
							atof(tokens[thisCohort._columnEffectAlleleFreq].c_str())<1)
						{
							myEaf = atof(tokens[thisCohort._columnEffectAlleleFreq].c_str());
						}
						else	//eaf out of range - report it
						{
								sprintf(sb, "E%.9d" , g._errNR); g._errNR++;
								LOG <<  g.studyList[fileNr] << " has problem with effect allele frequency for marker " << myMarker << "!(" << sb << ")" << endl;
								ERR << sb << " " << g.studyList[fileNr] << " has problem with effect allele frequency " << myMarker << "!"<< endl;
								ERR << sb << " " << g.studyList[fileNr] << " Given value: "<< string(tokens[thisCohort._columnEffectAlleleFreq]) << ". Value not used!" << endl;					
						}
					}
					//effect+se (or+ci)
					double myBeta=-9;
					double mySE=-9;
					if (g._binaryTrait)
					{
						double oratio = atof(tokens[thisCohort._columnOR].c_str());
						double oratio_95l = atof(tokens[thisCohort._columnOR_95L].c_str());
						if (oratio>0 && oratio_95l< oratio && oratio_95l>0)
						{
								myBeta=log(oratio);
								mySE = ((log(oratio)-log(oratio_95l))/1.96);
						}
						else // problem with or value
						{
								sprintf(sb, "E%.9d" , g._errNR); g._errNR++;
								LOG <<  g.studyList[fileNr] << " has problem with odds ratio and its 95_L for marker " << myMarker << "!(" << sb << ")" << endl;
								ERR << sb << " " << g.studyList[fileNr] << " has problem with odds ratio or its CI of " << myMarker << "!"<< endl;
								ERR << sb << " " << g.studyList[fileNr] << " Given values: OR="<< string(tokens[thisCohort._columnOR]) << " CI_95L="<<string(tokens[thisCohort._columnOR_95L]) << ". Marker not used!" << endl;					
								_markerProblems.problemEffect++;
								currentMarkerIsOK=false;

						}
					}
					else
					{
						myBeta = atof(tokens[thisCohort._columnBeta].c_str());
						mySE = atof(tokens[thisCohort._columnSE].c_str());
						if (mySE<=0)
						{
								sprintf(sb, "E%.9d" , g._errNR); g._errNR++;
								LOG <<  g.studyList[fileNr] << " has problem with beta and se for marker " << myMarker << "!(" << sb << ")" << endl;
								ERR << sb << " " << g.studyList[fileNr] << " has problem with odds ratio or its CI of " << myMarker << "!"<< endl;
								ERR << sb << " " << g.studyList[fileNr] << " Given values: BETA="<< string(tokens[thisCohort._columnBeta]) << " SE="<<string(tokens[thisCohort._columnSE]) << ". Marker not used!" << endl;					
								_markerProblems.problemEffect++;
								currentMarkerIsOK=false;
						}
					}
					//imputed
					bool myImputed = false;
					if (thisCohort._columnImputed!=-9)
					{
						if (string(tokens[thisCohort._columnImputed])=="1" || string(tokens[thisCohort._columnImputed])=="0")
						{ 
						if (atoi(tokens[thisCohort._columnImputed].c_str())==1){myImputed = true;}
						}
						else
						{
								sprintf(sb, "E%.9d" , g._errNR); g._errNR++;
								LOG <<  g.studyList[fileNr] << " has problem with imputation status for marker " << myMarker << "!(" << sb << ")" << endl;
								ERR << sb << " " << g.studyList[fileNr] << " has problem with imputation status for " << myMarker << "!"<< endl;
								ERR << sb << " " << g.studyList[fileNr] << " Given values: IMPUTED="<< string(tokens[thisCohort._columnBeta]) << "Value must be 1-imputed, 0-directly genotyped." << endl;					
						}
					}
					//n
					int myN = -9;
					if (thisCohort._columnN!=-9)
					{
						if (atoi(tokens[thisCohort._columnN].c_str())>0){myN = atoi(tokens[thisCohort._columnN].c_str());}
						else 
						{
								sprintf(sb, "E%.9d" , g._errNR); g._errNR++;
								LOG <<  g.studyList[fileNr] << " has problem with sample size for marker " << myMarker << "!(" << sb << ")" << endl;
								ERR << sb << " " << g.studyList[fileNr] << " has problem with sample size for " << myMarker << "!"<< endl;
								ERR << sb << " " << g.studyList[fileNr] << " Given values: N="<< string(tokens[thisCohort._columnN]) << "Value must be larger than zero." << endl;					
						}
					}
                    for (int j=0; j<_filterColumnNr.size();j++)
                    {
                        if (_filterColumnNr[j] != -1)
                        {
                            double valueX = atof(tokens[_filterColumnNr[j]].c_str());
                            if (isFiltered(valueX, g.filters[j]._equation, g.filters[j]._value))
                            {
                                _markerFilteredBy[j]++;
                                currentMarkerIsOK=0;
                            }
                        }
                    }
					//all data read - will add marker to marker list
					if (currentMarkerIsOK)
					{
						if (_markernames[myMarker]==0)
						{
							_markernames[myMarker]=_markerNamesCount;
							_markerNamesCount++;

							if (markerPosition[myMarker]==0)
							{
								markerlist.push_back(marker(myMarker, g._studyCount));
								if (markerlist[g._markerCount-1].addCohort(fileNr, g,
									thisCohort._directLambda, thisCohort._imputedLambda, 
									myMarker,myEffectAllele, myOtherAllele, myEaf,
									myStrand, myBeta, mySE, myImputed, myN, ERR, LOG, _markerProblems)==0)
		//										markerlist[markerCount-1].addStudy(currentimputed,
		//											currentstrand, currenteffectallele, currentnoneffectallele, currenteffectallelefreq,
		//											currentbeta, currentse, i, studies,n, ERR, errNR, LOG)	;
								{
									markerPosition[myMarker]=g._markerCount;
									g._markerCount++;
								}
							}
							else		//marker is already existing in database
							{
								int x = markerPosition[myMarker];	//this is the line number of current marker
								markerlist[x-1].addCohort(fileNr, g,
									thisCohort._directLambda, thisCohort._imputedLambda, 
									myMarker,myEffectAllele, myOtherAllele, myEaf,
									myStrand, myBeta, mySE, myImputed, myN, ERR, LOG, _markerProblems);

		//										markerlist[x-1].addStudy(currentimputed,
		//											currentstrand, currenteffectallele, currentnoneffectallele, currenteffectallelefreq,
		//											currentbeta, currentse, i, studies,n, ERR, errNR, LOG);

							}
						}
						else
						{
							_markerProblems.problemMulti++;
						}
						_countGoodMarkers++;
					}
					else
					{
						_countBadMarkers++;
					}
					
				}   //data lines end here
			
			}
			if (int(currentLine/10000)==currentLine/10000.0)	cout << "Line: " << currentLine << "\r";
				currentLine++;
				}//while lines
				F.close();
		} //if is open
		else
		{
			cerr << "Cannot open file. Exit program!"<< endl;
			exit(1);
		}
    }
		cout << "Marker count: " << _markerProblems.markersAll << " Markers passing sanity check (and filters): " << _markerProblems.markersOK << endl;
		cout << "Strand problems: " << _markerProblems.problemStrand << " Wrong alleles: " << _markerProblems.wrongAlleles << endl;
		cout << "Effect problems: " << _markerProblems.problemEffect << " Multiple occurances: " << _markerProblems.problemMulti << endl;
        for (int j=0; j<_filterColumnNr.size();j++)
        {   
            if (_filterColumnNr[j] != -1)
            {
                cout << "FILTER " << g.filters[j]._columnName << " " << g.filters[j]._equation << " "  << g.filters[j]._value << " :" << _markerFilteredBy[j] << endl;
            }
        }
    
		return true;
}

