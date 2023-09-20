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
#include <cctype> // std::toupper

#include "global.h"
#include "tools.h"

using namespace std;

void printVersion(string version)
{
   cout << "GWAMA version: " << version << endl;
   exit(0);
}
void printHelp(string version)
{
   string v;
   double len =17 -  version.size()-1;
   int f = (int) ((len/2.0)+0.6);
   int r = (int) ((len/2.0));
   for (int i=0; i < f;i++)v+=" ";
   v+='v';
   v+=version;
   for (int i=0; i < r;i++)v+=" ";
   cout << endl;
   cout << "@----------------------------------------------------------@" << endl;
   cout << "|      GWAMA        |" << v << "|    May, 2010       |" << endl;
   cout << "|----------------------------------------------------------|" << endl;
   cout << "|           (C) 2008 Reedik Magi & Andrew P Morris         |" << endl;
   cout << "|                    BSD 3-Clause License                  |" << endl;
   cout << "|----------------------------------------------------------|" << endl;
   cout << "|  For documentation, citation & bug-report instructions:  |" << endl;
   cout << "|              http://www.well.ox.ac.uk/GWAMA/             |" << endl;
   cout << "@----------------------------------------------------------@" << endl;
   cout << endl;	     
   cout << endl;	     
   cout << "Command line options:"<<endl;
   cout << endl;
   cout << "GWAMA --filelist {filename} or -i {filename} Specify     "<<endl;
   cout << "                    studies' result files. Default =     "<<endl;
   cout << "                    gwama.in                             "<<endl<<endl; 
   cout << "    --output {fileroot} or -o  {fileroot} Specify file   "<<endl;
   cout << "                    root for output of analysis. Default "<<endl;
   cout << "                    <gwama> (gwama.out)                  "<<endl<<endl;
   cout << "    --random or -r  Use random effect correction         "<<endl<<endl;;
   cout << "    --genomic_control or -gc  Use genomic control for    "<<endl;
   cout << "                    adjusting studies' result files      "<<endl<<endl;
   cout << "    --genomic_control_output or -gco      Use  genomic   "<<endl;
   cout << "                    control on meta-analysis summary.    "<<endl;
   cout << "                    (i.e. results  of meta-              "<<endl;
   cout << "                    analysis are corrected for gc)       "<<endl<<endl;
   cout << "    --quantitative or -qt    Use this option, if trait is"<<endl;
   cout << "                    quantitative (columns BETA & SE).    "<<endl;
   cout << "                    Default is binary trait (columns OR, "<<endl;
   cout << "                    OR95_U, OR_95_L)                     "<<endl<<endl;
   cout << "    --threshold {0-1} or -t {0-1} The p-value threshold  "<<endl;
   cout << "                    for showing direction. Default = 1   "<<endl<<endl;
   cout << "    --map {filename} or -m {filename} Marker position    "<<endl;
   cout << "                    file for chromosome and position info"<<endl<<endl;
   cout << "    --no_alleles    No allele information has been given."<<endl;
   cout << "                    Expecting always the same EA.        "<<endl<<endl;
   cout << "    --indel_alleles   Allele labes might contain more    "<<endl;
   cout << "                    than single letter. No strand checks."<<endl<<endl;    
   cout << "    --sex           Run gender-differentiated and gender-"<<endl;
   cout << "                    heterogeneity analysis. Gender info  "<<endl;
   cout << "                    must be provided in filelist file.   "<<endl;
   cout << "                    (second column after file names is   "<<endl;
   cout << "                    either M or F) More info in tutorial."<<endl<<endl;
   cout << "    --filter        Set a filtering based on column      "<<endl;    
   cout << "                    header. It needs 3 arguments: column "<<endl; 
   cout << "                    name, equation [>,<,>=,<=,==,!=],    "<<endl;  
   cout << "                    filter value. Multiple filters can be"<<endl; 
   cout << "                    set.                                 "<<endl<<endl;    
   cout << "    --name_marker   alternative header to marker name col"<<endl<<endl;
   cout << "    --name_strand   alternative header to strand column  "<<endl<<endl;
   cout << "    --name_n        alternative header to sample size col"<<endl<<endl;
   cout << "    --name_eaf      alternative header to EAF column     "<<endl<<endl;
   cout << "    --name_beta     alternative header to beta column    "<<endl<<endl;
   cout << "    --name_se       alternative header to std. err. col  "<<endl<<endl;
   cout << "    --name_or	    alternative header to OR column      "<<endl<<endl;
   cout << "    --name_or_95l   alternative header to OR 95L column  "<<endl<<endl;
   cout << "    --name_or_95u   alternative header to OR 95U column  "<<endl<<endl;
   cout << "    --help or -h    Print this help                      "<<endl<<endl;
   cout << "    --version or -v  Print GWAMA version number          "<<endl<<endl;   
   exit (0);
}


bool
readCommandline (int argc, char *argv[], global & _g)
{
	for (int i = 1; i < argc; i++)
	{
		if (string(argv[i]).compare("--version")==0 || string(argv[i]).compare("-v")==0)
		 {
			printVersion(_g._version);
			return 1;
		 }
		else if (string(argv[i]).compare("--help")==0 || string(argv[i]).compare("-h")==0)
		 {
			printHelp(_g._version);
			return 1;
		 }
		else if (string(argv[i]).compare("--genomic_control")==0 || string(argv[i]).compare("-gc")==0)
		 {
			_g._genomicControl=true;
			cout << "Genomic control enabled"  << endl;
		 }
		else if (string(argv[i]).compare("--no_alleles")==0)
		 {
			_g._noAlleles=true;
			cout << "All effects are expected to be according to the same allele."  << endl;
		 }

		else if (string(argv[i]).compare("--random")==0 || string(argv[i]).compare("-r")==0)
		 {
			_g._randomEffect=true;
			cout << "Using random effect correction"  << endl;
		 }
		else if (string(argv[i]).compare("--quantitative")==0 || string(argv[i]).compare("-qt")==0)
		 {
			 _g._binaryTrait = false;
			cout << "Quantitative trait (BETA+SE)"  << endl;

		 }
		else if (string(argv[i]).compare("--genomic_control_output")==0 || string(argv[i]).compare("-gco")==0)
		 {
			 _g._genomicControlOutput=true;
			cout << "Genomic control for meta-analysis output enabled"  << endl;
		 }
		else if (string(argv[i]).compare("--sex")==0)
		 {
			 _g._genderHet=true;
			cout << "Gender-specific heterogeneity analysis enabled"  << endl;
		 }
        else if (string(argv[i]).compare("--indel_alleles")==0)
        {
            _g._indel=true;
			cout << "Indel allele names enabled. No strand checks"  << endl;
        }
		else if (string(argv[i]).compare("--name_marker")==0)
		{
			if (i != argc-1) 
			{
				cout<<"Alternative marker column name: " << uc(string(argv[i+1])) << endl;; 
				_g._alternativeMarkerNames.push_back(uc(string(argv[i+1])));
			}
			if (i == argc-1 ) 
			{
				cout<<"ERROR: Missing a value for --name_marker\n"; 
				return 0;
			}	
			i++;
		}
		else if (string(argv[i]).compare("--name_ea")==0)
		{
			if (i != argc-1) 
			{
				cout<<"Alternative effect allele column name: " << uc(string(argv[i+1])) << endl;; 
				_g._alternativeEffectAlleles.push_back(uc(string(argv[i+1])));
			}
			if (i == argc-1 ) 
			{
				cout<<"ERROR: Missing a value for --name_ea\n"; 
				return 0;
			}	
			i++;
		}
		else if (string(argv[i]).compare("--name_nea")==0)
		{
			if (i != argc-1) 
			{
				cout<<"Alternative other allele column name: " << uc(string(argv[i+1])) << endl;; 
				_g._alternativeOtherAlleles.push_back(uc(string(argv[i+1])));
			}
			if (i == argc-1 ) 
			{
				cout<<"ERROR: Missing a value for --name_nea\n"; 
				return 0;
			}	
			i++;
		}
		else if (string(argv[i]).compare("--name_or")==0)
		{
			if (i != argc-1) 
			{
				cout<<"Alternative OR column name: " << uc(string(argv[i+1])) << endl;; 
				_g._alternativeORs.push_back(uc(string(argv[i+1])));
			}
			if (i == argc-1 ) 
			{
				cout<<"ERROR: Missing a value for --name_or\n"; 
				return 0;
			}	
			i++;
		}
		else if (string(argv[i]).compare("--name_or_95l")==0)
		{
			if (i != argc-1) 
			{
				cout<<"Alternative lower CI of OR column name: " << uc(string(argv[i+1])) << endl;; 
				_g._alternativeOR_95Ls.push_back(uc(string(argv[i+1])));
			}
			if (i == argc-1 ) 
			{
				cout<<"ERROR: Missing a value for --name_or_95l\n"; 
				return 0;
			}	
			i++;
		}
		else if (string(argv[i]).compare("--name_or_95u")==0)
		{
			if (i != argc-1) 
			{
				cout<<"Alternative upper CI of OR column name: " << uc(string(argv[i+1])) << endl;; 
				_g._alternativeOR_95Us.push_back(uc(string(argv[i+1])));
			}
			if (i == argc-1 ) 
			{
				cout<<"ERROR: Missing a value for --name_or_95u\n"; 
				return 0;
			}	
			i++;
		}
		else if (string(argv[i]).compare("--name_beta")==0)
		{
			if (i != argc-1) 
			{
				cout<<"Alternative beta column name: " << uc(string(argv[i+1])) << endl;; 
				_g._alternativeBetas.push_back(uc(string(argv[i+1])));
			}
			if (i == argc-1 ) 
			{
				cout<<"ERROR: Missing a value for --name_beta\n"; 
				return 0;
			}	
			i++;
		}
		else if (string(argv[i]).compare("--name_se")==0)
		{
			if (i != argc-1) 
			{
				cout<<"Alternative std. error column name: " << uc(string(argv[i+1])) << endl;; 
				_g._alternativeSEs.push_back(uc(string(argv[i+1])));
			}
			if (i == argc-1 ) 
			{
				cout<<"ERROR: Missing a value for --name_se\n"; 
				return 0;
			}	
			i++;
		}
		else if (string(argv[i]).compare("--name_n")==0)
		{
			if (i != argc-1) 
			{
				cout<<"Alternative sample size column name: " << uc(string(argv[i+1])) << endl;; 
				_g._alternativeNs.push_back(uc(string(argv[i+1])));
			}
			if (i == argc-1 ) 
			{
				cout<<"ERROR: Missing a value for --name_n\n"; 
				return 0;
			}	
			i++;
		}
		else if (string(argv[i]).compare("--name_eaf")==0)
		{
			if (i != argc-1) 
			{
				cout<<"Alternative effect allele frequency column name: " << uc(string(argv[i+1])) << endl;; 
				_g._alternativeEffectAlleleFreqs.push_back(uc(string(argv[i+1])));
			}
			if (i == argc-1 ) 
			{
				cout<<"ERROR: Missing a value for --name_eaf\n"; 
				return 0;
			}	
			i++;
		}
		else if (string(argv[i]).compare("--name_strand")==0)
		{
			if (i != argc-1) 
			{
				cout<<"Alternative strand column name: " << uc(string(argv[i+1])) << endl;; 
				_g._alternativeStrands.push_back(uc(string(argv[i+1])));
			}
			if (i == argc-1 ) 
			{
				cout<<"ERROR: Missing a value for --name_strand\n"; 
				return 0;
			}	
			i++;
		}

		else if (string(argv[i]).compare("--filelist")==0 || string(argv[i]).compare("-i")==0)
		{
			if (i != argc-1) 
			{
				_g._fileList = string(argv[i+1]);
			}
			if (i == argc-1 ) 
			{
				cout<<"ERROR: Missing a value for --filelist\n"; 
				return 0;
			}	
			i++;
		}
		else if (string(argv[i]).compare("--output")==0 || string(argv[i]).compare("-o")==0)
		{
			if (i != argc-1) 
			{
				_g._outputRoot = string(argv[i+1]);
			}
			if (i == argc-1 ) 
			{
				cout<<"ERROR: Missing a value for --output\n"; 
				return 0;
			}	
			i++;
		}
		else if (string(argv[i]).compare("--map")==0 || string(argv[i]).compare("-m")==0)
		{
			if (i != argc-1) 
			{
				_g._markerMap = string(argv[i+1]);
			}
			if (i == argc-1 ) 
			{
				cout<<"ERROR: Missing a value for --map\n"; 
				return 0;
			}	
			i++;
		}
        else if (string(argv[i]).compare("--filter")==0)
		{
            string a = "";
            string b = "";
            double c = -9999;
			if (i != argc-1) 
			{
				a = string(argv[i+1]);
			}
			if (i == argc-1 ) 
			{
				cout<<"ERROR: Missing a value for --filter\n"; 
				return 0;
			}	
            i++;
            if (i != argc-1) 
			{
				b = string(argv[i+1]);
			}
			if (i == argc-1 ) 
			{
				cout<<"ERROR: Missing a value for --filter\n"; 
				return 0;
			}	
            i++;
			if (i != argc-1) 
			{
				c = atof(argv[i+1]);
			}
			if (i == argc-1 ) 
			{
				cout<<"ERROR: Missing a value for --filter\n"; 
				return 0;
			}	
            i++;

            filter _filter(a,b,c);
            if (! _filter.valid())
            {
                cout<<"ERROR: Filter not valid. Please check the format\n"; 
				return 0; 
            }
            cout << "Setting filter for column: " << a << " to filter out samples with values " << b << " " << c << endl;
            _g.filters.push_back(_filter);
		}
		else if (string(argv[i]).compare("--threshold")==0 || string(argv[i]).compare("-t")==0)
		{
			if (i != argc-1)
			{
				if (atof(argv[i+1]) >0 && atof(argv[i+1]) <= 1)
				_g._thresholdPValDir = atof(argv[i+1]);
			}
			if (i == argc-1) 
			{
				cout<<"ERROR: Missing or errogeneous value for --threshold\n"; 
				return 0;
			}	
			i++;
		}
		else
		{
				cout<<"Unknown command line option: " << string(argv[i]) << ". Exit program.\n"; 
			return false;
		}
	}
	return true;
}
