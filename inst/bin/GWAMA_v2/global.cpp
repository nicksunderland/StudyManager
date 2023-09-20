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
#include "global.h"
#include <iostream>
#include <map>
#include <vector>
#include <stdlib.h>
#include <fstream>
using namespace std;
filter::filter(string columnName, string equation, double value)
{
    _columnName = columnName;
    _equation = equation;
    _value = value;
}
bool
filter::valid()
{
    if (_columnName == "" || _value == -9999) return false;
    if (_equation == "==" || _equation == "!=" || _equation == ">"
         || _equation == ">="  || _equation == "<"  || _equation == "<=") return true;
    return false;
}
global::global()
{
	_errNR =				1;				//count of errors and warnings found
	_thresholdPValDir =		1;				//p-value threshold for directions
	_markerCount =			1;
	_studyCount =			0;
	_genomicControl =		false;			//genomic control for input files
	_genomicControlOutput = false;			//genomic control for output
	_binaryTrait =			true;			//by default binary trait
	_randomEffect =			false;		    //use random effect
	_genderHet =			false;			//use gender-specific analysis
	_noAlleles =			false;			//dont use EA and NEA columns
    _indel =                false;          //dont use long allele names
	_fileList =				"gwama.in";		//gwama.in file
	_markerMap =			"N";			//map file name
	_outputRoot =			"gwama";	    //output file root (w/o file extension)
	_version =				"2.2.2";

	outputLambda=1;
	outputLambda_male=1;
	outputLambda_female=1;



//	_alternativeMarkerNames.push_back("SNP");
//	_alternativeMarkerNames.push_back("RS");
//	_alternativeMarkerNames.push_back("SNPID");
//	_alternativeMarkerNames.push_back("ID");
//	_alternativeMarkerNames.push_back("MARKER");
	_alternativeMarkerNames.push_back("MARKERNAME");
	_alternativeStrands.push_back("STRAND");
	_alternativeEffectAlleles.push_back("EA");
//	_alternativeEffectAlleles.push_back("EFFECT_ALLELE");
//	_alternativeOtherAlleles.push_back("NON_EFFECT_ALLELE");
//	_alternativeOtherAlleles.push_back("OTHER_ALLELE");
	_alternativeOtherAlleles.push_back("NEA");
//	_alternativeEffectAlleleFreqs.push_back("EFFECT_ALLELE_FREQUENCY");
	_alternativeEffectAlleleFreqs.push_back("EAF");
//	_alternativeEffectAlleleFreqs.push_back("EFFECT_ALLELE_FREQ");
	_alternativeBetas.push_back("BETA");
//	_alternativeBetas.push_back("EFFECT");
	_alternativeSEs.push_back("SE");
//	_alternativeSEs.push_back("STDERR");
	_alternativeORs.push_back("OR");
	_alternativeOR_95Ls.push_back("OR_95L");
	_alternativeOR_95Us.push_back("OR_95U");
	_alternativeNs.push_back("N");
	_alternativeImputeds.push_back("IMPUTED");
}
