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
#ifndef _READFILE_H_
#define _READFILE_H_

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


bool readFilelist(global & _g, ofstream & ERR,  ofstream & LOG);
bool
getLambda(string name, double & directLambda, double & imputedLambda, int & directCount, int & imputedCount, int columnBeta, int columnSE, int columnImputed, int columnCount, vector < int> filterColumns, global & g);
bool
readCohort(int fileNr, global & g,map<string, int> & markerPosition,vector <marker> & markerlist, ofstream & ERR, ofstream & LOG, ofstream & LAMBDA_FILE);
bool 
readMapFile(string markermap, vector <marker> & markerlist,
				 map<string, int> & markerPosition,
				 ofstream & ERR, int & errNR, ofstream & LOG);


#endif
