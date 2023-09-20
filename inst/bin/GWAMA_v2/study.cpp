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
#include "study.h"
#include <algorithm>

using namespace std;
study::study(std::string name)
{
	_name = name;
	_imputedLambda = 1;
	_directLambda = 1;
	_imputedCount = 0;
	_directCount = 0;
}

study::~study(void)
{
}

std::string
study::getName(void)
{
	return _name;

}


void 
study::addImputed(double x)
{
	_imputedChiStat.push_back(x);
	_imputedCount ++;
}

void 
study::addDirect(double x)
{
	_directChiStat.push_back(x);
	_directCount ++;

}
bool 
study::calculateLambdas(void)
{
	if (_imputedCount>0)
	{
		double _medianI = 0;
		study::sortVec( _imputedChiStat, _imputedCount);
		if (_imputedCount%2!=0) _medianI = _imputedChiStat[((_imputedCount+1)/2)-1];
		else _medianI = (_imputedChiStat[_imputedCount/2 -1] + _imputedChiStat[_imputedCount/2])/2;
		_imputedLambda = _medianI/0.4549364;		//median of chi-sq from R ... qchisq(0.5, df= 1)

	}

	if (_directCount>0)
	{
		double _medianD = 0;
		study::sortVec( _directChiStat, _directCount);
		if (_directCount%2!=0) _medianD = _directChiStat[((_directCount+1)/2)-1];
		else _medianD = (_directChiStat[_directCount/2 -1] + _directChiStat[_directCount/2])/2;
		_directLambda = _medianD/0.4549364;        //median of chi-sq from R ... qchisq(0.5, df= 1)
	}
	return true;
}
double 
study::getImputedLambda(void)
{
	return _imputedLambda;
}

double 
study::getDirectLambda(void)
{
	return _directLambda;
}

int
study::getImputedCount(void)
{
	return _imputedCount;
}

int 
study::getDirectCount(void)
{
	return _directCount;
}

void
study::printStudy()
{
	cout << "Study: " << _name.c_str() << endl;

}

void
study::sortVec(vector <double>& x, int size)
{
	 std::sort(x.begin(), x.end());
}
