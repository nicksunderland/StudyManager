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

using namespace std;
class study
{
private:
	int _imputedCount;
	int _directCount;
	std::string _name;
	vector <double> _imputedChiStat;
	vector <double> _directChiStat;
	double _imputedLambda;
	double _directLambda;

public:
	study(std::string name);
	~study(void);
	void addImputed(double x);
	void addDirect(double x);
	bool calculateLambdas(void);
	double getImputedLambda(void);
	double getDirectLambda(void);

	int getImputedCount(void);
	int getDirectCount(void);

	void printStudy(void);
	void sortVec(vector <double> &vec, int size);
	std::string getName(void);
};
