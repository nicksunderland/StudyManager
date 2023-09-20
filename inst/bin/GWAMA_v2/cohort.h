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

#ifndef _COHORT_H_
#define _COHORT_H_
#include <iostream>
#include <stdlib.h>
using namespace std;
class cohort
{
public:
	cohort(void);
	string _name;
	int _columnMarkerName;
	int _columnEffectAllele;
	int _columnOtherAllele;
	int _columnEffectAlleleFreq;
	int _columnStrand;
	int _columnBeta;
	int _columnSE;
	int _columnOR;
	int _columnOR_95L;
	int _columnOR_95U;
	int _columnN;
	int _columnImputed;
	int _columnCount;
	double _directLambda;
	double _imputedLambda;
	int _directCount;
	int _imputedCount;


};
#endif

