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

#include "cohort.h"
#include <iostream>
#include <stdlib.h>
using namespace std;
cohort::cohort(void)
{
	_name = "N";
	_columnMarkerName = -9;
	_columnEffectAllele = -9;
	_columnOtherAllele = -9;
	_columnEffectAlleleFreq = -9;
	_columnStrand = -9;
	_columnBeta = -9;
	_columnSE = -9;
	_columnOR = -9;
	_columnOR_95L = -9;
	_columnOR_95U = -9;
	_columnN = -9;
	_columnImputed = -9;
	_directLambda = 1;
	_imputedLambda = 1;
	_directCount = 0;
	_imputedCount = 0;
	_columnCount = -9;
}


