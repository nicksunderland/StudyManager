/*************************************************************************
GWAMA software:  May, 2009

Contributors:
    * Lauris Kaplinski lauris@kaplinski.com
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
#include <math.h>
#include <cctype> // std::toupper
#include <zlib.h>
#include <string>
#include <algorithm>
#include "tools.h"
using namespace std;



int Tokenize(const  string& str1,
                      vector<string>& tokens,
			const string& delimiters = " ")
{
    int cnt = 0;
    string emptyStr = "";
    string str = str1;
    std::replace( str.begin(), str.end(), '\r', ' ' );
    std::replace( str.begin(), str.end(), '\t', ' ' );



    // Skip delimiters at beginning.
    string::size_type lastPos = str.find_first_not_of(delimiters, 0);
    // Find first "non-delimiter".
    string::size_type pos     = str.find_first_of(delimiters, lastPos);

    while (string::npos != pos || string::npos != lastPos )
    {
	if (str.substr(lastPos, pos - lastPos) != emptyStr)
	{
        	// Found a token, add it to the vector.
        	tokens.push_back(str.substr(lastPos, pos - lastPos));
        	// Skip delimiters.  Note the "not_of"
        	lastPos = str.find_first_not_of(delimiters, pos);
        	// Find next "non-delimiter"
        	pos = str.find_first_of(delimiters, lastPos);
		//
		cnt++;
	}
    }
    return cnt;
}
void
sortVec(vector <double>& x, int size)
{
	 std::sort(x.begin(), x.end());
}



string uc(string s)
{
  const int length = s.length();
  for(int i=0; i!=length ; ++i)
  {
    s[i] = std::toupper(s[i]);
  }
  return s;
}

bool checkAlleles(string & s1, string & s2)
{
	if (s1.compare("1")==0) s1 = "A";
	if (s1.compare("2")==0) s1 = "C";
	if (s1.compare("3")==0) s1 = "G";
	if (s1.compare("4")==0) s1 = "T";
	if (s2.compare("1")==0) s2 = "A";
	if (s2.compare("2")==0) s2 = "C";
	if (s2.compare("3")==0) s2 = "G";
	if (s2.compare("4")==0) s2 = "T";
	if (!(s1.compare("A")==0 || s1.compare("C")==0 || s1.compare("G")==0 || s1.compare("T")==0))return false;
	if (!(s2.compare("A")==0 || s2.compare("C")==0 || s2.compare("G")==0 || s2.compare("T")==0))return false;
	if (s1==s2)return false;
	return true;
}
string flip(string s)
{
	if (s.compare("A")==0) return "T";
	if (s.compare("C")==0) return "G";
	if (s.compare("G")==0) return "C";
	if (s.compare("T")==0) return "A";
	return "N";
}
