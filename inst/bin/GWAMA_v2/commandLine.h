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
#ifndef _COMMANDLINE_H_
#define _COMMANDLINE_H_


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
using namespace std;
bool readCommandline(int argc, char *argv[], global & _g);
void printVersion(string version);
void printHelp(string version);

#endif

