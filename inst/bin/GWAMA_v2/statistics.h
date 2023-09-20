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


#define __BDM_STATISTICS_CPP__

//
// SNP Assistant
//
// Statistical methods
//
// Authors:
//   Reedik Mägi <reedik@ut.ee>
//   Lauris Kaplinski <lauris@kaplinski.com>
//
// Copyright (C) 2002-2003 BioData, Ltd.
//

//
// This is pure C
//

#define _USE_MATH_DEFINES
#include <math.h>


//static double bdm_prob_X (int x, int n11, int n12, int n22);
double bdm_chi2 (double X, int k);
double bdm_prob_norm (double x);
double bdm_prob_chi2 (double x2, int ndf);
double bdm_signif_min (int f11, int f12, int f21, int f22);
double bdm_signif_max (int f11, int f12, int f21, int f22);
static double bdm_prob_X (int x, int n11, int n12, int n22);
double MAX(double a, double b);
int MAX(int a, int b);
