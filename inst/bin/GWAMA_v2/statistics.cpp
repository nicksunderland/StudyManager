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

//#include "macros.h"
//#include "compat.h"
#include "statistics.h"


double 
MAX (double a, double b)
{
if (a>b) return a;
return b;
}

int 
MAX (int a, int b)
{
if (a>b) return a;
return b;
}





static double bdm_prob_X (int x, int n11, int n12, int n22);


double
bdm_t2x2_chi (double n11, double n12, double n21, double n22)
{
	double nt = n11 + n12 + n21 + n22;
	if (nt <= 1e-9) return -1.0;
	double rnt = 1.0 / nt;
	double e11 = rnt * (n11 + n12) * (n11 + n21);
	double e12 = rnt * (n11 + n12) * (n12 + n22);
	double e21 = rnt * (n21 + n22) * (n11 + n21);
	double e22 = rnt * (n21 + n22) * (n12 + n22);
	double x11 = (e11 != 0.0) ? ((n11 - e11) * (n11 - e11)) / e11 : 0;
	double x12 = (e12 != 0.0) ? ((n12 - e12) * (n12 - e12)) / e12 : 0;
	double x21 = (e21 != 0.0) ? ((n21 - e21) * (n21 - e21)) / e21 : 0;
	double x22 = (e22 != 0.0) ? ((n22 - e22) * (n22 - e22)) / e22 : 0;
	double xt = x11 + x12 + x21 + x22;
	if ((e11 < 5.0) || (e12 < 5.0) || (e21 < 5.0) || (e22 < 5.0)) return -1.0;
	return xt;
}

double
bdm_t2x3_chi (double n11, double n12, double n13, double n21, double n22, double n23)
{
    double N[3][4], E[2][3];
    double chi2;
    int r, c, lowc;
    // 2x3 Matrix
    N[0][0] = n11;
    N[0][1] = n12;
    N[0][2] = n13;
    N[1][0] = n21;
    N[1][1] = n22;
    N[1][2] = n23;
    // Rows
    for (r = 0; r < 2; r++) N[r][3] = N[r][0] + N[r][1] + N[r][2];
    // Columns
    for (c = 0; c < 3; c++) N[2][c] = N[0][c] + N[1][c];
    // Total
    N[2][3] = N[0][3] + N[1][3];
    // Expected
    lowc = 0;
    for (r = 0; r < 2; r++) {
        for (c = 0; c < 3; c++) {
            E[r][c] = (N[r][3] * N[2][c]) / N[2][3];
            if (E[r][c] <= 0.0) return -1.0;
            if (E[r][c] < 5.0) lowc += 1;
        }
    }
    if (lowc > 1) return -1.0;
    chi2 = 0;
    for (r = 0; r < 2; r++) {
        for (c = 0; c < 3; c++) {
            chi2 += (E[r][c] - N[r][c]) * (E[r][c] - N[r][c]) / E[r][c];
        }
    }
    return chi2;
}

double
bdm_t2x2_add (int n11, int n12, int n13, int n21, int n22, int n23)
{
	double x11 = n13 + n12 / 2.0;
	double x12 = n11 + n12 / 2.0;
	double x21 = n23 + n22 / 2.0;
	double x22 = n21 + n22 / 2.0;

	double x1_ = x11 + x12;
	double x2_ = x21 + x22;
	double x_1 = x11 + x21;
	double x_2 = x12 + x22;
	double x__ = x_1 + x_2;
	if (x__ != 0.0) {
		double expect_x11 = x_1 * x1_ / x__;
		double expect_x12 = x_2 * x1_ / x__;
		double expect_x21 = x_1 * x2_ / x__;
		double expect_x22 = x_2 * x2_ / x__;
		if ((expect_x11 >= 5.0) && (expect_x12 >= 5.0) && (expect_x21 >= 5.0) && (expect_x22 >= 5.0)) {
			double vahe = x11 * x22 - x21 * x12;
			double korrutis = x1_ * x2_ * x_1 * x_2;
			return (vahe * vahe * x__) / korrutis;
		}
	}
	return -1.0;
}

double
bdm_t2x2_hz (int n11, int n12, int n13, int n21, int n22, int n23)
{
	double x11 = n11 + n13;
	double x12 = n12;
	double x21 = n21 + n23;
	double x22 = n22;

	double x1_ = x11 + x12;
	double x2_ = x21 + x22;
	double x_1 = x11 + x21;
	double x_2 = x12 + x22;
	double x__ = x_1 + x_2;
	if (x__ != 0.0) {
		double expect_x11 = x_1 * x1_ / x__;
		double expect_x12 = x_2 * x1_ / x__;
		double expect_x21 = x_1 * x2_ / x__;
		double expect_x22 = x_2 * x2_ / x__;
		if ((expect_x11 >= 5.0) && (expect_x12 >= 5.0) && (expect_x21 >= 5.0) && (expect_x22 >= 5.0)) {
			double vahe = x11 * x22 - x21 * x12;
			double korrutis = x1_ * x2_ * x_1 * x_2;
			return (vahe * vahe * x__) / korrutis;
		}
	}
	return -1.0;
}

double
bdm_t2x2_nr (int n11, int n12, int n13, int n21, int n22, int n23)
{
	double x11 = n11;
	double x12 = n12 + n13;
	double x21 = n21;
	double x22 = n22 + n23;

	double x1_ = x11 + x12;
	double x2_ = x21 + x22;
	double x_1 = x11 + x21;
	double x_2 = x12 + x22;
	double x__ = x_1 + x_2;
	if (x__ != 0.0) {
		double expect_x11 = x_1 * x1_ / x__;
		double expect_x12 = x_2 * x1_ / x__;
		double expect_x21 = x_1 * x2_ / x__;
		double expect_x22 = x_2 * x2_ / x__;
		if ((expect_x11 >= 5.0) && (expect_x12 >= 5.0) && (expect_x21 >= 5.0) && (expect_x22 >= 5.0)) {
			double vahe = x11 * x22 - x21 * x12;
			double korrutis = x1_ * x2_ * x_1 * x_2;
			return (vahe * vahe * x__) / korrutis;
		}
	}
	return -1.0;
}

double
bdm_t2x2_dr (int n11, int n12, int n13, int n21, int n22, int n23)
{
	double x11 = n13;
	double x12 = n11 + n12;
	double x21 = n23;
	double x22 = n21 + n22;

	double x1_ = x11 + x12;
	double x2_ = x21 + x22;
	double x_1 = x11 + x21;
	double x_2 = x12 + x22;
	double x__ = x_1 + x_2;
	if (x__ != 0.0) {
		double expect_x11 = x_1 * x1_ / x__;
		double expect_x12 = x_2 * x1_ / x__;
		double expect_x21 = x_1 * x2_ / x__;
		double expect_x22 = x_2 * x2_ / x__;
		if ((expect_x11 >= 5.0) && (expect_x12 >= 5.0) && (expect_x21 >= 5.0) && (expect_x22 >= 5.0)) {
			double vahe = x11 * x22 - x21 * x12;
			double korrutis = x1_ * x2_ * x_1 * x_2;
			return (vahe * vahe * x__) / korrutis;
		}
	}
	return -1.0;
}

double
bdm_chi2 (double X, int k)
{
	if ((X <= 0.0) || (k < 1)) return 1.0;

	if (k == 1) {
		return 2 * bdm_prob_norm (sqrt (X));
	}

	if (k <= 30) {
		if (k & 1) {
			/* Odd */
			/* Calculate R */
			int k_1_2 = k >> 2;
			double DEN = 1.0;
			double R = 0.0;
			int r;
			for (r = 1; r <= (k_1_2); r++) {
				DEN *= (2 * r - 1);
				R += pow (X, r - 0.5) / DEN;
			}
			/* Calculate Q */
			return 2 * bdm_prob_norm (sqrt (X)) + R * sqrt (2 / M_PI) * exp (-X / 2.0);
		} else {
			/* Even */
			int k_2_2 = (k >> 2) - 1;
			double DEN = 1.0;
			double R = 0.0;
			int r;
			for (r = 1; r < k_2_2; r++) {
				DEN *= (2 * r);
				R += pow (X, r) / DEN;
			}
			/* Calculate Q */
			return exp (-X / 2.0) * (1.0 + R);
		}
	}

	if (X >= 150.0) {
		return 0.0;
	}

	double Z;

	if (X == (k - 1.0)) {
		Z = -(1 / 3.0 + 0.08 / k) / sqrt (2.0 * k - 2);
	} else {
		double d;
		d = X - k + 2.0 / 3.0 - 0.08 / k;
		Z = d * sqrt ((k - 1) * log ((k - 1) / X) + X - (k - 1)) / fabs (X - (k - 1));
	}

	if (Z < 0.0) {
//		return 1 - bdm_prob_norm (Z);
		return  bdm_prob_norm (Z);
	}

	return bdm_prob_norm (Z);
}

double
bdm_prob_norm_our (double x)
{
	double p;
	double xa = fabs (x);
	if (xa > 12.0) {
		p = 0;
	} else {
		p = 1 / (1 + 0.33267 * xa);
		p = 0.39894228 * exp (-0.5 * (xa * xa)) * p * (0.4361836 + p * (0.937298 * p - 0.1201676));
	}

	return (x >= 0.0) ? p : (1.0 - p);
}

double
bdm_prob_norm (double x)
{
        long double p;
        long double xa = fabs (x);
        if (xa > 120.0) {
                p = 0;
        } else {
                p = 1 / (1 + 0.33267 * xa);
                p = 0.39894228 * exp (-0.5 * (xa * xa)) * p * (0.4361836 + p * (0.937298 * p - 
0.1201676));
        }

        return (x >= 0.0) ? p : (1.0 - p);
}










double
bdm_prob_chi2 (double x2, int ndf)
{
	if ((x2 > 0.0) && (ndf > 0)) {
		if (ndf == 1) {
			double x = sqrt (x2);
			return 2 * bdm_prob_norm (x);
		} else {
			if (x2 > 169) {
				double ch = 2.0 / (9.0 * ndf);
				double x = (exp (log (x2 / ndf) / 3.0) - 1.0 + ch) / sqrt (ch);
				return bdm_prob_norm (x);
			} else {
				if (ndf == 2) {
					return exp (-0.5 * x2);
				} else {
					int N1 = (ndf - 1) / 2;
					int N2 = (ndf - 2) / 2;
					if (N1 == N2) {
						double sum = 1.0;
						double re = 0.5 * x2;
						double ch = 1.0;
						int i;
						for (i = 1; i <= N2; i++) {
							ch = ch * re / i;
							sum += ch;
						}
						return exp (-re) * sum;
					} else {
						double ch = sqrt (x2);
						double z = 0.39894228 * exp (-0.5 * x2);
						double p = bdm_prob_norm (ch);
						if (ndf == 3) {
							return 2 * (p + z * ch);
						} else {
							double chp = ch;
							double re = 1.0;
							int i;
							for (i = 2; i <= N1; i++) {
								re += 2;
								chp *= (x2 / re);
								ch += chp;
							}
							return 2 * (p + z * ch);
						}
					}
				}
			}
		}
	}

	return -1.0;
}

double
bdm_signif_min (int f11, int f12, int f21, int f22)
{
	static double *L = (double *) 0;
	static double *N = (double *) 0;
	static int size = 0;

	int len = 2 * (f11 + f12 + f21 + f22);
	if (len > size) {
		// Need to allocate more
		if (L) delete[] L;
		if (N) delete[] N;
		size = 2 * size;
		size = MAX (size, len);
		L = new double[size];
		N = new double[size];
	}

	int f_1 = f11 + f21;
	int f_2 = f12 + f22;
	int f1_ = f11 + f12;
	int f2_ = f21 + f22;
	int f__ = f_1 + f_2;

	double prob_f11 = 0.0;

	while ((f22 >= 0.0) && (f11 >= 0.0)) {
		double prob_f11_liidetav = 1.0;
		int g = 0;
		int k;
		for (k = 0; k < f1_; k++) L[g++] = k + 1;
		for (k = 0; k < f2_; k++) L[g++] = k + 1;
		for (k = 0; k < f_2; k++) L[g++] = k + 1;
		for (k = 0; k < f_1; k++) L[g++] = k + 1;

		g = 0;
		for (k = 0; k < f__; k++) N[g++] = k + 1;
		for (k = 0; k < f11; k++) N[g++] = k + 1;
		for (k = 0; k < f12; k++) N[g++] = k + 1;
		for (k = 0; k < f21; k++) N[g++] = k + 1;
		for (k = 0; k < f22; k++) N[g++] = k + 1;

		len = g;
		for (g = 0; g < len; g++) {
			double jagatis = L[g] / N[g];
			prob_f11_liidetav *= jagatis;
		}

		prob_f11 += prob_f11_liidetav;

		f11 -= 1;
		f12 += 1;
		f21 += 1;
		f22 -= 1;
	}

	return prob_f11;
}

double
bdm_signif_max (int f11, int f12, int f21, int f22)
{
	static double *L = (double *) 0;
	static double *N = (double *) 0;
	static int size = 0;

	int len = 2 * (f11 + f12 + f21 + f22);
	if (len > size) {
		// Need to allocate more
		if (L) delete[] L;
		if (N) delete[] N;
		size = 2 * size;
		size = MAX (size, len);
		L = new double[size];
		N = new double[size];
	}

	int f_1 = f11 + f21;
	int f_2 = f12 + f22;
	int f1_ = f11 + f12;
	int f2_ = f21 + f22;
	int f__ = f_1 + f_2;

	double prob_f11 = 0.0;

	while ((f12 >= 0.0) && (f21 >= 0.0)) {
		double prob_f11_liidetav = 1.0;
		int g = 0;
		int k;
		for (k = 0; k < f1_; k++) L[g++] = k + 1;
		for (k = 0; k < f2_; k++) L[g++] = k + 1;
		for (k = 0; k < f_2; k++) L[g++] = k + 1;
		for (k = 0; k < f_1; k++) L[g++] = k + 1;

		g = 0;
		for (k = 0; k < f__; k++) N[g++] = k + 1;
		for (k = 0; k < f11; k++) N[g++] = k + 1;
		for (k = 0; k < f12; k++) N[g++] = k + 1;
		for (k = 0; k < f21; k++) N[g++] = k + 1;
		for (k = 0; k < f22; k++) N[g++] = k + 1;

		len = g;
		for (g = 0; g < len; g++) {
			double jagatis = L[g] / N[g];
			prob_f11_liidetav *= jagatis;
		}

		prob_f11 += prob_f11_liidetav;

		f11 += 1;
		f12 -= 1;
		f21 -= 1;
		f22 += 1;
	}

	return prob_f11;
}

#ifdef BDM_DEBUG
#include <stdio.h>
#endif

static double
bdm_prob_X (int x, int n11, int n12, int n22)
{
	int N = n11 + n12 + n22;
	int N_2 = 2 * N;
	int c1 = n11 + (n12 - x) / 2;
	int c2 = n22 + (n12 - x) / 2;
	int c3 = 2 * n11 + n12;
	int c4 = N_2 - c3;

	double p = 1.0;
	for (int k = 1; k <= N_2; k++) {
		int e = -1;
		if (k <= N) e += 1;
		if (k <= c1) e -= 1;
		if (k <= x) e -= 1;
		if (k <= c2) e -= 1;
		if (k <= c3) e += 1;
		if (k <= c4) e += 1;
		if (e > 0) {
			for (int i = 0; i < e; i++) p *= k;
		} else if (e < 0) {
			for (int i = 0; i > e; i--) p /= k;
		}
		if (k <= x) p *= 2;
	}

#ifdef BDM_DEBUG
	double q = p;
	p = 1.0;
	for (int k = 1; k <= N_2; k++) {
		if (k <= N) p *= k;
		if (k <= c1) p /= k;
		if (k <= x) p /= k;
		if (k <= x) p *= 2;
		if (k <= c2) p /= k;
		p /= k;
		if (k <= c3) p *= k;
		if (k <= c4) p *= k;
	}
	if (fabs (1.0 - p/q) > 1e-9) {
		printf ("ERERERERE\n");
	}
#endif

	return p;
}



