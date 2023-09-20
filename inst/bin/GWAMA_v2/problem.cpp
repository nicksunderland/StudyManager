#include "problem.h"

problem::problem(void)
{
	markersAll = 0;		//count of all markers
	markersOK = 0;		//count of ok markers
	problemStrand = 0;  //strand problem
	wrongAlleles = 0;	//wrong alleles problem
	problemMulti = 0;	//multiple occurances of SNP
	problemEffect = 0;
}

problem::~problem(void)
{
}
