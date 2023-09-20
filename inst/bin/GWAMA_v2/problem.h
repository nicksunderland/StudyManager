#pragma once

class problem
{
public:
	int markersAll;		//count of all markers
	int markersOK;		//count of ok markers
	int problemStrand;  //strand problem
	int wrongAlleles;	//wrong alleles problem
	int problemMulti;	//multiple occurances of SNP
	int problemEffect;	//effect is missing or problematic

	problem(void);
	~problem(void);
};
