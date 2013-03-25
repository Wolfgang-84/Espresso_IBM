#include "vvolume.h"

double *CentVV = NULL;
int setvo = 1;
int vescnum = 0;
double *VVol = NULL;
double VVolo[200] = {0.0};

void SetCentVV() {
    int i;
	CentVV = (double *) malloc(4*vescnum*sizeof(double));
	for(i=0; i<4*vescnum; i++) {
		CentVV[i]=0.0;
	}
}

void SetVVol() {
	int i;
	VVol = (double *) malloc(vescnum*sizeof(double));
	for(i=0; i<vescnum; i++) {
		VVol[i]=0.0;
	}
}