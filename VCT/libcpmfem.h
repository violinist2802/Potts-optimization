
#ifndef EXAMPLES_H
#define EXAMPLES_H

int cpmfem(
	int NCX, int NCY, 
	double PART,
	double VOXSIZE,
	int NVX, int NVY,
	double GN_CM,
	double GN_FB,
	double TARGETVOLUME_CM,
	double TARGETVOLUME_FB,
	double DETACH_CM,
	double DETACH_FB,
	double INELASTICITY_FB,
	double INELASTICITY_CM,
	double JCMMD,
	double JFBMD,
	double JCMCM,
	double JFBFB,
	double JFBCM,
	double UNLEASH_CM,
	double UNLEASH_FB,
	double LMAX_CM,
	double LMAX_FB,
	double MAX_FOCALS_CM,
	double MAX_FOCALS_FB,
	int shifts,
	double distanceF,
	int NRINC,
	char* typ,
	int* cont_m,
	int* fibr,
	int* ctag_m
);

#endif