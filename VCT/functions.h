// file: functions.h
// #include "def.h"
#include "structures.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>

// init.c
int* 		alloc_attach(int NRc);
VOX*		init_voxels(int NVX, int NVY);
int 		init_cells(VOX* pv, int * types, BOX* pb, int NCX, int NCY, double PART, int shifts, double TARGETVOLUME_FB, double VOXSIZE, int NVX, int NVY);
FIBERS* 	set_fibers(double distanceF, double VOXSIZE, int NVX, int NVY);


// cellmoves.c
double 		CPM_moves(VOX* pv, short * CCAlabels,BOX* pb, FIBERS* pf, 
CM* CMs, int* attached, int* csize, double MAX_FOCALS_CM, double MAX_FOCALS_FB, double TARGETVOLUME_CM, double TARGETVOLUME_FB, double INELASTICITY_CM, double INELASTICITY_FB, double LMAX_CM, double LMAX_FB, double GN_CM, double GN_FB, double UNLEASH_CM, double UNLEASH_FB, double DETACH_CM, double DETACH_FB, double VOXSIZE, int NVX, int NVY, double JCMCM, double JCMMD, double JFBFB, double JFBMD,  double JFBCM);
double 		CH_moves(VOX* pv, CM* CMs, double k, double VOXSIZE, int NVX, int NVY);
BOOL 		splitcheckCCR(VOX* pv, short * CCAlabels, BOX* pb,  
int* csize, int xt, int ttag, int NVX, int NVY);

// CM.c
CM* 		allocCM(int NRc);
BOX*		allocBOX(int NRc);
void 		findCM(VOX* pv, CM* CMs, int NRc, int NVX, int NVY);

// CPM_dH.c
// distribute channels
double 		calcdH_CH(VOX* pv, CM* CMs, int xt, int xs, int NVX);
double 		calcdHborder(VOX* pv, int xt, int ttag, int NVX);
double 		calcdHdist(VOX* pv, CM* CMs, int xt, int xs, int ttag, int NVX);

//calculate H
double 		calcdH(VOX* pv, FIBERS* pf, CM* CMs, int* csize, int xt, int xs, int pick, int ttag, int stag, double TARGETVOLUME_CM, double TARGETVOLUME_FB, double INELASTICITY_CM, double INELASTICITY_FB, double LMAX_CM, double LMAX_FB, double GN_CM, double GN_FB, double UNLEASH_CM, double UNLEASH_FB, double DETACH_CM, double DETACH_FB, double VOXSIZE, int NVX, double JCMCM, double JCMMD, double JFBFB, double JFBMD,  double JFBCM);
double 		calcdHcontact(VOX* pv, int xt, int xs, int ttag, int stag, int NVX, double JCMCM, double JCMMD, double JFBFB, double JFBMD,  double JFBCM);
double 		contactenergy(int tag1, int tag2, int type1, int type2, double JCMCM, double JCMMD, double JFBFB, double JFBMD, double JFBCM);
double 		scaffoldenergy(int tag, int Q);
double 		calcdHvol(int* csize, int ttag, int stag, int ttype, int stype, double TARGETVOLUME_CM, double TARGETVOLUME_FB, double INELASTICITY_CM, double INELASTICITY_FB);
double 		calcdHprotrude(VOX* pv, CM* CMs, int xt, int xs, int ttag, int stag, int Qt, int Qs, double LMAX_CM, double LMAX_FB, double GN_CM, double GN_FB, double UNLEASH_CM, double UNLEASH_FB, double DETACH_CM, double DETACH_FB, int NVX);
double 		calcdHsyncytium(VOX* pv, CM* CMs, int xt, int xs, int ttag, int stag, int NVX);
double 		calcdHnuclei(VOX* pv, CM* CMs, int xt, int ttag, int stag, double DETACH_CM, double DETACH_FB, double VOXSIZE, int NVX);
//double      printdH(VOX* pv, FIBERS* pf, CM* CMs, int* csize, int xt, int xs, int pick, int ttag, int stag)

double 		findphi(CM* CMs, int xt, int tag);
double 		dist(CM* CMs, int xt, int tag, int NVX);

double		sqr(double x);

// write.c
void   		write_increment(int increment);
void 		write_cells(VOX* pv, int increment, int NVX, int NVY);
void 		write_types(int* types, int NRc);
void 		write_contacts(VOX* pv, int increment, int NVX, int NVY);
void 		write_fibers(FIBERS* pf, int NVX, int NVY);
void 		read_cells(VOX* pv, int* types, int NRc, char filename_ctag[40], char filename_cont[40], char filename_types[40]);

// mylib.c
void 		myitostr(int n, char s[]);
void 		myreverse(char s[]);
unsigned 	mystrlen(const char *s);

// mt.c
void 		mt_init();
unsigned long mt_random();


