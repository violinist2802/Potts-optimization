#include "functions.h"
#include <sys/time.h>
#include <unistd.h>
#include <getopt.h>
#include <math.h>
#include "libcpmfem.h"



int cpmfem(
	int NCX, int NCY, 
	double PART,
	double VOXSIZE,
	int NVX, int NVY,
    int NL,
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
	double JZ_CMCM,
	double JZ_FBFB,
	double JZ_CMFB,
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
	int* ctag_m,
	float* PART_matrix
)
{

	struct timeval tv;
    time_t time;

	int d;
	VOX *pv;
	FIBERS *pf;
	int NRc,c,v;
	int * types;
	int *csize;
	int incr, startincr;
	double acceptance, acceptance_phi;
	CONT = 0;
    CONT_INHIB = 1;
	
	if(!silence){
		printf("SEED = %d\n",SEED);
		printf("Sample size = %d x %d\n",NCX,NCY);
		printf("\n");
		printf("GN_CM = %.2f\n",GN_CM*SCALE);
		printf("TARGETVOLUME_CM = %.2f\n",TARGETVOLUME_CM*1000*VOXSIZE*VOXSIZE);
		printf("INELASTICITY_CM = %.2f\n",INELASTICITY_CM/(SCALE*SCALE*SCALE*SCALE));
		printf("DETACH_CM = %.2f\n",DETACH_CM*SCALE);
		printf("\n");
		printf("GN_FB = %.2f\n",GN_FB*SCALE);
		printf("TARGETVOLUME_FB = %.2f\n",TARGETVOLUME_FB*1000*VOXSIZE*VOXSIZE);
		printf("INELASTICITY_FB = %.2f\n",INELASTICITY_FB/(SCALE*SCALE*SCALE*SCALE));
		printf("DETACH_FB = %.2f\n",DETACH_FB*SCALE);
		printf("\n");
		printf("JMDMD = %.3f\n",JMDMD/VOXSIZE);
		printf("JCMMD = %.3f\n",JCMMD/VOXSIZE);
		printf("JFBMD = %.3f\n",JFBMD/VOXSIZE);
		printf("\n");
		printf("JCMCM = %.3f\n",JCMCM/VOXSIZE);
		printf("JFBFB = %.3f\n",JFBFB/VOXSIZE);
		printf("JFBCM = %.3f\n",JFBCM/VOXSIZE);
		printf("\n");
		printf("UNLEASH_CM = %.2f\n",UNLEASH_CM*SCALE);
		printf("UNLEASH_FB = %.2f\n",UNLEASH_FB*SCALE);
		printf("\n");
		printf("LMAX_CM = %.2f px\n",LMAX_CM);
		printf("LMAX_FB = %.2f px\n",LMAX_FB);
		printf("\n");
		printf("MAX_FOCALS_CM = %.2f\n",MAX_FOCALS_CM);
		printf("MAX_FOCALS_FB = %.2f\n",MAX_FOCALS_FB);
	}

	/// INITIALIZE ///
	srand(SEED); mt_init();

	int layer;
	int offsetNV;
	int offsetTypes;

	gettimeofday(&tv, NULL);
	time = tv.tv_sec;

	for(layer=0; layer<NL; layer++){
		pv = init_voxels(NVX, NVY);
		pf = set_fibers(distanceF, VOXSIZE, NVX, NVY);
		BOX * pb = allocBOX(NCX*NCY+1);

		write_fibers(pf, NVX, NVY);

        startincr = 0;
        types = calloc((NCX*NCY+1), sizeof(int));
        // For now, use regular init_cells - ready layer support would need additional parameters
        NRc = init_cells(pv,types,pb,NCX,NCY, PART, shifts, TARGETVOLUME_FB, VOXSIZE, NVX, NVY, PART_matrix + layer*NV);
        write_cells(pv,0, NVX, NVY);
		csize = calloc(NRc, sizeof(int)); for(c=0;c<NRc;c++) {csize[c]=0;}
		for(v=0;v<NV;v++) {if(pv[v].ctag) {csize[pv[v].ctag-1]++;}}

		CM* CMs = allocCM(NRc);
		int* attached = alloc_attach(NRc);

		short * CCAlabels = malloc(NV * sizeof(short));

		// START SIMULATION ///
		for(incr=startincr; incr<NRINC; incr++)
		{
			if (incr % STEP_PRINT == 0){
				if(!silence)
					printf("\nL%d START INCREMENT %d",layer,incr);
				write_cells(pv,incr, NVX, NVY);
				write_contacts(pv,incr, NVX, NVY);
			}

			findCM(pv,CMs,NRc, NVX, NVY);
			acceptance = CPM_moves(pv,CCAlabels,pb,pf,CMs,attached,csize,MAX_FOCALS_CM,MAX_FOCALS_FB, TARGETVOLUME_CM, TARGETVOLUME_FB, INELASTICITY_CM, INELASTICITY_FB, LMAX_CM, LMAX_FB, GN_CM, GN_FB, UNLEASH_CM, UNLEASH_FB, DETACH_CM, DETACH_FB, VOXSIZE, NVX, NVY, JCMCM, JCMMD, JFBFB, JFBMD, JFBCM, JZ_CMCM, JZ_FBFB, JZ_CMFB);

			if (incr % STEP_PRINT == 0 && !silence){
				printf("\nL%d Acceptance rate %.4f",layer,acceptance);
			}
		}

		// channel distribution
		findCM(pv,CMs,NRc, NVX, NVY);
		for(incr=startincr; incr<NRINC_CH; incr++)
		{
			if (incr % 100 == 0){
				if(!silence)
					printf("\nL%d START CHANNEL DISTRIBUTION %d",layer,incr);
				write_contacts(pv,incr+1, NVX, NVY);
			}

			acceptance = CH_moves(pv, CMs, 0.5 + 0.5*incr/NRINC, VOXSIZE, NVX, NVY);

			if (incr % 100 == 0 && !silence){
				printf("\nL%d Acceptance rate %.4f",layer,acceptance);
			}
		}

		// write outputs for this layer
		int i,j,k;
		offsetNV = layer * NV;
		for(i=0; i<NVX; i++) {
			for (j=0; j<NVY; j++) {
				k = i + j * NVX;
				ctag_m[offsetNV + k]=pv[k].ctag;
				cont_m[offsetNV + k]=pv[k].contact;
				fibr[offsetNV + k]=(char)pf[k].Q;
			}
		}
		int l;
		offsetTypes = layer * (NCX*NCY);
		for(l=0; l<NRc; l++) {
		  typ[offsetTypes + l] = types[1+l];
		}

		free(pv);
		free(pf);
		free(CCAlabels);
		free(attached);
		free(CMs);
		free(csize);
		free(types);
		free(pb);
	}

	if(!silence){
		gettimeofday(&tv, NULL);
		printf("\nSIMULATION FINISHED!\n");
		printf("Took %lds\n", tv.tv_sec - time);
	}

	return 0;
}

