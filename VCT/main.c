#include "def.h"
#include "libcpmfem.h"
int main(int argc, char *argv[]){
	PySys_SetArgv(argc, argv);
	parse_options(argc, argv);
	double PART = 0.1;
	cpmfem( GN_CM,
	TARGETVOLUME_CM,
	TARGETVOLUME_FB,
	DETACH_CM,
	GN_FB,
	INELASTICITY_FB,
	INELASTICITY_CM,
	DETACH_FB,
	JCMMD,
	JFBMD,
	JCMCM,
	JFBFB,
	JFBCM,
	UNLEASH_CM,
	UNLEASH_FB,
	LMAX_CM,
	LMAX_FB,
	MAX_FOCALS_CM,
	MAX_FOCALS_FB,
	NCX,
	NCY,
	PART);
	return 0;
}