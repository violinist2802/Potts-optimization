from utils.parse_config import parse_config
from libc.stdlib cimport malloc, free
import numpy as np

from Cython.Compiler import Options

Options.docstrings = True

cdef extern from "libcpmfem.h":
	int cpmfem(int NCX, int NCY, 
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
	double distanceF, int NRINC,
	char* typ,
	int* cont_m,
	int* fibr,
	int* ctag_m)
	
def py_cpmfem(int NCX, int NCY, PART, double VOXSIZE, double sizeX, double sizeY, scenario, NRINC, **kwargs):
	'''
	Simulates VCT model
	Args:
		NCX: int, number of cells by x-axis
		NCY: int, number of cells by y-axis
		PART: float, percentage of fibroblasts
		NRINC: int, number of simulation steps
		VOXSIZE: double, domain size in mm
		sizeX: double, horizontal simulation area size mm
		sizeY: double, vertical simulation area size mm
		scenario: str, one of 4('monolayer_on_fiber' , 'monolayer_without_fiber' , 'single_on_fiber' , 'single_without_fiber')
	Returns:
		(np.arrays for usage in draw.py)
		types, ctags, fibers, contacts
	'''
	cdef char* typ = <char*>malloc(NCX*NCY*sizeof(char))
	sizeMarginX = 0.1
	sizeMarginY = 0.1
	NVX = int(round((sizeX+2*sizeMarginX)/VOXSIZE))
	NVY = int(round((sizeY+2*sizeMarginY)/VOXSIZE))
	cdef int NV = NVX*NVY
	cdef int* cont_m = <int*>malloc(NVX*NVY*sizeof(int))
	cdef int* fibr = <int*>malloc(NVX*NVY*sizeof(int))
	cdef int* ctag_m = <int*>malloc(NVX*NVY*sizeof(int))

	cfg = parse_config('./utils/config.yaml', scenario)
	for arg, value in kwargs.items():
		cfg[arg]=value
	cpmfem(NCX, NCY, PART, VOXSIZE, NVX, NVY, cfg['GN_CM'], cfg['GN_FB'], cfg['TARGETVOLUME_CM'], cfg['TARGETVOLUME_FB'], cfg['DETACH_CM'], cfg['DETACH_FB'], cfg['INELASTICITY_FB'], cfg['INELASTICITY_CM'],  cfg['JCMMD'], cfg['JFBMD'], cfg['JCMCM'], cfg['JFBFB'], cfg['JFBCM'], cfg['UNLEASH_CM'], cfg['UNLEASH_FB'], cfg['LMAX_CM'], cfg['LMAX_FB'], cfg['MAX_FOCALS_CM'], cfg['MAX_FOCALS_FB'], cfg['shifts'], cfg['distanceF'], NRINC, typ, cont_m, fibr, ctag_m)
	types=[]
	ctags=[]
	fibers=[]
	contacts=[]
	for i in range(NVY):
		ctags.append([])
		fibers.append([])
		contacts.append([])
	for i in range(NVX):
		for j in range(NVY):
			k=i+j*NVX
			ctags[i].append(ctag_m[k])
			fibers[i].append(fibr[k])
			contacts[i].append(cont_m[k])
	for i in range(NCX*NCY):
		types.append(int(typ[i]))
	
	types=np.array(types)
	ctags=np.array(ctags)
	fibers=np.array(fibers)
	contacts=np.array(contacts)		
	free(typ)
	free(cont_m)
	free(fibr)
	free(ctag_m)
	return types, ctags, fibers, contacts

