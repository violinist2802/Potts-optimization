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
	double distanceF, int NRINC,
	char* typ,
	int* cont_m,
	int* fibr,
	int* ctag_m,
	float* PART_matrix)
	
def py_cpmfem(int NCX, int NCY, PART, double VOXSIZE, double sizeX, double sizeY, scenario, NRINC, PART_matrix, **kwargs):
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
	
	cdef float* matrix = <float*>malloc(NVX*NVY*sizeof(float))
	if type(PART_matrix)==float:
		for i in range(NVX*NVY):
			matrix[i]=PART_matrix
	if type(PART_matrix)==type(np.array([[0,1],[1,0]])):
		for i in range(NVY):
			for j in range(NVX):
				k=j+NVX*i
				matrix[k]=PART_matrix[i][j]
	cfg = parse_config('./utils/config.yaml', scenario)
	for arg, value in kwargs.items():
		cfg[arg]=value
	cpmfem(NCX, NCY, PART, VOXSIZE, NVX, NVY, 1, cfg['GN_CM'], cfg['GN_FB'], cfg['TARGETVOLUME_CM'], cfg['TARGETVOLUME_FB'], cfg['DETACH_CM'], cfg['DETACH_FB'], cfg['INELASTICITY_FB'], cfg['INELASTICITY_CM'],  cfg['JCMMD'], cfg['JFBMD'], cfg['JCMCM'], cfg['JFBFB'], cfg['JFBCM'], cfg.get('JZ_CMCM', 0.0), cfg.get('JZ_FBFB', 0.0), cfg.get('JZ_CMFB', 0.0), cfg['UNLEASH_CM'], cfg['UNLEASH_FB'], cfg['LMAX_CM'], cfg['LMAX_FB'], cfg['MAX_FOCALS_CM'], cfg['MAX_FOCALS_FB'], cfg['shifts'], cfg['distanceF'], NRINC, typ, cont_m, fibr, ctag_m, matrix)
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
			ctags[j].append(ctag_m[k])
			fibers[j].append(fibr[k])
			contacts[j].append(cont_m[k])
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
	free(matrix)
	return types, ctags, fibers, contacts

# Multi-layer interface
def py_cpmfem_layers(layers, double VOXSIZE, double sizeX, double sizeY, scenario, NRINC, **kwargs):
    """
    layers: list[dict] where each dict defines a layer configuration.
      Supported keys per layer:
        - type: 'basic' | 'prob' | 'ready'
        - NCX, NCY (for basic/prob), PART or PART_matrix
        - ready_ctag (NVY x NVX), ready_type (NCX*NCY) [not yet used]
    Returns per-layer arrays: types[L][NCX*NCY], ctags[L][NVY,NVX], fibers[L][NVY,NVX], contacts[L][NVY,NVX]
    """
    cdef int NL = len(layers)
    sizeMarginX = 0.1
    sizeMarginY = 0.1
    NVX = int(round((sizeX+2*sizeMarginX)/VOXSIZE))
    NVY = int(round((sizeY+2*sizeMarginY)/VOXSIZE))
    cdef int NV = NVX*NVY

    # Assume uniform NCX, NCY across layers unless overridden per layer; use first
    int_NCX = layers[0].get('NCX', 0)
    int_NCY = layers[0].get('NCY', 0)
    # Fallback to kwargs if provided
    if int_NCX == 0:
        int_NCX = kwargs.get('NCX', 0)
    if int_NCY == 0:
        int_NCY = kwargs.get('NCY', 0)

    cdef char* typ = <char*>malloc(int_NCX*int_NCY*NL*sizeof(char))
    cdef int* cont_m = <int*>malloc(NV*NL*sizeof(int))
    cdef int* fibr = <int*>malloc(NV*NL*sizeof(int))
    cdef int* ctag_m = <int*>malloc(NV*NL*sizeof(int))
    cdef float* matrix3d = <float*>malloc(NV*NL*sizeof(float))

    # Build per-layer probability matrix and handle ready layers
    ready_layers = []
    for li, lconf in enumerate(layers):
        ltype = lconf.get('type', 'basic')
        if ltype == 'prob':
            pm = lconf.get('PART_matrix', None)
            if pm is None:
                raise ValueError('Layer with type=prob requires PART_matrix')
            for i in range(NVY):
                for j in range(NVX):
                    k = j + NVX*i
                    matrix3d[li*NV + k] = pm[i][j]
        elif ltype == 'ready':
            # Store ready layer data for later use
            ready_ctag = lconf.get('ready_ctag', None)
            ready_type = lconf.get('ready_type', None)
            if ready_ctag is None or ready_type is None:
                raise ValueError('Layer with type=ready requires ready_ctag and ready_type')
            ready_layers.append((li, ready_ctag, ready_type))
            # Set matrix to 0 for ready layers (they use predefined config)
            for k in range(NV):
                matrix3d[li*NV + k] = 0.0
        else:  # basic
            part = lconf.get('PART', 0.0)
            for k in range(NV):
                matrix3d[li*NV + k] = part

    cfg = parse_config('./utils/config.yaml', scenario)
    for arg, value in kwargs.items():
        cfg[arg]=value

    # Use PART from first layer for init (C still expects PART scalar but uses matrix for placement)
    PART_scalar = layers[0].get('PART', 0.0)

    cpmfem(int_NCX, int_NCY, PART_scalar, VOXSIZE, NVX, NVY, NL, cfg['GN_CM'], cfg['GN_FB'], cfg['TARGETVOLUME_CM'], cfg['TARGETVOLUME_FB'], cfg['DETACH_CM'], cfg['DETACH_FB'], cfg['INELASTICITY_FB'], cfg['INELASTICITY_CM'],  cfg['JCMMD'], cfg['JFBMD'], cfg['JCMCM'], cfg['JFBFB'], cfg['JFBCM'], cfg.get('JZ_CMCM', 0.0), cfg.get('JZ_FBFB', 0.0), cfg.get('JZ_CMFB', 0.0), cfg['UNLEASH_CM'], cfg['UNLEASH_FB'], cfg['LMAX_CM'], cfg['LMAX_FB'], cfg['MAX_FOCALS_CM'], cfg['MAX_FOCALS_FB'], cfg['shifts'], cfg['distanceF'], NRINC, typ, cont_m, fibr, ctag_m, matrix3d)

    # Pack outputs per layer with proper numpy array shapes
    types_out = []
    ctags_out = np.zeros((NL, NVY, NVX), dtype=np.int32)
    fibers_out = np.zeros((NL, NVY, NVX), dtype=np.int32)
    contacts_out = np.zeros((NL, NVY, NVX), dtype=np.int32)
    
    for li in range(NL):
        # types per layer
        t_layer = []
        for i in range(int_NCX*int_NCY):
            t_layer.append(int(typ[li*int_NCX*int_NCY + i]))
        types_out.append(np.array(t_layer))

        # grids with shape (NL, NVY, NVX)
        for i in range(NVX):
            for j in range(NVY):
                k = i + j*NVX
                base = li*NV
                ctags_out[li, j, i] = ctag_m[base + k]
                fibers_out[li, j, i] = fibr[base + k]
                contacts_out[li, j, i] = cont_m[base + k]

    free(typ)
    free(cont_m)
    free(fibr)
    free(ctag_m)
    free(matrix3d)
    return types_out, ctags_out, fibers_out, contacts_out

