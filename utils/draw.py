import numpy as np
import sys
from numpngw import write_png

def make_image(types,ctags,fibs,conts,cont_param):
	'''
	Draws result of simulation

	Args:
		types: np.ndarray, array with type of each cell
		ctasg: np.ndarray, array with id of cell in each point 
		conts: np.ndarray, array with cell contacts
		fibs: np.ndarray, array with fibers
		cont_param: int, 0-not to and 1-to show contacts on image
		

	Returns:
		img: np.ndarray, array with image 
	'''

	CMs_ind = np.where(types==1)[0]+1

	table=ctags
	img = np.ones(table.shape+(3,), dtype=np.uint8)*255

	edges = np.zeros(table.shape, dtype = np.uint8)
	vert = table[:,1:] != table[:,:-1]
	hor  = table[1:,:] != table[:-1,:]
	edges[:,1:]   = vert
	edges[:,:-1] += vert
	edges[1:,:]  += hor
	edges[:-1,:] += hor
	edges = np.uint8(edges != 0)

	is_CM = np.vectorize(lambda x: x in CMs_ind)

	null = np.uint8(table == 0)
	CMs = is_CM(table).astype('uint8')
	FBs = np.ones_like(table) - null - CMs

	f = 0;
	cont_edges = edges
	if cont_param != 0:		#if not 0, show contacts/attachments at least on the edge of the cell
		f = 0.5
	if cont_param == 2:		#if 2, show all of the attachment sites, even those under the cell
		cont_edges = 1 

	img[:,:,0] = conts*255
	img[:,:,1] = edges*255
	#write_png("./imgs/conts.png", img)

	Fibs = fibs * (1-CMs) * (1-FBs) * 255 * 0.1

	img[:,:,0] = CMs*(1-0.5*edges)*255 + FBs*(1-0.5*edges)*255 + f*conts*cont_edges*CMs*255 + f*conts*cont_edges*FBs*255 + Fibs
	img[:,:,1] = FBs*(1-0.5*edges)*255 + 2*f*conts*cont_edges*CMs*255 + f*conts*cont_edges*FBs*255 + Fibs
	img[:,:,2] = FBs*(1-0.5*edges)*255 - f*conts*cont_edges*FBs*255 + Fibs

	areas = np.zeros(4)
	conv = (2.5/1000)**2 #to mm
	areas = (np.sum(CMs)*conv,np.sum(FBs)*conv)

	#write_png("./imgs/exampleof.png", img)

	return img



