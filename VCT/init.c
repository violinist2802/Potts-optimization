// file: init.c
#include "functions.h"

////////////////////////////////////////////////////////////////////////////////
int* alloc_attach(int NRc)
{
	int i;
	int * attached;

	attached = calloc(NRc+1, sizeof(int));

	for(i=0;i<NRc+1;i++)
		attached[i]=0;
	
	return attached;
}

////////////////////////////////////////////////////////////////////////////////
VOX* init_voxels(int NVX, int NVY)
{
	VOX* pv;
	int v, vx, vy;
	int i;

   	pv = calloc(NV, sizeof(VOX));

	// set voxel information
   	for(vy=0; vy<NVY; vy++)
   	for(vx=0; vx<NVX; vx++)
   	{
   		v = vx + vy*NVX;

		//pv[v].x = vx * VOXSIZE; pv[v].y = vy * VOXSIZE;
		pv[v].ctag = 0;
		pv[v].type = 0;
		pv[v].contact = 0;
		pv[v].bond = 0;
	}
	return pv;
}

////////////////////////////////////////////////////////////////////////////////
int init_cells(VOX* pv, int * types, BOX* pb, int NCX, int NCY, double PART, int shifts, double TARGETVOLUME_FB, double VOXSIZE, int NVX, int NVY, float* PART_matrix)
{
	int v, vx, vy, i, j, ix, iy, dxi, dyi, count, k;
	int NRc;
	double prob, prob_norm;
	double r01;
	double d;
	double dx, dy,dvx,dvy; // distance to center
	int r;
	NRc = 0;

	dx = (double) (NVX - 2 * MARGINX) / (NCX);
	dy = (double) (NVY - 2 * MARGINY) / (NCY);
	r  = (int) (sqrt(STARTVOLUME)/2);
	if(dx<2*r || dy<2*r)
		printf("Too dense!");
	for (iy = 0; iy < NCY; iy++){
		for (ix = 0; ix < NCX; ix++){
			dvx = (mt_random()%((int) dx-2*r+1)) -(dx/2 - r);
			dvy = (mt_random()%((int) dy-2*r+1)) -(dy/2 - r);
			vx = MARGINX + (int) (((double) ix + 0.5) * dx + shifts*dvx);
			vy = MARGINY + (int) (((double) iy + 0.5) * dy + shifts*dvy);
			count=0;
			prob=0;
			prob_norm=0;
			for (dxi = MARGINX+ix*(int)dx+1; dxi < MARGINX+(ix+1)*(int)dx-1; dxi++){
				for(dyi = MARGINY+iy*(int)dy+1; dyi < MARGINY+(iy+1)*(int)dy-1; dyi++){
					k=dxi+NVX*dyi;
					prob=prob+PART_matrix[k];
					count=count+1;

				}
			}
			prob_norm=prob/count;
			
			NRc++;
			types[NRc] = (prob_norm<(rand()/(double)RAND_MAX) ? 1 : 2);
			pb[NRc].x1 = vx-r;
			pb[NRc].x2 = vx+r;
			pb[NRc].y1 = vy-r;
			pb[NRc].y2 = vy+r;
			for(i = -r; i<=r; i++){
				for (j = -r; j<=r; j++){
					v = vx + i + (vy + j)*NVX;
					if (v<NV){
						pv[v].ctag = NRc;
						pv[v].type = types[NRc]; 
					}
					else
						printf("Cell out of area: (%d,%d)\n",vx+i-NVX,vy+j-NVY);
				}
			}
		}	
	}

	return NRc;
}

////////////////////////////////////////////////////////////////////////////////
FIBERS* set_fibers(double distanceF, double VOXSIZE, int NVX, int NVY)
{
	FIBERS* pf;
	
	int v, vx, vy, vd, fd, kc;
	int i;
	double dx,dy;
	double k, k0;

	dx = F_ANGLE!=0    ? F_DISTANCE / sin(F_ANGLE) : 0.0;
	dy = F_ANGLE!=PI/2 ? F_DISTANCE / cos(F_ANGLE) : 0.0;

   	pf = calloc(NV, sizeof(FIBERS));

	// set voxel information
    for(v=0; v<NV; v++)
    	pf[v].Q = 0;

   	if(F_ANGLE!=PI/2){
   		fd = round(fiberD / VOXSIZE / cos(F_ANGLE));
	   	fd = fd<1 ? 1 : fd;
	   	for(vx=0; vx<NVX; vx++){
	   		k0 = fmod(vx*tan(F_ANGLE),dy);
	   		k = k0;
	   		for(vy=0; vy<=(NVY-k0)/dy; vy++){
	   			for(vd=0; vd<fd; vd++){
	   				kc = round(k)+vd;
	   				v = vx + kc*NVX;
	   				if(kc<NVY)
	   					pf[v].Q = 1;
	   			}
	   			k += dy;
	   		}
	   	}
	}

	if(F_ANGLE!=0){
		fd = round(fiberD / VOXSIZE / sin(F_ANGLE));
	   	fd = fd<1 ? 1: fd;
	   	printf("fiberdX: %d\n", fd);
	   	for(vy=0; vy<NVY; vy++){
	   		k0 = fmod(vy/tan(F_ANGLE),dx);
	   		k = k0;
	   		for(vx=0; vx<=(NVX-k0)/dx; vx++){
	   			for(vd=0; vd<fd; vd++){
	   				kc = round(k)+vd;
		   			v = kc + vy*NVX;
		   			if(kc<NVX)
		   				pf[v].Q = 1;
	   			}
	   			k += dx;
	   		}
	   	}
	}

	return pf;
}

////////////////////////////////////////////////////////////////////////////////
int init_ready_cells(VOX* pv, int * types, BOX* pb, int NCX, int NCY, int* ready_ctag, int* ready_type, double VOXSIZE, int NVX, int NVY)
{
	int v, vx, vy, i, j, ix, iy, NRc;
	int r;
	
	NRc = 0;
	r = (int) (sqrt(STARTVOLUME)/2);
	
	// Initialize all voxels as empty
	for(vy=0; vy<NVY; vy++)
		for(vx=0; vx<NVX; vx++)
		{
			v = vx + vy*NVX;
			pv[v].ctag = 0;
			pv[v].type = 0;
			pv[v].contact = 0;
			pv[v].bond = 0;
		}
	
	// Set up cells from ready configuration
	for(ix=0; ix<NCX; ix++){
		for(iy=0; iy<NCY; iy++){
			if(ready_type[ix*NCY + iy] > 0) { // if cell exists
				NRc++;
				types[NRc] = ready_type[ix*NCY + iy];
				
				// Find bounding box for this cell
				int x1=NVX, x2=0, y1=NVY, y2=0;
				for(vy=0; vy<NVY; vy++){
					for(vx=0; vx<NVX; vx++){
						if(ready_ctag[vy*NVX + vx] == NRc) {
							if(vx < x1) x1 = vx;
							if(vx > x2) x2 = vx;
							if(vy < y1) y1 = vy;
							if(vy > y2) y2 = vy;
						}
					}
				}
				pb[NRc].x1 = x1; pb[NRc].x2 = x2;
				pb[NRc].y1 = y1; pb[NRc].y2 = y2;
				
				// Set voxels
				for(vy=0; vy<NVY; vy++){
					for(vx=0; vx<NVX; vx++){
						if(ready_ctag[vy*NVX + vx] == NRc) {
							v = vx + vy*NVX;
							pv[v].ctag = NRc;
							pv[v].type = types[NRc];
						}
					}
				}
			}
		}
	}
	
	return NRc;
}