/*    Copyright (C) <2003>  Centre National de la Recherche Scientifique, CNRS.
    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program; if not, write to the Free Software
    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

Program written by A.C. Maggs, tony@turner.pct.espci.fr
Laboratoire de physico-chimie theorique.
CNRS UMR 7083
10 rue Vauquelin 75005 Paris.
*/

#include "simul.h"

int right(int i,int L) {return (i+1)%L;  };
int left(int i,int L){    return (i-1+L)%L;  };
int inx(int i,int j,int k,int L){ return i+L*j+L*L*k;}

void Net::initfce(){
  for(int i=0;i<L;i++){
    for(int j=0;j<L;j++){
      for(int k=0;k<L;k++){
	int p0 = inx(i,j,k,L);
	//X-Y plane*****************************
	int p1= inx(i,j,k,L);
	//	fc[p0].vx[0]=vx+p1;
	fc[p0].lk[0]= lk + p0; //first X oriented link
	fc[p0].lk[3] = lk + p0 +L3; //second Y oriented link

	p1= inx( right(i,L), j, k, L);
	//fc[p0].vx[1]=vx +p1;
	fc[p0].lk[1] = lk +p1 + L3;  //first Y oriented link

	p1= inx( right(i,L), right(j,L), k, L);
	//fc[p0].vx[2]=vx +p1;

	p1= inx( i, right(j,L), k, L);
	//fc[p0].vx[3]=vx +p1;
	fc[p0].lk[2] = lk +p1;// second X oriented link

	//Y-Z Plane********************
	p1 = inx(i,j,k,L);
	//fc[p0].vx[0] = vx +p1;
	fc[p0+L3].lk[0]= lk + p0 +L3; //first Y oriented link
	fc[p0+L3].lk[3] = lk +p0 + 2*L3; //second Z oriented link

	p1= inx( i, right(j,L), k, L);
	//fc[p0].vx[1]=vx +p1;
	fc[p0+L3].lk[1] = lk +p1 + 2*L3;  //first Z oriented link

	p1= inx( i, right(j,L), right(k,L), L);
	//fc[p0].vx[2]=vx +p1;

	p1= inx( i, j, right(k,L), L);
	//fc[p0].vx[3]=vx +p1;
	fc[p0+L3].lk[2] = lk +p1 +L3;// second Y oriented link

	//Z-X plane*************************
	p1 = inx(i,j,k,L);
	//fc[p0].vx[0] = vx +p1;
	fc[p0+2*L3].lk[0]= lk + p0 +2*L3; //first Z oriented link
	fc[p0+2*L3].lk[3] = lk +p0; //second X oriented link

	p1= inx( i, j, right(k,L), L);
	//fc[p0].vx[1]=vx +p1;
	fc[p0+2*L3].lk[1] = lk +p1;  //first X oriented link

	p1= inx( right(i,L), right(j,L), k, L);
	//fc[p0].vx[2]=vx +p1;

	p1= inx( right(i,L), j, k, L);
	//fc[p0].vx[3]=vx +p1;
	fc[p0+2*L3].lk[2] = lk +p1 +2*L3;// second Z oriented link
      }
    }
  }
}


void Net::initvtx(){
  // Vertices are on a simple cubic lattice numbered 0...(L-1)
  // in each axis. There are L3 such sites.

  // There are 3*L3 links. These are stocked in the following manner
  // X oriented links in (0... L3-1)
  // Y oriented kinks in (L3 ... 2*L3-1)
  // Z oriented links in (2*L3... 3*L3-1)
  for(int i=0;i<L;i++){
    for(int j=0;j<L;j++){
      for(int k=0;k<L;k++){
	int p0 = inx(i,j,k,L);
	vx[p0].i=i;
	vx[p0].j=j;
	vx[p0].k=k;
	// Three sites in positive directions 
	//X
	int p1= inx( right(i,L),j,k,L);
	vx[p0].vx[0] =  vx + p1;  
	vx[p0].lk[0] =  lk + p0;

	lk[p0].vx[0] =  vx + p0;
	lk[p0].vx[1] =  vx + p1;
	lk[p0].dir=0;
	//Y
	p1= inx( i,right(j,L),k,L);
	vx[p0].vx[1] =  vx + p1;  
	vx[p0].lk[1] =  lk + p0 + L3;
	lk[p0+L3].vx[0] =  vx + p0;
	lk[p0+L3].vx[1] =  vx + p1;
	lk[p0+L3].dir=1;
	//Z
	p1= inx( i,j,right(k,L),L);
	vx[p0].vx[2] =  vx + p1;  
	vx[p0].lk[2] =  lk + p0 + 2 * L3;
	lk[p0+2*L3].vx[0] =  vx + p0;
	lk[p0+2*L3].vx[1] =  vx + p1;
	lk[p0+2*L3].dir=2;

	//Three sites in negative directions 
	//X
	p1= inx( left(i,L),j,k,L);
	vx[p0].vx[3] =  vx + p1;  
	vx[p0].lk[3] =  lk + p1;
	//Y
	p1= inx( i,left(j,L),k,L);
	vx[p0].vx[4] =  vx + p1;  
	vx[p0].lk[4] =  lk + p1 + L3;
	//Z
	p1= inx( i,j,left(k,L),L);
	vx[p0].vx[5] =  vx + p1;  
	vx[p0].lk[5] =  lk + p1 + 2 * L3;
      }
    }
  }
 }
