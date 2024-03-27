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
#include "Anal.h"

extern Anal * a;

void 
Net::sq(double field[2],double particle[2],double ch[2],const int k[3]){
  double q0=2*p.pi()/L;
  double qx=k[0]*q0;
  double qy=k[1]*q0;
  double qz=k[2]*q0;

  field[0]=running_energy;
  field[1]=running_energy;
  particle[0]=0;
  particle[1]=0;
  ch[0]=0;
  ch[1]=0;

  for(int l=0;l<p.np();l++){
    Ptl *pp = pt + l;
    Vtx *vp = pp->vx;
    int q = pp->q();
    int i=vp->i;
    int j=vp->j;
    int k=vp->k;
    double inx= qx * i +qy * j + qz * k;
    double c=  cos(inx);
    double s=  sin(inx);
    particle[0] += c;
    particle[1] += s;
    ch[0] += q * c;
    ch[1] += q * s;
  }
}

void Net::record(){
  static int firstcall=1;
  double field[2],particle[2],ch[2];
  const static int wv[NMODE][3]= {
    {0,0,0}, {1,0,0}, {1,1,0} , {1,1,1}, 
      {2,0,0}, {2,1,0},  {2,1,1}, {1,2,1},
	{2,2,0},{0,2,2}, {1,2,2}, {2,2,1},
	  {2,2,2}, {3,0,0}, {3,1,0}, {3,1,1},
	    {3,2,1}, {3,2,2},{3,3,3}, {4,0,0},
	      {4,0,4},{4,1,1}, {1,1,4},{3,2,0}, 
		{5,0,5}, {5,5,5} 
  } ;


 if(firstcall){
    firstcall=0;
    ofstream modes("qq.mat");
    for(int i=0;i<NMODE;i++){
      modes<< wv[i][0]* wv[i][0]  +wv[i][1]* wv[i][1] + wv[i][2]* wv[i][2] <<" ";
    }
    modes<<endl;
    modes.close();
  }

  for(int i=0;i<NMODE;i++){
    sq(field,particle,ch,wv[i]);
    a->analyse(field,particle,ch,i);
  }

}
