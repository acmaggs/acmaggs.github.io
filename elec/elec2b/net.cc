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

void 
Net::place_particles(){
  int nf=0;
  Lnk *lp;
  //place particles by choosing a link. Put one particle
  // at each end of the link. Random initial conditions.
  // Update the field between the two links

  while (nf<np){
    int r= int(RANDOM *3 *L3); //Chose a random link
    //cout<<"R "<<r<<endl;
    lp = lk + int(r);
    if (lp->vx[0]->pt==0 && lp->vx[1]->pt==0){//If neigbouring pair of sites empty
      //      cout<<"V1 "<< lp->vx[0] - vx<<"\t V2\t"<<lp->vx[1]-vx<<endl;;
      //cout<<"Nf  "<<nf<<endl;
      int sn= (RANDOM - .5)> 0 ? 1:-1 ; //charge of first particle
      pt[nf].q(sn); //set charge
      lp->vx[0]->pt = pt + nf ; //place particle of first vertex
      pt[nf].vx= lp->vx[0];
      nf++;

      pt[nf].q(-sn); //set charge
      lp->vx[1]->pt = pt + nf; //place particle on second vertex
      pt[nf].vx= lp->vx[1];
      nf++;

      lp->field_ = sn; //set field for connecting link
    }
  }
  cout<<"Particles placed "<<nf<<" "<<np<<endl;
}


Ptl::Ptl(){
  charge=0;
  eps=0;
  vx=0;
}

Net::Net(int lin,int n){
  L=lin;
  L3=L*L*L;
  np=n;
  vx= new Vtx[L3];
  lk= new Lnk[3*L3];
  fc= new Fce[3*L3];
  pt= new Ptl[np];
  initvtx();
  initfce();
  place_particles();
  for(int i=0;i<3;i++){
    tries[i]=0;
    accept[i]=0;
  }
  cout<<"Set energy in constructor "<<endl;
  running_energy=0;
  energy();
  running_energy=full_energy;
}

Net::~Net()
{
  delete [] vx;
  delete [] lk;
  delete [] fc;
  delete [] pt;
}


Vtx::Vtx(){
  int l;
  i=j=k=-1;
  pt=0;
  for(l=0;l<6;l++){
    lk[l]=0;
    vx[l]=0;
  }
}

double Lnk::ebar[3]={0,0,0};
double Lnk::oldebar[3]={0,0,0};
double Lnk::esum[3]={0,0,0};
double Lnk::oldesum[3]={0,0,0};

Lnk::Lnk(){
  field_=0;
  oldfield=0;
  dir=-1;
  vx[0]=vx[1]=0;
  //  for(int i=0;i<4;i++) fc[i]=0;
}

Fce::Fce(){
  for (int i=0;i<4;i++){
    lk[i]=0;
    //    vx[i]=0;
  }
}
