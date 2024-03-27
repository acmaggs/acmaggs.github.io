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
#include <math.h>
#include <unistd.h>
extern Param p;
char *s[3]={"field","particle","global"};

inline double dot(double * a, double *b){
  return a[0]*b[0] + a[1]*b[1] +a[2]*b[2];
}

void
Net::mc(int n){
  // Perform the three possible MC moves with probabilites controlled by f1
  // and f2 which should be between 0 and 1.
  double f1=.5;
  double f2=.95;
  int i;
  double r;

  for(i=0;i<n;i++){
    r= RANDOM ;
    if (r < f1 ) {
      field_try();
      tries[0]++;
    }
    else if(  r > f1  && r <f2 ){
      particle_try();
      tries[1]++;
    }
    else{
      tries[2]++;
      global_try();
    }
  }
}

void
Net::field_try(){
  int acc=0;
  //Choose a plaquette and try to update circulation
  //A pure rotational motion does not change the sum
  //of the fields so we do not have to bother about Lnk::ebar
  //or Lnk::esum
  Fce *fp = fc + int(RANDOM * 3 *L3);
  double inc= ( RANDOM -.5)*p.fieldstep();
  double e1,e2;
  e1= fp->energy();
  fp->increment(inc);
  e2=fp->energy();
  double Delta_e=e2-e1;
  if(Delta_e<0) {acc=1;}
  else if( RANDOM < exp ( -Delta_e/p.t()) ){ acc=1; }
  if (acc==1){
    running_energy += Delta_e;
    accept[0]++;    
  }
  else{
    fp->refuse();
  }
}

void
Net::particle_try(){
  Ptl *pp;
  int acc=0;
  double e1,e2;
  int rn= int(RANDOM *np);
  pp = pt + rn; //Choose a particle to move
  Vtx * vp= pp->vx; //Find the site occupied by particle.
  int direc= int( RANDOM *6);
  double sgn= (direc<3) ? 1.:-1.; // are we moving in the positive or negative direction
  Vtx *v2 = vp->vx[direc]; // new site we shall try to move to
  if(v2->pt !=0) return;  // stop if site occupied

  Lnk *lp =  vp->lk[direc]; //modified link
  e1= lp->energy();
  int q=pp->q();

  lp->fieldinc(-q*sgn); //source is minus the current
  lp->globalshift(-q*sgn); 
  e2=lp->energy();

  double Delta_e=e2-e1;
  if(Delta_e < 0) { acc=1; }
  else if ( RANDOM < exp(-(Delta_e)/p.t())){acc=1; }
  if(acc==1){
    running_energy += Delta_e ; 
    pp->vx = v2;  // New vertex of the particle
    vp->pt=0;  // old site is now empty
    v2->pt = pp; //new site is now occupied.
    accept[1]++;
  }
  else{
    lp->refuse();
    lp->refuse_globalshift();
  }
}

void 
Net::global_try(){
  int acc=0;
  double e1,e2,e3,e4;
  double step[3];
  int i;
  double inc=3*sqrt( p.t())/sqrt(double(L3));
  double Delta_e;

  for(i =0; i<3;i++){step[i]=inc*( RANDOM -.5); }
  Delta_e= 2*dot(Lnk::esum,step)+L3*(2*dot(Lnk::ebar,step)+dot(step,step));
  if (Delta_e<0){acc=1;}
  else if ( RANDOM  < exp(- Delta_e/p.t()  ) ) { acc=1;  }  
  if(acc==1){
    lk[0].update_ebar(step);    
    accept[2]++;
    running_energy += Delta_e;
  }
  else {
    // This time we are up to date.
  }
}


void
Net::check(){
  cout<<"Checking "<<endl;
  int i,j,q;
  Lnk *lp;
  Vtx *vp;
  for(i=0;i<L3;i++){
    vp = vx + i;
    q =  (vp->pt==0 ? 0 : vp->pt->q());
    double sum=0;
    for(j=0;j<3;j++){
      sum += vp->lk[j]->field();
    }
    for(j=3;j<6;j++){
      sum -= vp->lk[j]->field();
    }
    //     cout<<"i\t"<<i<<"\t q "; 
    //cout<<q<<"\t "<< vp->i<<" "<<vp->j<<" "<<vp->k<<"\tsum\t "<<sum;
    // cout<< " q-sum \t"<<q-sum<<endl;
    assert( fabs(q -sum < 1.e-10));
  }
}


void 
Net::stats(){
  cout<<endl;
  for(int i=0;i<3;i++){
    if(tries[i] >0){
      cout<<s[i]<<" "<< 100*float(accept[i])/float(tries[i])<<"\t";
    }
  }
  cout<<endl;
}
