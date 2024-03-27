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

#ifndef _Simul_h
#define _Simul_h
#include <stdlib.h>
#include "Param.h"
#include "MersenneTwister.h"
#include <assert.h>


const int NMODE=26;
extern Param p;
extern MTRand r_gen;

#define RANDOM ( r_gen() )
//#define RANDOM ( drand48() )
#include <iostream.h>

class Vtx;
class Lnk;
class Fce;
class Net;
class Ptl;

class Ptl{
  friend class Net;
 private:
  Vtx *vx;
  double eps;
  int charge;
 public:
  Ptl();
  void epsilon(double e){ eps=e;}
  double epsilon(){ return eps; }
  void q(int ch) { charge=ch;}
  int q(){ return charge;}
};

class Vtx{
 friend class Net;
 private:
 int i,j,k;
  Vtx  *vx[6];
  Lnk  *lk[6];
  Ptl  *pt;
 public:
  Vtx();
};

class Lnk{
  friend class Net;
 private:
  Vtx *vx[2];
  //Fce *fc[4];
  static double oldebar[3];
  static double oldesum[3];
  double field_;
  double oldfield;
  short dir;
  static double ebar[3];
  static double esum[3];
 public:
  void globalshift(double d){
    oldesum[dir] = esum[dir];
    esum[dir] += d;
  }
  void refuse_globalshift(){
    esum[dir] = oldesum[dir];
  };
  void update_ebar(double d[3]){
    for(int i=0;i<3;i++){
      oldebar[i]=ebar[i];
      ebar[i] += d[i];
    }
  }

  void fieldinc(double d){
    oldfield=field_;
    field_ +=d;
  }

  double energy(){
    double f = field_+ebar[dir];    
    return f*f;
  }

  void refuse() {
    field_=oldfield;
  }
  double field(){
    return field_;
  }

  void static refuse_ebar(){
    for(int i=0;i<3;i++){
      ebar[i]=oldebar[i];
    }
  }
  Lnk();
};

class Fce{
  friend class Net;
 private:
  //Vtx *vx[4];
  Lnk *lk[4];
 public:
  double energy(){
    double e=0;
    for (int i=0;i<4;i++){
      e += lk[i]->energy();
    }
    return e;
}

  void increment(double d){
    lk[0]->fieldinc(d);
    lk[1]->fieldinc(d);
    lk[2]->fieldinc(-d);
    lk[3]->fieldinc(-d);
  }
  void refuse(){
    for (int i=0;i<4;i++){
      lk[i]->refuse();
    }
}

Fce();
};


class Net{
 private:
  double running_energy;
  double full_energy;
  int L;
  int L3;

  int np;
  Vtx *vx;
  Lnk *lk;
  Fce *fc;
  Ptl *pt;
  void initvtx();
  void initfce();
  int tries[3];
  int accept[3];
  void place_particles();
 public:
  void stats();
  void print();
  void sq(double *, double*, double *, const int *);
  void record();
  double energy();
  void check();
  void field_try();
  void particle_try();
  void global_try();
  void mc(int);
  Net(int,int);
  ~Net();
};

#endif
