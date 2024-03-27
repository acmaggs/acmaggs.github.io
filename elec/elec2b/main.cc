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

Param p;
MTRand r_gen;
Anal *a;

void setparms(){
  p.L(10);
  p.minus(30);
  p.plus(30);
  p.neutral(0);
  p.nrecord(1000);
  int nstep= p.L()*p.L()*p.L()*4;
  p.nmc(nstep);
  p.t(.25);
}

void loop(Net *n){
  n->mc((p.nrecord()*p.nmc())/15);
  cout<< n->energy()<<endl;
  for(int i=0;i<p.nrecord();i++){
    if( (50*i) % p.nrecord() ==0 ) {
      cout<< 100.*double(i)/double(p.nrecord() ) <<"% of simulation\t";
      n->stats();
      n->energy();
      n->check();
    }
    n->mc( p.nmc() );
    n->record();
  }
  cout<< n->energy()<<endl;
}

int main(){
  setparms();
  Net n( p.L() , p.np() );
  a= new Anal(p.nrecord(),NMODE);
  n.energy();
  n.check();
  loop(&n);
  n.print();
  delete a;
  return 0;
}
