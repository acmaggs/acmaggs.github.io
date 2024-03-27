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
#include "fstream.h"

void
Net::print(){
  ofstream out("net.out");
  int i,j;
  out<<"L "<<L<<endl;
  out<<"L3 "<<L3<<endl;
  out<<"3*L3 "<<3*L3<<endl;

  out<<"\n Vtx-Vtx\n";
  for(i=0;i<L3;i++) {
    out<<i<<"\t\t";
    for(j=0;j<6;j++){
      out<< vx[i].vx[j] - vx<<"\t";
    }
    out<<endl;
  }

  out<<"\n Vtx-Lnk\n";
  for(i=0;i<L3;i++) {
    out<<i<<"\t\t";
    for(j=0;j<6;j++){
      out<< vx[i].lk[j] - lk<<"\t";
    }
    out<<endl;
  }

  out<<"\n Vtx-Ptl\n"<<endl;
  for(i=0;i<L3;i++) {
    if(vx[i].pt != 0){
     out<<"i "<<i<<"\t\t";
    out<< vx[i].pt - pt<< "\t"<<vx[i].pt->q()<< endl;
    } 
  }
  out<<endl;

  out<<"\n Lnk-Vtx\n";
  for(i=0;i<3*L3;i++){
    out<<lk[i].vx[0] - vx<<" "<<lk[i].vx[1]-vx<<"\t"<<lk[i].dir<<endl;
  }

  out<<"\n Field\n";
  for(i=0;i<3*L3;i++){
    out<<lk[i].field()<<"\t";
    if( (i+1) % L ==0) out<<endl;
    if( (i+1) % (L*L) ==0) out<<endl;
  }


 out<<"\nPtl->Vtx\n"<<endl;
 for(i=0;i<np;i++){
   out<<i<<"\t"<< pt[i].vx-vx<<"\t"<< pt[i].q()<<endl;
 }

 out<<"\nFce->Lnk\n"<<endl;
 for(i=0;i<3*L3;i++){
   out<<i<<"\t\t";
   for(j=0;j<4;j++){
     out<< fc[i].lk[j]-lk<<"\t";
   }
   out<<endl;
 }
}
