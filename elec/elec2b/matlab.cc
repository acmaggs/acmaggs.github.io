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
#include "Anal.h"
#include <fstream.h>
#include "Param.h"
#include <assert.h>

int Anal::inx(int n1,int n2,int n3,int nq){
  /*cout<< "n "<<endl;
  cout<< n1<<" ";
  cout<< n2<<" ";
  cout<< n3<<" ";
  cout<< nq<<" ";*/
  int tr;
  assert( n1>0 && n2>0 && n3>0 && nq>0);
  assert( n1 <= dims[0]);  assert( n2 <= dims[1]);
  assert( n3 <= dims[2]);  assert( nq <= dims[3]);
  tr=(n1-1) + (n2-1)*dims[0] + (n3-1)*dims[0]*dims[1] +  (nq-1)*dims[0]*dims[1]*dims[2];
  assert( tr <= dims[0] *dims[1] *dims[2]*dims[3] );
  return tr;
}

Anal::Anal(int nrecords, int nq){
  dims[0]=nrecords;
  if(nrecords==0) {
    cout<<"*********************Analyse is OFF********************"<<endl;
    return;
  }
  dims[1]=nq; // number of wavevectors
  dims[2]=3; //Density, field , charge //Field not used anymore 
  dims[3]=2; //Real-Imaginary
  cout<<"Dims "<<dims[0]<<" ";
  cout<< dims[1]<<" ";   
  cout<< dims[2]<<" ";   
  cout<< dims[3]<<" "; 
  cout<<endl;
  pm=matOpen("net.mat","w");
  p=mxCreateNumericArray(ndim,  dims, mxSINGLE_CLASS, mxREAL);
  assert(p);
  mxSetName(p,"sq");
  matlab=   (float *) mxGetPr(p);
  ncall=0;
}

Anal::~Anal(){
  cout<<"Anal finished\t"<<dims[0]<<endl;
  if(dims[0] !=0){
    cout<<"Ncall = "<<ncall<<endl;
    anal<<";"<<endl;
    memcpy( (char *) mxGetPr(p), (char * ) matlab,dims[0]*dims[1]*dims[2]*dims[3]*sizeof(float));
    matPutArray(pm,p);
    anal.close();
  }
}

void
Anal::analyse(double f[2], double rho[2],double ch[2], int nq ){
   if(nq==0) ncall++;
  nq++; //Fortran labeling
  if(dims[0] ==0 ) return;

  matlab[ inx(ncall, nq, 1, 1) ] = float(f[0]);
  matlab[ inx(ncall, nq, 1, 2) ] = float(f[1]);

  matlab[ inx(ncall, nq, 2, 1) ] = float(rho[0]);
  matlab[ inx(ncall, nq, 2, 2) ] = float(rho[1]);

  matlab[ inx(ncall, nq, 3, 1) ] = float(ch[0]);
  matlab[ inx(ncall, nq, 3, 2) ] = float(ch[1]);
}

