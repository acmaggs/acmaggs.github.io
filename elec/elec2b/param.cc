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
#include <iostream.h>
#include <fstream.h>
#include <math.h>
Param::Param(){
  _L=5;
  _plus=3;
  _minus=3;
  _neutral=0;
  _nrecord=4000;
  _t=10.;
  _e1=1;
  _e2=1;
  _nmc=100;
  _avoid=1;
  _every=2;
  _globalstiff=.0;
  _preset_e=0;
  _net=1;
  _fieldstep=2*sqrt(_t);
}


Param::~Param(){
  ofstream out("data.out");
  out<<"L "<<_L<<endl;;
  out<<"plus " <<_plus<<endl;;
  out<<"minus "<<_minus<<endl;;
  out<<"neutral "<<_neutral<<endl; 
  out<<"nrecord "<<_nrecord<<endl;
  out<<"t "<<_t<<endl;
  out<<"e1 "<<_e1<<endl;
  out<<"e2 "<<_e2<<endl;
  out<<"nmc "<<_nmc <<endl;
  out<<"avoid "<<_avoid <<endl;
  out<<"globalstiff "<<_globalstiff <<endl;
  out<<"preset_e "<<_preset_e <<endl;
  out<<"net "<<_net <<endl;
}
