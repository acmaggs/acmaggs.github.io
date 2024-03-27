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

#ifndef _Param_h
#define _Param_h
#include <math.h>
class Param{
 private:
  double _t;
  double _e1,_e2;
  int _L;
  int _plus;
  int _minus;
  int _neutral;
  int _nrecord;
  int _nmc;
  int _avoid;
  double _lambda;
  int _every;
  double _globalstiff;
  int _preset_e;
  double _net;
  double _fieldstep;
 public:
  int nmc() { return _nmc;}
  void  nmc(int n){_nmc=n;} 

  int avoid() { return _avoid;}
  void avoid(int n) {_avoid=n;}

  int nrecord(){return _nrecord;}
  void nrecord(int n){_nrecord=n;}

  int preset_e(){return _preset_e;}
  void preset_e(int n){_preset_e=n;}

  double t(){return _t;};
  void t(double a){ _t=a; _fieldstep=2.*sqrt(_t);};

  double globalstiff(){return _globalstiff;};
  void globalstiff(double a){ _globalstiff=a;};

  double net(){return _net;};
  void net(double a){ _net=a;};

  double e1(){ return _e1;}
  double e2() {return _e2;}

  int L() {return _L;}
  void L(int n) { _L=n;}

  double fieldstep() {return _fieldstep;}

  double lambda(){return _lambda;}
  void lambda(double x){_lambda=x;;}
  int plus() {return _plus;}
  int minus() {return _minus;}
  void minus(int i) { _minus=i;}
  void plus(int i) { _plus=i;}
  void neutral(int i) { _neutral=i;}
  int neutral() {return _neutral;}
  int np() {return _plus+_minus +_neutral;}
  double pi(){return 4*atan2(1.,1.);}
  int every() {return _every;}
  void every(int i){ _every=i;}
  Param();
  ~Param();
};
#endif
