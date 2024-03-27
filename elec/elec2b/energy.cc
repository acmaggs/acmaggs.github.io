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

double
Net::energy(){//The energy is in the links
  double e=0;
  int i;
  for(i=0;i<3*L3;i++){
    e += lk[i].energy();
  }
  //  e/=1.;

  cout<<"Esum1\t";
  for(i=0;i<3;i++){
    cout<<Lnk::esum[i]<<"\t";
    Lnk::esum[i]=0;
  }
  cout<<endl;

  for(i=0;i<L3;i++){
    Lnk::esum[0] += lk[i].field(); 
    Lnk::esum[1] += lk[i+L3].field(); 
    Lnk::esum[2] += lk[i+2*L3].field(); 
  }
  cout<<"Esum2\t";
  for(i=0;i<3;i++){
    cout<<Lnk::esum[i]<<"\t";
  }
  cout<<endl;

  cout<<"Ebar\t";
  for(i=0;i<3;i++){
    cout<<Lnk::ebar[i]<<"\t";
  }
  cout<<endl;

  cout<<"Total fields\t";
  for(i=0;i<3;i++){
    cout<<Lnk::ebar[i]*L3 + Lnk::esum[i] <<"\t";
  }
  cout<<endl;

  cout<<endl;
  full_energy=e;
  cout<<"Full, running energies "<<full_energy/L3<<"\t"<<running_energy/L3<<"\t "
      <<(full_energy-running_energy)<<endl;
  return (e);
}


