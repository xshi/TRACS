/*
 * @ Copyright 2014-2017 CERN and Instituto de Fisica de Cantabria - Universidad de Cantabria. All rigths not expressly granted are reserved [tracs.ssd@cern.ch]
 * This file is part of TRACS.
 *
 * TRACS is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public License as published by the Free Software Foundation,
 * either version 3 of the Licence.
 *
 * TRACS is distributed in the hope that it will be useful , but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * See the GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License along with TRACS. If not, see <http://www.gnu.org/licenses/>
 */

#include "TMeasHeader.h"
#include <string.h>  //strcpy

ClassImp(TMeasHeader)


TMeasHeader::TMeasHeader( ) {
     Nav=0;
     NV=Nx=Ny=Nz=0 ;
     At=Ax=Ay=Az=Cend=tann=Etann=iann=Fluence=0. ;
     Temp=Position=0. ;
     comment=TString("");
     vVbias = 0.0 ;
     Polarity=0;
     Setup = 1;
     Hum = Illum = Power = Gain = 0.0;	 
     Lambda =1060.0 ;    
     Freq = 200 ;	  

}

TMeasHeader::TMeasHeader( Int_t NV) {
     Nx=Ny=Nz=0 ;
     At=Ax=Ay=Az=Cend=tann=Etann=iann=Fluence=0. ;
     Temp=0. ;
     comment=TString("");
     vVbias.ResizeTo(NV) ; //http://root.cern.ch/root/roottalk/roottalk09/1231.html
     vIleak.ResizeTo(NV) ; //http://root.cern.ch/root/roottalk/roottalk09/1231.html
     vTemp.ResizeTo(NV) ;  //http://root.cern.ch/root/roottalk/roottalk09/1231.html
     Polarity=0;
     Setup = 1;
     Hum = Illum = Power = Gain = 0.0;	 
     Lambda =1060.0 ;    
     Freq = 200 ;	  
}

TMeasHeader::~TMeasHeader()
{
  //Emptiness
  vVbias.Clear() ;
  
}

