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
#include "TMeas.h"

ClassImp(TMeas)

TMeas::TMeas( ) {
     NV=Nx=Ny=Nz=0;
     At=Ax=Ay=Az=0. ;	    
     event=Ntevent=0 ;	    
     Itot=x=y=z=Vbias=0. ;	    
     Temp=999.0;
     ix=iy=iz=0 ;
     
     Setup = 1 ;
     Nt=0 ;
     volt=0 ;
     time=0 ;
     Qt=0;

     pWidth = pAmpl = 0. ;
}

TMeas::~TMeas()
{
  delete [] volt ;
  delete [] time ;
  delete [] Qt ;
}

//____________________________________________________________________________


