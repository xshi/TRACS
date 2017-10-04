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

#ifndef TMEASHEADER_H
#define TMEASHEADER_H

#include <TObject.h>
#include <TString.h>
#include <TDatime.h>
#include <TVectorD.h>

using namespace std ;       //Evita usar std::cin;

class TMeasHeader: public TObject {

   public:
     TMeasHeader() ;
     TMeasHeader( Int_t NV ) ;
     ~TMeasHeader();
     
     UShort_t   NV ;	          //Number of voltage scans	   
     UShort_t   Nx ;	          //Number of different x positions 		   
     UShort_t   Ny ;	          //Number of different y positions 			    
     UShort_t   Nz ;	          //Number of different z positions 		    
     Short_t    Polarity ;	  //Polarity of the FILE! 		    
     
     Double_t   At ;              //Time step [ns]
     Double_t   Ax ;              //X step [mum]
     Double_t   Ay ;              //Y step [mum]                 
     Double_t   Az ;              //Z step [mum]                 
     Double_t   Nt ;              //Number of steps in time
     
     TString    comment ;         //File comment
     Double_t   Temp ;            //Temperature
     Double_t   Hum ;             //Humidity
     Double_t   Lambda ;          //Lambda
     Double_t   Power ;           //Pulse Power
     Double_t   Gain ;            //Amp. Gain
     Double_t   Freq ;            //Frequency
     Int_t      Nav ;             //Ner of averages
     Double_t   t0 ;              //Scope t0
     Double_t   Position ;        //Position of the sensor in mm (mostly for TCT standard)
  
     Double_t   Cend ;            //Capacitance used for RC deconvolution of the file
     UShort_t   Setup ;           //0=OldTCT , 1=eTCT (default), 2=TCT+
                                  //3=Hamburg, 4=IFCA , 5=TRACS, 6=TPA
     
     Int_t      Illum ;           //Illumination 1 = top , 0=side, -1 = bottom

     TDatime date0 ;              //Start date
     TDatime datef ;              //End date
     
     Double_t   iann ;            //annealing step
     Double_t   tann ;            //Specific annealing time at TempAnn
     Double_t   Etann ;           //Accumulated annealing time at 60C for this sample
     Double_t   TempAnn ;         //Annealing temperature
     Double_t   Fluence ;         //Fluence

     TVectorD   vVbias ;          //Vector with bias voltages
     TVectorD   vIleak ;          //Vector with leakage current
     TVectorD   vTemp ;           //Vector with temperature readings
     
     //Methods
     UShort_t GetNV() { return NV ;}
     UShort_t GetNx() { return Nx ;}
     UShort_t GetNy() { return Ny ;}
     UShort_t GetNz() { return Nz ;}
     UShort_t GetNav(){ return Nav ;}

     Double_t GetAt() { return At ;}
     Double_t GetAx() { return Ax ;}
     Double_t GetAy() { return Ay ;}
     Double_t GetAz() { return Az ;}
     Double_t GetGain()    { return Gain ;}
     Double_t Getiann()    { return iann ;   }  
     Double_t Gettann()    { return tann ;   }  
     Double_t GetEtann()   { return Etann ;  } 
     Double_t GetTempAnn() { return TempAnn ;} 
     Double_t GetFluence()  { return Fluence ;} 
     Double_t GetFrequency()    { return Freq;}
     Double_t GetFilePolarity() { return Polarity ;}

     Double_t GetTemperature()    {return Temp ;} 	    
     TString  GetComment() { return comment; } 

    ClassDef(TMeasHeader,3)  //Edge-TCT data header class
} ;

#endif
