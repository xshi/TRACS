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

#ifndef TMEAS_H
#define TMEAS_H

#include "TObject.h"
#include "TDatime.h"


//CONSTANTS      ----------------------
#define NXCHAR 512


using namespace std ;       //Evita usar std::cin;

class TMeas: public TObject {

   public:
     TMeas() ;
     ~TMeas() ;
     UShort_t   NV ;	          //! (DNS) Number of voltage scans	   
     UShort_t   Nx ;	          //! (DNS) Number of different x positions 		   
     UShort_t   Ny ;	          //! (DNS) Number of different y positions 			    
     UShort_t   Nz ;	          //! (DNS) Number of different z positions 		    
     
     Double_t   At ;              //! (DNS) Time step [ns]
     Double_t   Ax ;              //  X step [mum]
     Double_t   Ay ;	          //  Y step [mum]		   
     Double_t   Az ;	          //  Z step [mum]		   

     TDatime    utc   ;           //  Time in UTC 
     UInt_t     event ;           //  DATA-START-event iterative
     
     Double_t   Itot ;            //  Current from HV supply [A]
     Double_t   x    ;            //  X [mm]
     Double_t   y    ;            //  Y [mm]
     Double_t   z    ;            //  Z [mm]

     Int_t   ix    ;              // Slice in X 
     Int_t   iy    ;	          // Slice in Y 
     Int_t   iz    ;	          // Slice in Z 

     Double_t   Temp ;            // Temperature
     Double_t   Vbias ;           //  Bias voltage [V]
     
     //I-DLTS specific
     Double_t pWidth ;            //I-DLTS pulse width
     Double_t pAmpl ;             //I-DLTS pulse amplitude
     

     char       comment[NXCHAR] ; //! (DNS) File comment
     
     Int_t      Nt       ;	  //Number of bins in scope		    
     Double_t   *volt	 ;        //[Nt]
     Double_t   *time	 ;        //[Nt]
     //One possible change: Qt is not neccessary to be declared as a pointer.
     Double_t * Qt	 ;            //[Nt]

     
     UShort_t   Setup    ;          //! 0=OldTCT , 1=eTCT (default), 2=TCT+
                                    //! 3=Hamburg, 4=IFCA , 5=TRACS, 6=IDLTS
             
     UInt_t     Ntevent ;          //! Total number of DATA-START

    ClassDef(TMeas,1)  //Edge-TCT data class
} ;


#endif
