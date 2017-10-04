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

#ifndef TWAVEFORM_H
#define TWAVEFORM_H

#include <TObject.h>
#include <TH1.h>


#define RTIME 0.8



#include "TMeas.h"

using namespace std ;       //Evita usar std::cin;


class TWaveform : public TObject {
		
  public:
	  	  
       TWaveform() ;
       TWaveform( TMeas *em ) ;
       TWaveform( int nel , double *tin , double *vin , double Bias  ) ;
	  ~TWaveform() ;
	  
	  /*------------------- DATA MEMBERS ---------------------*/
      Double_t *volt   ;			        //!
      Double_t *time  ;			            //!
	  double   Vbias   ;			        //!

	  /* Histograms */
	  TH1D *hvt    ;                        //! Histo of Voltage vs time
	  TH1D *hbl    ;		        //! Only bline
	  TH1D *hpbl   ;			//! Projected bline (vplot of bline)
	  TH1D *hvtbl  ;			//! Voltage vs time, baseline corrected
	  
	  double LPower ;                       //Laser power calculated as: Ipd [A]/0.77 [A/W], 
	                                        //where Ipd=Photocurrent is simply Vmax [V]/50 [Ohm]
	  double LNph ;                         //Number of photons injected in photodiode: Int(Vpd*dt)/(50*0.77*Ephoton)
	  
	  int      GetNbins()	 ;  //Ner of time points      //!
	  int      GetPolarity() ;			      //!
	  Double_t GetAbsVmax()  ;			      //!
	  double   GetVmax()	 ;			      //!
	  double   GetVmin()	 ;			      //!
	  double   GettVmax()	 ;			      //!
	  double   GettVmin()	 ;			      //!
	  double   BlineGetMean();			      //!
	  double   BlineGetRMS() ;			      //!
	  double   GetTleft( )   ;			      //!
	  double   GetTright( )  ;			      //!
	  double   GetTrms()	 ;			      //!
	  double   GetRiseTime( double fraction=0.9 ) ;	      //!
	  Double_t GetCharge( double tinf , double tsup ) ;   
	  Double_t RGetCharge( double tinf , double tsup ) ;   
          void     CalcRunningCharge( TMeas *em ) ;           //!
	  
	  

  private:						     
          int    Nbins     ;				    //!
	  int    Polarity  ;				    //!
	  
	  Double_t Vmax      ;
	  Double_t Vmin      ;
	  int      iVmax     ;
	  int      iVmin     ;
	  Double_t tVmax     ;
	  Double_t tVmin     ;
	   
	  double tleft     ;
	  double tright    ;
	  double trms      ;
	  int    itleft    ;
	  int    itright   ;
	  
	  double BlineMean ; 
	  double BlineRMS  ;
	  
	  Double_t Q50       ;   //From [tleft,tright+10 ns]
	  Double_t Qtot      ;   //From [time[0],time[Nt]]
	  
	  double RiseTime  ;
	  
	  
	  void   CalcVmaxmin()  ;						       //!
	  void   CalcPolarity() ;						       //!
	  void   CalcBline()  ; 						       //!
	  double CalcRiseTime( double fraction=RTIME )  ;			       //!
	  void   CreateHistos() ;						       //!		       //!
   
  /*protected:  //classes that inherit from TWaveform can access these methods*/
    
	  void   CalcSignalTimeL( double fraction , double &TimeL , int &iTimeL )  ;   //!
	  void   FitSignalTimeL( double &TimeL , int &iTimeL )  ;   //!
	  void   CalcSignalTimeR( double fraction , double &TimeR , int &iTimeR )  ;   //!
          void   CalcQTimeR(  double fraction , double &TimeR , int &iTimeR ) ;        //!
	  void   FitSignalTimeR( double &TimeR , int &iTimeR )  ;   //!
	  void   CalcSignalTimeLR( )  ;			                               //!
	  
  ClassDef(TWaveform,3) ; //ROOT RTTI
	
};

#endif
