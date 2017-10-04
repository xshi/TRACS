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

/* 
  Known problems:
  File: 20121108094858_Micron_FZ2852-23_PonN.ZVscan.root
  tleft evaluation using a 5%-25% fit of risetime does not work
    tree->Draw("volt:time","Vbias==10 && event==473","l")
  tleft ends up being:
     1) too big->BlineMean catches signal which does not go to Q50->Q50 drops
     2) or too small (maybe 0) and Q50 has no meaning 

*/
#include <stdlib.h>
#include <iostream>

#include "TWaveform.h"
#include <TMath.h>
#include <TVectorD.h>
#include <TGraph.h>
#include <TF1.h>

#define VMAX -999999.9
#define CHARSZ 256

Double_t T4BL=0.5; ;

ClassImp(TWaveform)

const int LOW  = -1;
const int HIGH = 15;

/* Default constructor, see ROOT manual page 270 */

TWaveform::TWaveform(  ) {
  Nbins = 0 ;
  time  = nullptr ;
  volt  = nullptr ;
  Vbias = 0.;
  itleft= Nbins ;
  hvt  = 0 ;
  hbl  = 0 ; 
  hpbl = 0 ; 
  hvtbl= 0 ;
  LPower=0.;
  LNph  =0.;
  tVmax = 0.;
  Qtot = 0.;
  Q50 = 0.;
  BlineRMS = 0.;
  tright = 0.;
  iVmin = 0.;
  iVmax = 0.;
  RiseTime = 0.;
  BlineMean = 0.;
  itright = 0.;
  trms = 0.;
  tleft = 0.;
  tVmin = 0.;
  Vmin = 0.;
  Vmax = 0.;
  Polarity = 0.;


}

TWaveform::TWaveform( TMeas *em ) {
   time = new Double_t [em->Nt] ;
   volt = new Double_t [em->Nt] ;

   for ( int i=0 ; i < em->Nt ; i++ ) {
     time[i]= em->time[i];
     volt[i]= em->volt[i];
   }

   Nbins=em->Nt ;

   Vbias=em->Vbias ;

   CalcVmaxmin( ) ;

   CalcPolarity () ;

   /* 
   Note that we need a rough estimation of itleft to calculate BlineMean. 
   After knowing BlineMean we can properly calculate the left and right extremes of the
   signal, and actually the BlineMean
   
   Taking  itleft=iVmax/2 does not work for very wide pulses (like in N bulk, low Vbias)
   Take 5 ns worth of data
   
   */
   
   T4BL = (em->Setup==5)? 0.5 : 10. ;
   if ( Polarity == 1 )   //NEW: Sept 2012
          itleft=( 0<iVmax && iVmax<=Nbins )? TMath::Nint( T4BL/(time[1]-time[0]) ) : Nbins ;
   else 
          itleft=( 0<iVmin && iVmin<=Nbins )? TMath::Nint( T4BL/(time[1]-time[0]) ) : Nbins ;
    
   CalcBline( ) ;

   CalcSignalTimeLR( ) ;

   CreateHistos() ;

   CalcRiseTime( RTIME ) ;
   
   Q50 = GetCharge( tleft   , tleft+25.  ) ;
   Qtot= GetCharge( time[0] , time[Nbins-1] ) ;

   CalcRunningCharge( em ) ;
   
}

TWaveform::TWaveform( int nel , double *tin , double *vin , double Bias  ) {

   time = new double [nel] ;
   volt = new double [nel] ;
   for ( int i=0 ; i < nel ; i++ ) {
     time[i]=tin[i];
     volt[i]=vin[i];
   }

   Nbins=nel ;

   Vbias=Bias ;

   CalcVmaxmin( ) ;

   CalcPolarity () ;

   /* Note that we need a rough estimation of itleft to calculate BlineMean. 
   After knowing BlineMean we can properly calculate the left and right extremes of the
   signal, and actually the BlineMean 
   */

   //itleft=( 0<iVmax && iVmax<=Nbins )? (int) iVmax/2 : Nbins ;
   if ( Polarity == 1 )   //NEW: Sept 2012
          itleft=( 0<iVmax && iVmax<=Nbins )? TMath::Nint(  T4BL/(time[1]-time[0]) ) : Nbins ;
   else 
          itleft=( 0<iVmin && iVmin<=Nbins )? TMath::Nint(  T4BL/(time[1]-time[0]) ) : Nbins ;

   CalcBline( ) ;

   CalcSignalTimeLR( ) ;

   CreateHistos() ;

   CalcRiseTime( RTIME ) ;

   Q50 = GetCharge( tleft   , tleft+25.  ) ;
   Qtot= GetCharge( time[0] , time[Nbins-1] ) ;
      
}

//Destructor
TWaveform::~TWaveform(){
  //if (volt != (Double_t *) 0) delete [] volt ; //19Sept2012
  //if (time != (Double_t *) 0) delete [] time ;
  delete [] volt ;
  delete [] time ;
  delete hvt   ;  
  delete hbl   ; 
  delete hpbl  ; 
  delete hvtbl ;
}

/*------------------- PUBLIC METHODS -----------------------------*/ 
/*------------------- GETTERS         -----------------------------*/ 
int    TWaveform::GetNbins(){
  return Nbins ;
}

Double_t TWaveform::GetVmax(){
  return Vmax ;
}

Double_t TWaveform::GetVmin(){
  return Vmin ;
}

Double_t TWaveform::GettVmax(){
  return tVmax ;
}

Double_t TWaveform::GettVmin(){
  return tVmin ;
}

Double_t   TWaveform::GetAbsVmax(){
  Double_t val = (fabs(Vmax)>fabs(Vmin)) ? Vmax : Vmin  ;
  return val ;
}

int    TWaveform::GetPolarity(){
  return Polarity ;
}

double TWaveform::GetTleft() { 
  return tleft ; 
}

double TWaveform::GetTright(){ 
  return tright; 
}

double TWaveform::GetTrms()  { 
  return trms  ; 
}

double TWaveform::BlineGetMean() {
  return BlineMean ;
}

double TWaveform::BlineGetRMS() {
  return BlineRMS ;
}

double TWaveform::GetRiseTime( double fraction )  { 
  double rt=CalcRiseTime( fraction );
  return rt  ; 
}

Double_t TWaveform::GetCharge(double tinf , double tsup ) {
  if (hvtbl==0) { 
    cout << "No histogram"<<endl;
    exit(0) ;
  }
  int imin=hvtbl->FindBin(tinf);
  int imax=hvtbl->FindBin(tsup);
  //cout << "imin="<<imin<<" imax" << imax << endl ;
  return hvtbl->Integral(imin,imax,"width") ;
}

Double_t TWaveform::RGetCharge(double tinf , double tsup ) {
//Estudiar como acceder a los datos llamando a este metodo desde root:
//    tree->Draw("RGetCharge(100.,200.)")

  cout << "Nt="<< Nbins <<endl;
//   CreateHistos();
//   cout << "Histos Created" << endl ;
//   int imin=hvtbl->FindBin(tinf);
//   int imax=hvtbl->FindBin(tsup);
//   cout << "imin="<<imin<<" imax" << imax << endl ;
//   return hvtbl->Integral(imin,imax,"width") ;
  return 0.;
}


/*------------------- PRIVATE METHODS -----------------------------*/

void   TWaveform::CalcVmaxmin() {
  iVmax=TMath::LocMax(Nbins,volt) ;
  iVmin=TMath::LocMin(Nbins,volt) ;
  Vmax=volt[iVmax];
  Vmin=volt[iVmin];
  tVmax=time[iVmax];
  tVmin=time[iVmin];
}

#ifdef OLD
void   TWaveform::CalcPolarity(){

  //Gives +1 if signal is positive, -1 otherwise
  if ( Vmax==VMAX ) CalcVmaxmin() ;
  double avx=fabs(Vmax) ;
  double avn=fabs(Vmin) ;
  Polarity = ( avx > avn ) ? 1 : -1 ;
}
#endif

#ifdef OLD1 //Failed for 2014-05-20_09-00-10_Liv-2935-7-1-15-L
void   TWaveform::CalcPolarity(){

  /*Gives +1 if signal is positive, -1 otherwise */
  
  //Estimate baseline
  Int_t nbl=TMath::Nint(  T4BL/(time[1]-time[0]) );
  TVectorD myblv(nbl,volt) ;
  Double_t mybl = myblv.Sum()/nbl ;
  
  TVectorD voltbl(Nbins,volt) ;
  voltbl.AddSomeConstant(-1.0*mybl,voltbl) ;
  
  double sum = voltbl.Sum() ;
  Polarity = ( sum > 0. ) ? 1 : -1 ;
  //cout  << "Polarity" << Polarity << endl;

}
#endif

void   TWaveform::CalcPolarity(){

  /*Gives +1 if signal is positive, -1 otherwise */
    
  //Estimate baseline
  Int_t nbl=TMath::Nint(  T4BL/(time[1]-time[0]) ) ;
  TVectorD myblv(nbl,volt) ;
  Double_t mybl = myblv.Sum()/nbl ;
  
  TVectorD voltbl(Nbins,volt) ;
  voltbl.AddSomeConstant(-1.0*mybl,voltbl) ;
  
  double avx = voltbl.Max( ) ;
  double avn = voltbl.Min( ) ;
  Polarity = ( fabs(avx) > fabs(avn) ) ? 1 : -1 ;
  //cout  << "Polarity" << Polarity << endl;

}


void TWaveform::CalcSignalTimeL( double fraction , double &TimeL , int &iTimeL ) {

  /* Identify entries that are above a fraction of the maximum. Discriminate
     those voltages (bline corrected) above a threshold given as fraction*Vmax.
     If abovethold=0, we have the left extreme of the figure
  */

  Double_t thold , abovethold ;
  thold = (Polarity==1)? BlineMean+fraction*(Vmax-BlineMean) : BlineMean+fraction*(Vmin-BlineMean) ;
  int iPos = (Polarity==1)? iVmax : iVmin ;
  
  TimeL=0 ; iTimeL=0 ; //For the case where the loop is not visited  
  for (int i=iPos-1 ; i>0 ; i--) {     //Very few iterations needed
    //Double_t diff=Polarity*(volt[i]-BlineMean-thold) ;
    Double_t diff=Polarity*(volt[i]-thold) ;
    abovethold=( diff>0.)? 1:0 ;
    if ( abovethold==0 ) { 
      TimeL = ( i+1 < Nbins && (volt[i+1]-volt[i])!=0. )?
              time[i] + (time[i+1]-time[i])/(volt[i+1]-volt[i])*(thold-volt[i]):
              time[i] ;
      iTimeL=i;
      break ;
    }
  } 
    
}

void TWaveform::FitSignalTimeL( double &TimeL , int &iTimeL ) {

  /* Identify entries that are above a fraction of the maximum. Discriminate
     those voltages (bline corrected) above a threshold given as fraction*Vmax.
     If abovethold=0, we have the left extreme of the figure
  */

  Double_t thold , abovethold ;
  int iPos = (Polarity==1)? iVmax : iVmin ;
  
  //Look for the moment we cross the Baseline+25% and Baseline +5%
  //thold = (Polarity==1)? BlineMean+0.85*(Vmax-BlineMean) : BlineMean+0.85*(Vmin-BlineMean) ;
  thold = (Polarity==1)? BlineMean+0.25*(Vmax-BlineMean) : BlineMean+0.25*(Vmin-BlineMean) ;
  for (int i=iPos-1 ; i>0 ; i--) {     //Very few iterations needed
    //Double_t diff=Polarity*(volt[i]-BlineMean-thold) ;
    Double_t diff=Polarity*(volt[i]-thold) ;  //NEW: 9 Nov 2012
    abovethold=( diff>0.)? 1:0 ;
    if ( abovethold==0 ) { 
      TimeL=time[i] ;
      iPos=i;
      break ;
    }
  }
  
  thold = (Polarity==1)? BlineMean+0.05*(Vmax-BlineMean) : BlineMean+0.05*(Vmin-BlineMean) ;

  //Protection against cases where thold falls below the BlineRMS
  if (Polarity*thold<BlineRMS) thold = (Polarity==1)? BlineMean+0.15*(Vmax-BlineMean) : BlineMean+0.15*(Vmin-BlineMean);
  
  for (int i=iPos-1 ; i>0 ; i--) {     //Very few iterations needed
    //Double_t diff=Polarity*(volt[i]-BlineMean-thold) ;
    Double_t diff=Polarity*(volt[i]-thold) ;  //NEW: 9 Nov 2012
    abovethold=( diff>0.)? 1:0 ;
    if ( abovethold==0 ) { 
      TimeL=time[i] ;
      iTimeL=i;
      break ;
    }
  }
  
  //These should be "no signal" cases. I leave the "naive" estimation as iPos/2 as the BlineMean
  if ( (iTimeL>iPos) || ( iPos-iTimeL+1<3 ) ) {
    //TimeL  = 0.0  ; iTimeL = 0.0  ;
    return ;
  }
  
  //Fit the slope of volt:time from the Baseline to the middle of the Vmax
  TF1 *fpol1 = new TF1("fpol1","pol1",time[iTimeL],time[iPos]);
  TGraph *gr = new TGraph( iPos-iTimeL+1, &time[iTimeL] , &volt[iTimeL] );
  #define LPLOTFIT 0
  #if LPLOTFIT==1
    gr->Fit("fpol1","RQ");
  #else
    gr->Fit("fpol1","RQ0");
  #endif
  Double_t offset=fpol1->GetParameter(0);
  Double_t  slope=fpol1->GetParameter(1);
 
  TimeL  = (slope!=0.)? (BlineMean-offset)/slope : 0.0 ;
  iTimeL = TMath::Nint( (TimeL-time[0])/(time[1]-time[0]) ) ;
  
  #if LPLOTFIT==1
    TGraph *grt = new TGraph( Nbins, &time[0] , &volt[0] );
    grt->Draw("awlp") ; fpol1->Draw("lsame") ;
    grt->Write(); fpol1->Write();
    delete grt;
  #endif

  delete fpol1;
  delete gr;
  
  //cout << "tleft="<<TimeL<<endl;
    
}

void TWaveform::FitSignalTimeR( double &TimeR , int &iTimeR ) {

  /* Identify entries that are above a fraction of the maximum. Discriminate
     those voltages (bline corrected) above a threshold given as fraction*Vmax.
     If abovethold=0, we have the left extreme of the figure
  */

  Double_t thold , abovethold ;
  int iPos = (Polarity==1)? iVmax : iVmin ;
  
  //Look for the moment we cross the Baseline+25% and Baseline +5%
  //thold = (Polarity==1)? BlineMean+0.45*(Vmax-BlineMean) : BlineMean+0.45*(Vmin-BlineMean) ;
  thold = (Polarity==1)? BlineMean+0.5*(Vmax-BlineMean) : BlineMean+0.5*(Vmin-BlineMean) ;
  for (int i=iPos+1 ; i<Nbins ; i++) {     //Very few iterations needed
    //Double_t diff=Polarity*(volt[i]-BlineMean-thold) ;
    Double_t diff=Polarity*(volt[i]-thold) ;  //NEW: 9 Nov 2012
    abovethold=( diff>0.)? 1:0 ;
    if ( abovethold==0 ) { 
      TimeR = ( i-1 >= 0 && (volt[i-1]-volt[i])!=0. )?
              time[i] + (time[i-1]-time[i])/(volt[i-1]-volt[i])*(thold-volt[i]):
              time[i] ;
      iPos=i;
      break ;
    }
  }
  thold = (Polarity==1)? BlineMean+0.05*(Vmax-BlineMean) : BlineMean+0.05*(Vmin-BlineMean) ;

  //Protection against cases where thold falls below the BlineRMS
  if (Polarity*thold<BlineRMS) thold = (Polarity==1)? BlineMean+0.15*(Vmax-BlineMean) : BlineMean+0.15*(Vmin-BlineMean);

  for (int i=iPos-1 ; i<Nbins ; i++) {     //Very few iterations needed
    //Double_t diff=Polarity*(volt[i]-BlineMean-thold) ;
    Double_t diff=Polarity*(volt[i]-thold) ;  //NEW: 9 Nov 2012
    abovethold=( diff>0.)? 1:0 ;
    if ( abovethold==0 ) { 
      TimeR = ( i-1 >= 0 && (volt[i-1]-volt[i])!=0. )?
              time[i] + (time[i-1]-time[i])/(volt[i-1]-volt[i])*(thold-volt[i]):
              time[i] ;
      iTimeR=i;
      break ;
    }
  }
  
  //These should be "no signal" cases. I leave the "naive" estimation as iPos/2 as the BlineMean
  if ( (iTimeR<iPos) || ( iTimeR-iPos+1<3 ) || (iTimeR>Nbins-1) ) {
    return ;
  }
  
  //Fit the slope of volt:time from the Baseline to the middle of the Vmax
  TVectorD voltbl(Nbins,volt) ;
  TVectorD timed(Nbins,time) ;
  voltbl.AddSomeConstant(-1.0*BlineMean,voltbl) ;
  TGraph *gr = new TGraph( timed , voltbl );
  TF1 *fpol1 = new TF1("fpol1","pol1",time[iPos],time[iTimeR]);
  
  #define RPLOTFIT 0
  #if RPLOTFIT==1
    gr->Fit("fpol1","RQ");
  #else
    gr->Fit("fpol1","RQ0");
  #endif
  Double_t offset=fpol1->GetParameter(0);
  Double_t  slope=fpol1->GetParameter(1);
 
  TimeR  = (slope!=0.)? (BlineMean-offset)/slope : 0.0 ;
  iTimeR = TMath::Nint( (TimeR-time[0])/(time[1]-time[0]) ) ;
  
  //cout << "tleft="<<TimeL<<endl;

  #if RPLOTFIT==1
    gr->Draw("awlp") ; fpol1->Draw("lsame") ;
    gr->Write(); fpol1->Write();
  #endif
  
  delete fpol1;
  delete gr;
  timed.Clear();
  voltbl.Clear();
    
}

void TWaveform::CalcSignalTimeR( double fraction , double &TimeR , int &iTimeR) {

  /* Identify entries that are above a fraction of the maximum. Discriminate
     those voltages (bline corrected) above a threshold given as fraction*Vmax.
     From iVmax to the right, once abovethold=0, we have the right extreme of the figure
  */

  Double_t thold , abovethold ;
  thold = (Polarity==1)? BlineMean+fraction*(Vmax-BlineMean) : BlineMean+fraction*(Vmin-BlineMean) ;
  int iPos = (Polarity==1)? iVmax : iVmin ;
    
  for (int i=iPos+1 ; i<Nbins ; i++) { //Very few iterations here
    //Double_t diff=Polarity*(volt[i]-BlineMean-thold) ;
    Double_t diff=Polarity*(volt[i]-thold) ;
    abovethold=( diff>0.)? 1:0 ;
    if ( abovethold==0 ) { 
      TimeR = ( i-1 >= 0 && (volt[i-1]-volt[i])!=0. )?
              time[i] + (time[i-1]-time[i])/(volt[i-1]-volt[i])*(thold-volt[i]):
              time[i] ;
      iTimeR=i;
      break ;
    }
    
  }
       
}

void TWaveform::CalcQTimeR(  double fraction , double &TimeR , int &iTimeR ) {

  /* Calculate the total charge of the pulse, then look for the time at which 
     Q=fraction*Qtot. Example: fraction is 98%
     NEW: sometime [for instance simulated cases], the pulses are longer than 25 ns
     Therefore use an estimation of tright of the signal, coming from CalcSignalTimeR, 
     for instance.
  */

 
  //Voltage baseline corrected
  TVectorD voltbl(Nbins,volt) ;
  voltbl.AddSomeConstant(-1*BlineMean,voltbl) ;
    
  // We calculate Qtot=int(tleft,tleft+tright), where tright is pre-estimated
  //Int_t Nb25   = TMath::Nint( 20. / (time[1] - time[0]) ) ;  
  Int_t Nb25   = (iTimeR-itleft+10>0) ? iTimeR-itleft+10 : 2 ;  
  Int_t iqlast = ( itleft + Nb25 > Nbins-1 ) ? Nbins-1 : itleft +  Nb25 ;
  Nb25 = iqlast - itleft ;
  TVectorD voltblr( Nb25 ) ; 
  voltblr = voltbl.GetSub( itleft ,  iqlast-1 ) ;

  Double_t Qt = voltblr.Sum() ;
  Double_t Qfrac = fraction * Qt ;  
    
  Double_t Qsum = Qt ;
  Int_t i = Nb25 ;  
  //cout << Nb25 << " " << iqlast << " " << itleft << endl ;
  while ( TMath::Abs(Qsum) > TMath::Abs(Qfrac) ) {
     
    i-- ;
    if ( (Polarity==1 && voltblr[i]>0.) || (Polarity==-1 && voltblr[i]<0.) ) Qsum-=voltblr[i] ;
    //Qsum-=voltblr[i] ;
    if (i==0) break ;
  }
  TimeR=time[itleft + i] ;
  iTimeR=itleft + i;
  
  voltblr.Clear() ; voltbl.Clear() ;
       
}

void TWaveform::CalcSignalTimeLR( ) {

  /* Identify entries that are above a fraction (normally half) of the maximum. 
     Coming from small times, the 1st time we reach this value, we will be
     in the rising of the pulse. The second time will be the falling edge. In PAW
     I also fitted a P0+G and then extracted the RMS.
  */
  
    //For the RMS of the signal, it is better to study half of the width of the signal
    CalcSignalTimeL(0.5 , tleft  , itleft  ) ;
    CalcSignalTimeR(0.5 , tright , itright ) ;
    trms   = (tright-tleft)/2.35 ;
    //tleft  = tleft - 2.0*trms ;
    //tright = tright+ 6.0*trms ;
    
    //Now, for the extremes of the signal, we look for the place where it crosses the 
    //baseline
    FitSignalTimeL( tleft  , itleft  ) ;
    Double_t tstep = time[1]-time[0];
    itleft  = TMath::Nint(  tleft/tstep ) ;
    
    //Protection against no signal cases
    if (itleft<0 || itleft>Nbins) {
      itleft=0 ; tleft=-1. ;
      //itleft=Nbins/2 ; tleft=time[Nbins/2] ;
    }
    
    //Right extreme
    if (tleft>=0) {
    
      CalcSignalTimeR(0.1 , tright , itright ) ;
      CalcQTimeR(0.98 , tright , itright ) ;
      //FitSignalTimeR( tright  , itright  ) ;
      itright = TMath::Nint( tright/tstep ) ;

      if (itright<0 || itright>Nbins) {
	itright=0 ; tright=-2.0 ;
      }
      
    } else { 
      itright=0 ; tright=-2.0 ;
    }    

    //Now that we know where the signal is, calculate the Bline properly
    CalcBline() ;
    
}

void TWaveform::CalcBline() {

  int N5ns = TMath::Nint( T4BL/(time[1]-time[0]) );
  if ( itleft>0 ) {
    //If the pulses are not monotonous (increasing or decreasing), itleft can be wrong. 
    //If itleft falls in the signal, then BlineMean is wrong. We use BlineMean to correct 
    //many magnitudes.
    BlineMean = ( N5ns < itleft ) ? TMath::Mean( N5ns , volt ) : TMath::Mean( itleft , volt );
    BlineRMS  = ( N5ns < itleft ) ? TMath::RMS( N5ns , volt )  : TMath::RMS(  itleft , volt );
  } else {
    BlineMean = TMath::Mean( (int) Nbins , volt );
    BlineRMS  = TMath::RMS(  (int) Nbins , volt );    
  }
  //cout << "BlineMean["<<itleft<<"]=" << BlineMean ;
  //cout << " BlineRMS[" <<itleft<<"]=" << BlineRMS << endl ;
}

double TWaveform::CalcRiseTime( double fraction ) {
  double rtl  , rtr  ;
  int    irtl , irtr ;
  CalcSignalTimeL( 1.0-fraction  , rtl , irtl );
  CalcSignalTimeL( fraction      , rtr , irtr );
  RiseTime = rtr - rtl ;
  return RiseTime ;
}

void TWaveform::CreateHistos( ) {

  double tmin = TMath::MinElement( Nbins , time ) ;
  double tmax = TMath::MaxElement( Nbins , time ) ;
  double At   = (tmax-tmin)/(Nbins-1.0) ;
  tmin = -0.5*At ;
  tmax = tmin + Nbins*At ;
  
  //double xpos  = xyz.X() , ypos=xyz.Y(), zpos=xyz.Z() ;
  double xpos  = 0 , ypos=0, zpos=0 ;
  
  char hid[CHARSZ] , htxt[CHARSZ] ;
  
  //Voltage vs time, as seen in the scope
  sprintf( hid , "(%f,%f,%f),Vb=%f" ,xpos,ypos,zpos,Vbias ) ;
  sprintf( htxt , "V(t),(%f,%f,%f),Vb=%f" ,xpos,ypos,zpos,Vbias ) ;
  hvt=new TH1D( hid , htxt , Nbins , tmin , tmax ) ;
  hvt->FillN(Nbins, time , volt ) ;

  //Only baseline
  sprintf( hid , "BL (%f,%f,%f),Vb=%f" ,xpos,ypos,zpos,Vbias ) ;
  sprintf( htxt , "Baseline,(%f,%f,%f),Vb=%f" ,xpos,ypos,zpos,Vbias ) ;
  hbl=new TH1D( hid , htxt , itleft , tmin , tmin + itleft*At ) ; //As long as it starts in bin 1
  hbl->FillN (itleft, time, volt ) ;
  
  //Projected baseline (v/hfill equivalent)
  //http://www.yolinux.com/TUTORIALS/LinuxTutorialC++STL.html
  vector<double> weight(itleft);
  weight.assign(itleft,1.0);
  double bmin=TMath::MinElement( itleft , volt ) ;
  double bmax=TMath::MaxElement( itleft , volt ) ;
  double Ab=0.1*(bmax-bmin) ;
  sprintf( hid , "PBL (%f,%f,%f),Vb=%f" ,xpos,ypos,zpos,Vbias ) ;
  sprintf( htxt , "1D Baseline,(%f,%f,%f),Vb=%f" ,xpos,ypos,zpos,Vbias ) ;
  hpbl=new TH1D(hid , htxt , 100, bmin-Ab , bmin+Ab );
  hpbl->FillN (itleft, volt , &weight[0] ) ;
  weight.clear() ;

  //Voltage vs time, baseline corrected
  TVectorD voltbl(Nbins,volt) ;
  voltbl.AddSomeConstant(-1*BlineMean,voltbl) ;
  sprintf( hid , "V(t)-bl (%f,%f,%f),Vb=%f" ,xpos,ypos,zpos,Vbias ) ;
  sprintf( htxt , "BaselineCorr,(%f,%f,%f),Vb=%f" ,xpos,ypos,zpos,Vbias ) ;
  hvtbl=new TH1D( hid , htxt , Nbins , tmin , tmax ) ;
  //hvtbl=(TH1D*)hbl->Clone() ;
  //hvtbl->Reset();
  hvtbl->FillN(Nbins, time , &voltbl[0] );
  
    
  //hvt->Write() ; hbl->Write() ; hpbl->Write() ; hvtbl->Write() ;
  
}

void TWaveform::CalcRunningCharge( TMeas *em ) {
  
  Int_t ir = 1 ;
  em->Qt[0] = volt[0]- BlineMean ;
  while ( ir < Nbins ) { 
    em->Qt[ir]= em->Qt[ir-1] + volt[ir] - BlineMean ;
    ir++;
  }
  
    
}
