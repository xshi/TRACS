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

   Function TTreeChi2 assumes that both, measurement and simulation share
   the same (X,Y,Z) range. The time of the waveforms does not need to be the same

*/


#ifndef TRACSFIT_H
#define TRACSFIT_H

#include <stdlib.h>
#include <iostream>
#include <thread>
/* Note: With the simulation we should be able to estimate the Neh */
#include <Minuit2/FCNBase.h>
#include <Minuit2/MnConfig.h>
#include <Minuit2/GenericFunction.h>
#include <Fit/FitResult.h>
#include <TVirtualFitter.h>
#include <TMath.h>
#include <TTreeFormula.h>
#include <TEventList.h>
#include <TEntryList.h>
#include <TEntryListArray.h>
#include <Math/Interpolator.h>

#include <TMeas.h>
#include <TWaveform.h>
#include <Threading.h>
#include <TRACSInterface.h>


//using namespace std;       //Avoids using std::cin;

//template <class Function>; //Chapter 5.7.5.2 of ROOT User's Guide

class TRACSFit : public ROOT::Minuit2::FCNBase {

  public:

     TRACSFit();

     TRACSFit ( TString FileMeas , TString FileConf , TString howstr ) ;

     ~TRACSFit(  ) ;

     Double_t MyTRACSFit( TH1D *hsimc , TH1D *hmeasc );

     Double_t LeastSquares( ) ;

     virtual  Double_t operator( )(const std::vector<Double_t>& par  ) const ;

     virtual  Double_t Up() const {return 1.;}

     void MinimumCommonHistogram( TH1D *hsim, TH1D *hmeas ) ;

#ifdef FUTURE

     void SetParameter(  Int_t ipar, Double_t val ) ;

     void SetParameters( vector<Double_t> vpar ) ;

     void FixParValue() ;

     void SetGlobalCondition( TString str_gcond ) { how = str_gcond } ;

     //Limit Fit to this range of coordinate(s)
     void SetCoordRange( Int_t icoord  , Double_t cmin , Double_t cmax ) ;
     void SetCoordRange( TString coord , Double_t cmin , Double_t cmax ) ;

     //Limit Fit to event range
     void SetEventRange( Int_t imin , imax ) ;

     //Ban a region in space for Fit (maybe too close to boundaries?)
     void BanCoordRange( TString coord , Double_t cmin , Double_t cmax );


     //Minimization inside ROOT
     Int_t FitTCTSignal( TH1D *h1  ) ;

     Int_t FitTCTSignal( TTree *tree , TString how ) ;

     Int_t SetScaleToken( int_t val ) { SameScale = val ; } ;
     Int_t SetInputToken( int_t val ) { H1orTree = val ; } ;
     Int_t GetScaleToken() { return SameScale ; } ;
     Int_t GetInputToken() { return H1orTree ; } ;

#endif


  private:

     TString how ;          //Something like "Vbias==-80 && (z>0.05 && z<0.25)"

     //std::string file_conf = "Config.TRACS";
     //std::vector<TRACSInterface> * sim = nullptr ;  //Object holding the simulation of TRACS

     /* Measurement   and  Simulation TREEs info */
     TTree * tmeas, * tsim;
     TMeas * em, * ems;
     TMeasHeader * emh, * emhs;
     TWaveform * wv, * wvs;



     TEntryListArray * listm, *lists; //Sets of events fulfilling condition "how"
     Int_t Nevm, Nevs; //Number of events fulfilling "how" condition
     Int_t ntm, nts   ; //Dimension of tmeas and tsim
     Int_t imins, imaxs, iminm, imaxm ; //Common maximum and minimum indexes

     Double_t theChi2, fitNorm;

     vector<Double_t> neffArray;

     Int_t  sameScale;       //Token that informs if input histograms are or not
     			      //in the same X axis range

     Int_t  h1orTree;        //Token that informs if input is histo or tree

     TH1D * hsim; //TRACS simulated induced current
     TH1D * hmeas; //Measurement
     TH1D * hsimc; //hsim in a time range compatible with the measurement
     TH1D * hmeasc; //Measurement in a time range compatible with the simulation

};
#endif
