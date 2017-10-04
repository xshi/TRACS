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

#include <TRACSFit.h>

#include <cmath>
#include <linear.h>

//ClassImp(TRACSFit)
ClassImp(TMeas)
ClassImp(TMeasHeader)
ClassImp(TWaveform)

extern std::vector<TRACSInterface*> TRACSsim;
extern std::vector<std::thread> t;
/* Default constructor, see ROOT manual page 270 */

/**
 *
 */
TRACSFit::TRACSFit() {
	imaxs = 0;
	imaxm = 0;
	iminm = 0;
	imins = 0;
	sameScale = 0 ;
	how = "";
	tsim = nullptr;
	hmeas = nullptr;
	hsimc = nullptr;
	hsim = nullptr;
	hmeasc = nullptr;
	h1orTree = 0;
	theChi2 = 0;
	Nevs = 0;
	Nevm = 0;
	listm = nullptr;
	ems = 0;
	em = 0;
	wvs = nullptr;
	wv = nullptr;
	lists = nullptr;
	tmeas = nullptr;
	nts = 0;
	ntm = 0;
	emh = 0;
	emhs = 0;
	fitNorm = 0;
}
/**
 *
 * @param FileMeas
 * @param FileConf
 * @param howstr
 */
TRACSFit::TRACSFit( TString FileMeas , TString FileConf ,  TString howstr ) {

	//Create simulation tree
	ntm = 0;
	nts = 0;
	how = howstr ;
	tsim = nullptr;
	hmeas = nullptr;
	hsimc = nullptr;
	hsim = nullptr;
	hmeasc = nullptr;
	sameScale = 0;
	h1orTree = 0;
	theChi2 = 0;
	fitNorm = 0;
	//em = 0;
	//ems = 0;
	/*-------------  M E A S U R E M E N T   TREE --------------------------*/
	//Get MEASUREMENT tree
	//Good except tmeas->GetEntry( iev );
	//Normalization constant
	//fitNorm = TRACSsim[0]->get_fitNorm();

	TFile *fmeas = new TFile( FileMeas.Data() );
	tmeas = (TTree*) fmeas->Get("edge");

    em   = 0 ;
    TBranch *raw  = tmeas->GetBranch("raw") ;
    raw->SetAddress(&em) ;

    raw->GetEntry(0) ;

    wv = new TWaveform( em ) ;
    TBranch *proc = tmeas->GetBranch("proc") ;
    proc->SetAddress(&wv) ;
    emh = (TMeasHeader *) tmeas->GetUserInfo()->At(0) ;
    tmeas->GetEntry(0) ;  //Commented 19Sept2012


	//MEASUREMENT: Subset of entries fulfilling "how" condition
	tmeas->Draw( ">>myListMeas" , (char*) how.Data() , "entrylistarray" ) ;
	listm=(TEntryListArray*)gDirectory->Get("myListMeas") ;
	tmeas->SetEntryList( listm ) ; //Use tree->SetEventList(0) to switch off

	//MEASUREMENT: Time vector
	Nevm = listm->GetN() ;
    if ( Nevm == 0 ) {
       std::cout<<"Error: No events in measurement file satisfy condition " << how.Data() << std::endl ;
       std::cout<<"Exiting!"<< std::endl ;
       exit(-1);
    }
    Int_t iev = listm->GetEntry(0) ;
	tmeas->GetEntry( iev );
	ntm = wv->GetNbins() ;
	std::vector<Double_t> timem(ntm);

	//Dump time information into timem
    for ( Int_t ii = 0 ; ii< ntm ; ii++ ) timem[ii] = em->time[ii] ;

	//delete em;



	/*-------------  S I M U L A T I O N    TREE --------------------------*/
	//Create objects TWaveform and TMeas for SIMULATION
	//Convert simulation data into a TTree
	tsim = new TTree("edge","TRACS simulation");
	TRACSsim[0]->GetTree( tsim );

    ems = 0 ;
	TBranch *raws  = tsim->GetBranch("raw") ;
	raws->SetAddress(&ems) ;
	raws->GetEntry(0) ;
	wvs = new TWaveform( ems ) ;
	TBranch *procs  = tsim->GetBranch("proc") ;
	procs->SetAddress(&wvs) ;
	emhs=0;

	//SIMULATION: Subset of entries fulfilling "how" condition
	//tsim->Draw("volt-BlineMean:time","event==0","l"); gPad->Modified();gPad->Update();
	tsim->Draw( ">>myListSim" ,  (char*) how.Data() , "entrylistarray" ) ;
	//tsim->Draw( ">>myListSim" ,  (char*) how.Data()  ) ;
	lists=(TEntryListArray*) gDirectory->Get("myListSim") ;
	tsim->SetEntryList( lists ) ; //Use tree->SetEventList(0) to switch off

	//SIMULATION: Get time vector
	iev = lists->GetEntry(0) ;
	tsim->GetEntry( iev );
	nts = wvs->GetNbins() ;
	vector<Double_t> tims(nts);
	Nevs  = lists->GetN() ;

	//Dump time information into tims
    for ( Int_t ii = 0 ; ii< nts ; ii++ ) tims[ii] = ems->time[ii] ;

	/*-------------  F I N D   C O M M O N   T I M E   R A N G E  --------------------------*/

	if (Nevs != Nevm ) {
		std::cout << "We need both the simulation and the measurement to be in the same (X,Y,Z;V) ranges" << std::endl ;
	}

	//... from here on, we assume Nev = Nevm
	//Double_t *volts;
	//Double_t *voltm;


	/*
  Find minimum common time range between sim and meas. Simulation Isim(t) does
  not have to be evaluated at the same time "t" as the measurement Imeas(t). We
  can interpolate the simulation but we need simulation and measurement defined
  in the same time range [tmin,tmax] (to avoid extrapolation).
	 */


	Double_t tmin = (tims[0]<=timem[0])?         timem[0] : tims[0]   ;
	Double_t tmax = (tims[nts-1]<=timem[ntm-1])? tims[nts-1] : timem[ntm-1] ;
	tsim->GetEntry( lists->GetEntry(0) );
	tmeas->GetEntry( listm->GetEntry(0) );
	Double_t Ats=wvs->time[2]-wvs->time[1] , Atm=wv->time[2]-wv->time[1];

	//Common maximum and minimum indexes
	imins = TMath::Nint( (tmin-tims[0])/Ats ) , imaxs =TMath::Nint( (tmax-tims[0])/Ats );
	iminm = TMath::Nint( (tmin-timem[0])/Atm ) , imaxm =TMath::Nint( (tmax-timem[0])/Atm );


}
//---------------------------------------------------------------------------

//Destructor
TRACSFit::~TRACSFit(){

	//delete sim   ;
	delete wv ;
	delete wvs ;
	delete tsim ;
	//delete [] times ; delete [] timem ;

}


//---------------------------------------------------------------------------

Double_t TRACSFit::LeastSquares( ) {

	/*
     Chi2 between events of the simulation and events of the measurement
     Both trees need to have common (X,Y,Z,V) coordinates.
     They can differ in the time resolution of the waveforms, that is, Nt
     (number of points in the waveform).
	 */

	tsim = new TTree("edge","TRACS simulation");
	TRACSsim[0]->GetTree( tsim );
    ems = 0 ;
	TBranch *raws  = tsim->GetBranch("raw") ;
	raws->SetAddress(&ems) ;
	raws->GetEntry(0) ;
	TWaveform *wvit = new TWaveform( ems ) ;
	TBranch *procs  = tsim->GetBranch("proc") ;
	procs->SetAddress(&wvit) ;
	emhs=0;
	tsim->SetEntryList(  lists ) ;
	tmeas->SetEntryList( listm ) ;


	vector<Double_t>  timem, voltm , tims, volts ;
    tims.resize(nts);
    volts.resize(nts);
    timem.resize(ntm);
    voltm.resize(ntm);
	Double_t chi2 = 0.;
	//Int_t Nt_iterm = imaxm - iminm + 1 ;
	for ( Int_t ii=0 ; ii < Nevm ; ii++ ) {

		//Retrieve simulated waveform fulfilling condition "how"
		Int_t iev = lists->GetEntry(ii) ;
		tsim->GetEntry( iev );

		for ( Int_t iv = 0 ; iv< nts ; iv++ ) {
			volts[iv] = ems->volt[iv] ;
			tims[iv]  = ems->time[iv] ;
		}

		//Prepare interpolator for simulation data in case it has not been simulated in the right time grid
		ROOT::Math::Interpolator itp1( tims , volts , ROOT::Math::Interpolation::kLINEAR);
		itp1.SetData(tims , volts);

		//Retrieve measured waveform fulfilling condition "how"
		iev = listm->GetEntry(ii) ;
		tmeas->GetEntry( iev );
		fitNorm = TRACSsim[0]->get_fitNorm();
		for ( Int_t iv = 0 ; iv< ntm ; iv++ ) {
			voltm[iv] = em->volt[iv] - wv->BlineGetMean();
			voltm[iv] = -1 *voltm[iv]; //Change sign of Meas *******
 			timem[iv] = em->time[iv] ;
		}
		double simulation; //, norm=TRACSsim[0]->get_fitNorm();
		for ( Int_t iv = iminm ; iv< imaxm ; iv++ ) {
			simulation = itp1.Eval(timem[iv]);

		if (simulation != 0) //For not to fit the 0's part!.********
			//std::cout << "Scan="<<ii<<" Iteration= "<<iv<<" m=" << voltm[iv] << " s=" << simulation << " X2=" << chi2 << std::endl;
			chi2+=( voltm[iv]-fitNorm*simulation )*(voltm[iv]-fitNorm*simulation) ;
		}


	}

	chi2 = chi2/(TRACSsim[0]->GetchiFinal()*TRACSsim[0]->GetchiFinal());

	//thechi2??
	delete wvit ;
	return chi2 ;

}

//---------------------------------------------------------------------------
/**
 *
 * @param hsim
 * @param hmeas
 */
void TRACSFit::MinimumCommonHistogram( TH1D *hsim, TH1D *hmeas ) {

	/*
      Moves hsim and hmeas to a common X axis.
      It has to be the minimum range in common
	 */
	Double_t Xsm = hsim->GetXaxis()->GetXmin()  , Xsx = hsim->GetXaxis()->GetXmax() ;
	Double_t Xmm = hmeas->GetXaxis()->GetXmin() , Xmx = hmeas->GetXaxis()->GetXmax();
	Double_t xmin  = (Xsm>=Xmm)? Xsm : Xmm , xmax = (Xsx<=Xmx)? Xsx : Xmx ;

	//Characteristics of histograms
	Double_t Axs = hsim->GetBinCenter(2)  - hsim->GetBinCenter(1)  ;
	Double_t Axm = hmeas->GetBinCenter(2) - hmeas->GetBinCenter(1) ;
	Int_t Nxs = TMath::Nint( (xmax-xmin)/Axs ) ;
	Int_t Nxm = TMath::Nint( (xmax-xmin)/Axm ) ;

	hsimc  = new TH1D( "hsimc" , hsim->GetTitle()  ,Nxs , xmin , xmax );
	hmeasc = new TH1D( "hmeasc", hmeas->GetTitle() ,Nxm , xmin , xmax );

	//Pack old histograms into new ones
	for ( Int_t i=1 ; i<=Nxs ; i++ ) {
		Double_t xval = hsimc->GetBinCenter( i ) ;
		Double_t val  = hsim->GetBinContent( hsim->FindBin( xval ) ) ;
		hsimc->SetBinContent( i , val );
	}

	for ( Int_t i=1 ; i<=Nxm ; i++ ) {
		Double_t xval = hmeasc->GetBinCenter( i ) ;
		Double_t val  = hmeas->GetBinContent( hmeas->FindBin( xval ) ) ;
		hmeasc->SetBinContent( i , val );
	}

	sameScale = 1 ;

}



/*------------------- PRIVATE METHODS -----------------------------*/
