/* 
   Usage:
   		edge_tree fnm [Cend]


   Give a measurement file and it will create a ROOT tree
   of the data. If Cend is given, then it calculates the RC deconvoluted
   signal and stores it, instead of the measured one.

   It uses 2 classes:
      TMeas     contains info specific of the measurement
      TWaveform  contains info relative to each waveform

   The objects of each class are:
      em:  edge_measurement is an object of TMeas
      wv:  an object of TWaveform. There are many waveforms stored within 1 em

   The tree created has 2 branches
      raw:   contains the raw data (em object)
      proc:  contains processed measurement information (wv object) 

   The script can be compiled outside of ROOT using:

         ./compile.sh edge_tree

   Geneve, May 2012
   Marcos Fernandez Garcia

   TO BE DONE:

   tree->SetUserInfo() and tree->GetUserInfo()

 */


#include <algorithm> //For the count method
#include <cstdarg>
#include <iostream>
#include <fstream>
#include <sstream> //ostringstream
#include <string>  //strcpy

#include <TFile.h>
#include <TTree.h>
#include <TMath.h>
#include <TH1.h>

#include <TMeas.h>
#include <TMeasHeader.h>
#include <TWaveform.h>
#include <TRACSInterface.h>

//ClassImp(TMeas)
//ClassImp(TMeasHeader)
//ClassImp(TWaveform)


/* ---------------------------------------------------------------- */

TTree * TRACSInterface::GetTree( ) {


	// Create a ROOT Tree
	TTree *tree = new TTree("edge","TRACS simulation");
	tree->SetDirectory(0);

	// Create a pointer to an raw data object
	TMeas *em = new TMeas( );

	em->volt = new Double_t [n_tSteps] ;
	em->time = new Double_t [n_tSteps] ;
	em->Qt   = new Double_t [n_tSteps] ;

	// Create branches
	tree->Branch("raw", &em,32000,1);

	//Read RAW file
	DumpToTree( em , tree ) ;


	delete em ;

	return tree;

}
/*----------------------------------------------------------*/

void TRACSInterface::DumpToTree( TMeas *em , TTree *tree ) {


	//Time information
	TString sdate = TString( line ) ;
	UShort_t  dd, mm, yy , hh , mn, ss;
	yy=GetYear();
	mm=GetMonth();
	dd=GetDay();
	hh=GetHour();
	mn=GetMinute();
	ss=GetSecond();
	TDatime date ;
	date.Set(yy, mm, dd, hh, mn, ss) ;
	em->utc = date ;

	//Frequency
	Int_t Freq = 200 ;

	//Total number of Scans

	Int_t NumOfScans = 1 ;
	if ( n_vSteps !=0 ) NumOfScans = NumOfScans * n_vSteps ;
	if ( n_zSteps !=0 ) NumOfScans = NumOfScans * n_zSteps ;
	if ( n_ySteps !=0 ) NumOfScans = NumOfScans * n_ySteps;


	// NVoltages
	em->NV = n_vSteps ;

	//Create a TMeasHeader object to store info that does not depend on event
	TMeasHeader *emh=new TMeasHeader( em->NV ) ;
	emh->Lambda = waveLength ;
	emh->NV = em->NV;
	emh->comment=TString( "TRACS simulated data" ) ;
	emh->Setup = 5 ;
	emh->Fluence = fluence ;
	emh->Nav  = 1 ;
	emh->Gain = 1.0 ;
	emh->iann = 0. ;
	emh->Illum = -2 ;

	//Cend
	emh->Cend = 0. ;

	//Vbias vector
	for ( int i=0 ; i< em->NV ; i++) emh->vVbias[i] = voltages[i] ;

	//Nominal power
	emh->Power = 1.0 ;

	//Ax
	emh->Ax = 0. ;
	emh->Nx = 1 ;

	//Ay
	em->Ay  = emh->Ay = deltaY  ;
	em->Ny  = emh->Ny = n_ySteps ;

	//Az
	em->Az  = emh->Az = deltaZ  ;
	em->Nz  = emh->Nz = n_zSteps ;

	em->Temp = temp ;

	int  iRead=0 , iactual = 1 ;
	int  Polarity = 0 ;
	Double_t ImA = 0., x0=0., y0=0., z0=0., t0s;
	em->x = 0. ;
	for (Int_t iloop = 0 ; iloop < NumOfScans ; iloop++ ) {

	    Int_t tms = iloop * 4000 ;
	    em->Vbias = detector->get_vbias();
        em->y = beamy ;
        em->z = beamz ;

		em->Itot = 0. ;

		//Calculate the bin slice along each coordinate
		if (iloop==0) { x0=em->x ;y0=em->y ;z0=em->z ; }
		em->ix = (em->Ax!=0.)? 1 + TMath::Nint((em->x - x0)/(em->Ax)):1 ;
		em->iy = (em->Ay!=0.)? 1 + TMath::Nint((em->y - y0)/(em->Ay)):1 ;
		em->iz = (em->Az!=0.)? 1 + TMath::Nint((em->z - z0)/(em->Az)):1 ;

		//Calculate the time increment wrt t0s
		if (iloop==0) t0s=tms ;
		Double_t Nseconds = 0.001*(tms - t0s) ;
		UShort_t ddi , hhi , mni, ssi , ndd, nhh, nmn , nss  ;
		ndd = TMath::Floor(Nseconds/86400) ;
		ddi = dd + ndd ;
		nhh = TMath::Floor((Nseconds - ndd*86400)/3600.) ;
		hhi = hh+nhh;
		nmn = TMath::Floor((Nseconds - ndd*86400-nhh*3600)/60) ;
		mni = mn+nmn;
		nss = TMath::Floor(Nseconds - ndd*86400-nhh*3600-nmn*60) ;
		ssi = ss + nss;
		if (ssi>60) { mni++ ; ssi=ssi-60; }
		if (mni>60) { hhi++ ; mni=mni-60; }
		if (hhi>24) { ddi++ ; hhi=hhi-24; }
		date.Set(yy, mm, ddi, hhi, mni, ssi) ;
		em->utc = date ;

		em->At = dt*1.e9 ; emh->At = em->At ;
		valarray I_tot = vectorObjetosSimuladosAlmacenadosCorrelativamente.Get_Itot[iloop] ;
		for ( int i=0 ; i< em->Nt ; i++) {
			em->volt[i] = I_tot[i] ;
			em->time[i] = i*em->At ;
			i_tot
		}


		//Estimate polarity
		if ( TMath::Abs(TMath::MaxElement(em->Nt,em->volt)) > TMath::Abs(TMath::MinElement(em->Nt,em->volt)) ) Polarity++ ;
		else Polarity-- ;

		em->event=iRead ;

		//Now postprocess this entry (find out baseline, rtime and so on
		TWaveform *wv = new TWaveform( em ) ;

		if (iRead==0) tree->Branch("proc" , &wv , 32000 , 1 );


		tree->Fill() ;
		iRead++   ;

		delete wv ;

		iactual++ ;

	}

	emh->Polarity = (Polarity>1) ? 1 : -1 ;
	tree->GetUserInfo()->Add( emh ) ;
	em->Ntevent=iRead ; //It does not go into the tree, only in the class!

	//delete emh ;


}


/*----------------------------------------------------------*/

int parse_TRACS( const char *filename, TMeas *em , TTree *tree , Double_t Cend ) {

	//Note that there is a difference in the coordinates. He calls X to my Y.
	//Open file
	ifstream in( filename );

	string line, what ;
	for (Int_t iloop = 1 ; iloop<= 3 ; iloop ++ ) getline(in, line);
	stringstream myStream(line);
	Double_t version ; myStream >> what >> version ;


	//StartTime
	Int_t istart=4 , iend = 6 ;
	for (Int_t iloop = istart ; iloop<= iend ; iloop ++ ) getline(in, line);
	TString sdate = TString( line ) ;
	UShort_t  dd, mm, yy , hh , mn, ss;
	yy=atoi( sdate(11,4).Data() );
	mm=atoi( sdate(16,2).Data() );
	dd=atoi( sdate(19,2).Data() );
	hh=atoi( sdate(22,2).Data() );
	mn=atoi( sdate(25,2).Data() );
	ss=atoi( sdate(28,2).Data() );
	TDatime date ;
	date.Set(yy, mm, dd, hh, mn, ss) ;
	em->utc = date ;


	//Wavelength
	getline(in, line);
	myStream.str(""); myStream.clear() ; myStream << line ;
	Double_t wvlength ; myStream >> what >>  wvlength ;

	//Top/Bottom/edge
	//Try to decide on the type of scan this is 0=tct std, 1=eTCT (default), 2=TCT+
	Int_t Illum ;
	istart=iend+1; //11
	iend=istart;   //11
	for (Int_t iloop = istart ; iloop<= iend ; iloop ++ ) getline(in, line);
	myStream.str(""); myStream.clear() ; myStream << line ;
	string direction ; myStream >> what >>  direction ;
	if      ( direction.find("edge") == 1  ) Illum = 1 ;
	if      ( direction.find("edge") == 0  ) Illum = 0 ;
	if      ( direction.find("edge") == -1 ) Illum = -1 ;

	//Amp gain
	Double_t Gain ;
	getline(in, line);
	myStream.str(""); myStream.clear() ; myStream << line ;
	myStream >> what >> Gain ;


	//Fluence and annealing
	Double_t Phi =0., iann =0. ;

	//Fluence
	getline(in, line);
	myStream.str(""); myStream.clear() ; myStream << line ;
	myStream >> what >> Phi ;

	//Annealing
	getline(in, line);
	myStream.str(""); myStream.clear() ; myStream << line ;
	myStream >> what >> iann ;


	//Total number of Scans
	getline(in, line);
	myStream.str(""); myStream.clear() ; myStream << line ;
	Int_t NumOfScans ; myStream >> what >> NumOfScans ;


	// NVoltages
	Int_t Ival , t0s , tms ;
	for (Int_t iloop = 1 ; iloop<= 9 ; iloop ++ ) getline(in, line);
	myStream.str(""); myStream.clear() ; myStream << line ;
	myStream >> what ; myStream >> Ival ;
	em->NV = Ival ;

	//Create a TMeasHeader object to store info that does not depend on event
	TMeasHeader *emh=new TMeasHeader( em->NV ) ;
	emh->Lambda = wvlength ;
	emh->NV = em->NV;
	emh->comment="" ;
	emh->Setup=5 ;
	emh->Fluence = Phi ;
	emh->Nav   = 1 ;
	emh->Gain  = Gain ;
	emh->iann  = iann ;
	emh->Illum = Illum ;


	//Vbias vector
	getline(in, line); myStream.str(""); myStream.clear() ; myStream << line ;
	myStream >> what ; for ( int i=0 ; i< em->NV ; i++) myStream >> emh->vVbias[i] ;

	//Ax
	istart = 1;
	iend   = istart+12;
	Double_t Dval ;
	for (Int_t iloop = istart ; iloop<= iend ; iloop ++ ) getline(in, line); myStream.str(""); myStream.clear() ; myStream << line ;
	myStream >> what >> Dval ;
	em->Ax  = Dval *1000. ;
	emh->Ax = em->Ax;

	//Nx
	istart = iend+1;
	iend   = istart;
	getline(in, line); myStream.str(""); myStream.clear() ; myStream << line ;
	myStream >> what >> em->Nx ;
	emh->Nx = em->Nx;

	//Ay
	istart = iend+1;
	iend   = istart+2;
	for (Int_t iloop = istart ; iloop<= iend ; iloop ++ ) getline(in, line); myStream.str(""); myStream.clear() ; myStream << line ;
	myStream>> what  >> Dval ;
	em->Ay  = Dval *1000. ;
	emh->Ay = em->Ay;

	//Ny
	istart = iend+1;
	iend   = istart;
	getline(in, line); myStream.str(""); myStream.clear() ; myStream << line ;
	myStream >> what >> em->Ny ;
	emh->Ny = em->Ny;

	//Az
	istart = iend+1;
	iend   = istart+3;
	for (Int_t iloop = istart ; iloop<= iend ; iloop ++ ) getline(in, line); myStream.str(""); myStream.clear() ; myStream << line ;
	myStream>> what  >> Dval ;
	em->Az  = Dval *1000. ;
	emh->Az = em->Az;

	//Nz
	istart = iend+1;
	iend   = istart;
	getline(in, line); myStream.str(""); myStream.clear() ; myStream << line ;
	myStream >> what >> em->Nz ;
	emh->Nz = em->Nz;

	//At
	getline(in, line); getline(in, line); myStream.str(""); myStream.clear() ; myStream << line ;
	myStream >> what >> em->At ;
	em->At  = em->At*1.e9;
	emh->At = em->At ;

	//Capacitance
	getline(in, line); myStream.str(""); myStream.clear() ; myStream << line ;
	myStream >> what >> Cend ;
	emh->Cend = Cend*1.e12 ;

	//remaining info (useful for logging, not in the trees)
	for (Int_t iloop = 1 ; iloop<= 7 ; iloop ++ ) getline(in, line) ;

	if ( emh->Illum == 1  ) cout << "TRACS: top illumination"   << endl ;
	if ( emh->Illum == 0  ) cout << "TRACS: edgeTCT" << endl ;
	if ( emh->Illum == -1 ) cout << "TRACS: bottom illumination"  << endl ;

	for (Int_t iloop = 1 ; iloop<= 4 ; iloop ++ ) getline(in, line);

	int  iRead=0 , iactual = 1 ;
	int  Polarity = 0 ;
	Double_t ImA,x0,y0,z0;
	for (Int_t iloop = 0 ; iloop < NumOfScans ; iloop++ ) {

		getline(in, line); myStream.str(""); myStream.clear() ; myStream << line ;
		if ( in.eof() ) break;
		myStream >> em->Nt >> em->Temp >> em->Vbias >> em->x >> em->y >> em->z  ;

		//Calculate the bin slice along each coordinate
		if (iloop==0) { x0=em->x ;y0=em->y ;z0=em->z ; }
		em->ix = (em->Ax!=0.)? 1 + TMath::Nint((em->x - x0)/em->Ax):1 ;
		em->iy = (em->Ay!=0.)? 1 + TMath::Nint((em->y - y0)/em->Ay):1 ;
		em->iz = (em->Az!=0.)? 1 + TMath::Nint((em->z - z0)/em->Az):1 ;

		for ( int i=0 ; i< em->Nt ; i++) {
			myStream >> em->volt[i] ;
			em->time[i]=i*em->At ;
		}

		//Estimate polarity
		if ( TMath::Abs(TMath::MaxElement(em->Nt,em->volt)) > TMath::Abs(TMath::MinElement(em->Nt,em->volt)) ) Polarity++ ;
		else Polarity-- ;

		em->event=iRead ;

		//Now postprocess this entry (find out baseline, rtime and so on
		TWaveform *wv = new TWaveform( em ) ;

		if (iRead==0) tree->Branch("proc" , &wv , 32000 , 1 );

		tree->Fill() ;
		iRead++   ;
		//if (iRead%1000 == 0) tree->AutoSave();

		if (iRead%100==0) cout << "Read " << iRead << flush <<"\r" ;
		delete wv ;

		iactual++ ;

	}

	emh->Polarity = (Polarity>1) ? 1 : -1 ;
	tree->GetUserInfo()->Add( emh ) ;
	cout << endl ;
	cout << "Total read:" << iRead << endl ;
	em->Ntevent=iRead ; //It does not go into the tree, only in the class!

	in.close() ;

	//delete emh ;

	return(1) ;

}
