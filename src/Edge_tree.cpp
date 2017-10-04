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

#include "TMeas.h"
#include "TMeasHeader.h"
#include "TWaveform.h"

#include <cstdarg>
#include <iostream>
#include <fstream>
#include <sstream> //ostringstream
#include <string>  //strcpy

#include "TFile.h"
#include "TTree.h"
#include "TMath.h"
#include "TH1.h"

#include <algorithm> //For the count method

#define DEBUG  0
#define NXCHAR 512
#define ROSC   50

#define SKIPLINES 0    //if not 0, then skip these many lines in each scan. If SKIPLINES=1, then it skips 1 out of 2
		       //lines of input, which is equivalent to setting Az=2. So the new AZ=SKIPLINES+1

#define SMOOTHING 0

#define READONLYPDWVF 0      //1= skip sensor waveform and read only the photodiode waveform (tct+ only)

#define AWIR  0.77 //A/W  

#define AWRED 0.4  //A/W  

#define JeV      1.60217657e-19  //J
#define QELEC    1.60217646e-19  //C

#define BL0SUBS    0 //Correct baseline by substracting measurement at Vbias=0
#define XYTEMP_INFO 0   //0=IFCA data does not have any coordinate  info
                        //1=In case we add (X,Y,temp) info to IFCA TCT
			//2=Adding (X,Y,Z,Temp) to IFCA data

#define SPACORR 0  //Only for TPA files. Correct all waveforms by substracting
		   //1st waveform Vi(t)=[Vi(t)-Pi/P1*V1(t)]/Pi**2

#define SPAEVENT 0      //Number of event (starting by 0) taken as reference for SPA correction

#define DISCARDEVENT -1  //We can skip 1 event
                         //-1 Discard no event
			 //>=0 discard event

ClassImp(TMeas)
ClassImp(TMeasHeader)
ClassImp(TWaveform)

// P R O T O T Y P I N G .....

int       parse_file( const char *filename, TMeas *em , TTree *tree , Double_t Cend  ) ;
int    parse_tctplus( const char *filename, TMeas *em , TTree *tree , Double_t Cend  ) ;
int        parse_TPA( const char *filename, TMeas *em , TTree *tree , Double_t Cend  ) ;
int     parse_tctUHH( const char *filename, TMeas *em , TTree *tree , Double_t Cend  ) ;
int       parse_ifca( const char *filename, TMeas *em , TTree *tree , Double_t Cend  ) ;
int       parse_TRACS( const char *filename, TMeas *em , TTree *tree , Double_t Cend  ) ;
int       parse_IDLTS( const char *filename, TMeas *em , TTree *tree , Double_t Cend  ) ;
int    ReadNerOfBins( const char *filename ) ;

//rootcint mytest_dictionary.cxx -c mytest.h LinkDef.h

/* ---------------------------------------------------------------- */
/**
 *
 * @param argc
 * @return
 */
int tree( int argc ,  ... ) {
  
  va_list argv;
  
  // Initializing argv to store all values after argc
  va_start( argv, argc );  
    
  char *pfnm = va_arg( argv , char *) ;
  Double_t Cend = 0. ;
  if ( argc==2 ) Cend = atof( va_arg( argv , char *) ) ;
  Cend = Cend / 1.e3;  //Input in pF, then changing to nF, so RC=ns

  /* Build root filename */
  string fnm , bname=pfnm , ext=".root"; 
  
  //Pick original extension (ZVscan, YZscan, ...)
  int idot    = bname.find(".") ;
  string sext = bname.substr(idot) ;
  bname       = bname.substr(0,idot) ;
  
  if (Cend!=0) {   
    ostringstream sCend ;
    Int_t icend = TMath::Nint(Cend*1.e3) ;
    sCend << icend <<"pF" ;
    bname =  bname+ "_deconv" + sCend.str() ;
  }
  
  if (SKIPLINES!=0) {
    ostringstream sskip ;
    Int_t iskip = SKIPLINES ;
    sskip << iskip ;   
    bname = bname+ "_skip" + sskip.str() ;
  }
  
  fnm = bname + sext + ext ;
  cout<<"Output file: " << fnm.c_str() << endl ;
  
  
  //create a Tree file tree4.root
  TFile froot( fnm.c_str() , "RECREATE" );

  // Create a ROOT Tree
  TTree *tree = new TTree("edge","eTCT measurement");

  // Create a pointer to an raw data object
  TMeas *em = new TMeas( );
  em->Nt=ReadNerOfBins( pfnm ) ;

  em->volt = new Double_t [em->Nt] ;
  em->time = new Double_t [em->Nt] ;
  em->Qt   = new Double_t [em->Nt] ;

  // Create branches
  tree->Branch("raw", &em,32000,1);

  //Read RAW file
  parse_file( pfnm , em , tree , Cend ) ;

  // Write the file header
  froot.Write();
  delete tree;

  froot.Close();

  delete em ;
    
  va_end ( argv );   // Cleans up the list
  
  return 0;
  
}
/*----------------------------------------------------------*/
/**
 *
 * @param filename
 * @param em
 * @param tree
 * @param Cend
 * @return
 */
int parse_file( const char *filename, TMeas *em , TTree *tree , Double_t Cend ) {


    //Open file
    ifstream in( filename );
    if (!in || in.bad()) {
      cout << "Error opening " << filename << endl ;
      return(0);
    }
    
    cout << "Creating tree for " << filename << endl ;
    

    //Read header
    string what , sval , unit , line ;
    
    
    getline(in,line) ; if (DEBUG) cout << line.c_str() <<endl ; //Skip line
    if ( line.find("================") == 0 ) {
      getline(in,line);
      if ( line.find("SSD simulation") == 0 ) {      
	in.close() ;
	em->Setup = 5;
	parse_TRACS( filename , em , tree , Cend ) ; 
	return(1) ;
     } else if ( line.find("SSD measurement") == 0 ) { 
        for (Int_t i=0;i<3;i++)getline(in,line) ; 
	if ( line.find("scanType: TCT+") == 0 ) {
	  in.close() ;
	  em->Setup = 2;
	  parse_tctplus( filename , em , tree , Cend ) ; 
	} else if ( line.find("scanType: DLTS") == 0 ) {
	  in.close() ;
	  em->Setup = 6;
	  parse_IDLTS( filename , em , tree , Cend ) ; 
        }
	return(1) ;
      }
    }
    if ( line.find("----------------") == 0 ) {
      in.close() ;
      em->Setup = 3;
      parse_TPA( filename , em , tree , Cend ) ; 
      return(1) ;
    }
    if ( line.find("MTCT Header") == 0 ) {
      in.close() ;
      em->Setup = 3;
      parse_tctUHH( filename , em , tree , Cend ) ; 
      return(1) ;
    }
    if ( line.find("* t0      :") == 0 ) {
      in.close() ;
      em->Setup = 4;
      parse_ifca( filename , em , tree , Cend ) ; 
      return(1) ;
    }
    
    //This is eTCT type of file
    em->Setup = 1;
    
    getline(in,line) ; if (DEBUG) cout << line.c_str() <<endl ; //Skip line
    getline(in,line) ; if (DEBUG) cout << line.c_str() <<endl ; //Skip line


    // NVoltages:      1
    in >> what >> em->NV ;  
    em->NV = ( em->NV==0 ) ? 1 : em->NV ;  //Lvb bug
    if (DEBUG) cout << what.c_str() << em->NV <<endl ;

    //Create a TMeasHeader object to store info that does not depend on event
    TMeasHeader *emh=new TMeasHeader( em->NV ) ;
    emh->NV = em->NV;
    emh->Setup = 1;
    
    //Cend
    emh->Cend = Cend ;
    
    // NXpos:  0
    in  >> what >> em->Nx ; 
    if (DEBUG) cout<<what.c_str() << em->Nx <<endl ;
    emh->Nx = em->Nx;
    
    // X0:     -1.700000E+1 mm
    in >> what >> sval >> unit  ;  //Skip line
  
    // dX:     0.000000E+0 mm
    in  >> what >> em->Ax >> unit;  
    if (DEBUG) cout << what.c_str() << em->Ax <<endl ; 
    em->Ax=em->Ax*1000. ;
    emh->Ax = em->Ax;
    
    // NYpos:  0
    in  >> what >> em->Ny ;  
    if (DEBUG) cout << what.c_str()<< em->Ny <<endl ;
    emh->Ny = em->Ny;
    
    // Y0:     -7.300000E+0 mm
    in >> what >> sval >> unit  ;  //Skip line
    
    // dY:     0.000000E+0 mm
    in  >> what >> em->Ay >> unit;  
    em->Ay=em->Ay*1000. ;
    if (DEBUG) cout<<what.c_str() << em->Ay <<endl ;
    emh->Ay = em->Ay;

    // NZpos:  401
    in  >> what >> em->Nz;
    if (DEBUG) cout << what.c_str() << em->Nz<<endl ;
    #if SKIPLINES!=0
      Int_t Nzorig = em->Nz ;
      em->Nz = ( (double)em->Nz/2.0 - (int)(em->Nz/2.0) == 0.5 )? 
               1+em->Nz/(SKIPLINES+1) : em->Nz/(SKIPLINES+1) ;
    #endif
    emh->Nz = em->Nz;
    
    // Z0:     7.830000E+0 mm
    in >> what >> sval >> unit  ;  //Skip line

    // dZ:     2.000000E-3 mm
    in  >> what >> em->Az >> unit;  
      em->Az = em->Az*1000. ;
    #if SKIPLINES!=0
      em->Az = (SKIPLINES+1.0)*em->Az ;         //Already in mum
    #endif
    if (DEBUG) cout << what.c_str() << em->Az<<endl ;
    emh->Az = em->Az;
    
    
    // NTimePoints:    10003
    in  >> what >> em->Nt;	     //14 NTimePoints:  10003
    if (DEBUG) cout << what.c_str() << em->Nt<<endl ;

    // t0:     -1.611310E-8 s
    in >> what >> sval >> unit  ;  //Skip line

    // dt:     5.000000E-11 s
    in  >> what >> em->At >> unit;  //16 dt:	5.000000E-11 s
    em->At=em->At*1.e9 ;
    if (DEBUG) cout << what.c_str() << em->At<<endl ;
    emh->At = em->At ;

    //Comment
    string comm ;
    in.ignore(); getline(in,comm)  ;  //Skip line
    if (DEBUG) cout << comm << endl ;
    emh->comment=comm ;
        
    //HEADER-END
    string hend ;
    getline(in,hend) ;
    if (DEBUG) cout << hend << endl  ;  //Skip line
     
    UShort_t  dd, mm, yy , hh , mn, ss;
    TDatime date ;

    char   c ;
    int    iRead=0 , iactual = 1 ;
    int    ivb= 0 , Polarity = 0 ;
    Double_t VbiasOld=-999999.9 ;
    Double_t x0 , y0 , z0 ;
    in >> what >> dd >> c>> mm >> c >> yy>> hh >> c >> mn >>c >> ss ;
    while (!in.eof() && !in.bad()) {
      
      //Date
      if (DEBUG) cout <<what<<" "<<dd<<c<< mm<<c<< yy<<" "<< hh<<c << mn<<c << ss ;
      date.Set(yy, mm, dd, hh, mn, ss) ;
      em->utc = date ;
      
      //Itot , Vbias , X,Y,Z
      in >> em->Itot >> em->Vbias >> em->x >> em->y >> em->z  ;
      if (DEBUG) cout<<" " << em->Itot<<" " << em->Vbias<<" " << em->x<<" " << em->y<<" " << em->z  ;
            
      //Calculate the bin slice along each coordinate
      if (iRead==0) { x0=em->x ;y0=em->y ;z0=em->z ; }
      em->ix = (em->Ax!=0.)? 1 + TMath::Nint((em->x - x0)/em->Ax):1 ;
      em->iy = (em->Ay!=0.)? 1 + TMath::Nint((em->y - y0)/em->Ay):1 ;
      em->iz = (em->Az!=0.)? 1 + TMath::Nint((em->z - z0)/em->Az):1 ;
      
      if (em->Vbias!=VbiasOld) {
        emh->vVbias[ivb] = em->Vbias ;
	VbiasOld = em->Vbias ;
	ivb++;
      }
      
      //All Voltages
       for ( int i=0 ; i< em->Nt ; i++) {
        in >> em->volt[i] ;
        em->time[i]=i*em->At ; 
      }
      
      if ( Cend ) {
        
	/* Deconvolute the pulse */
        Double_t *dvolt = new Double_t [em->Nt] ;
        for ( int i=0 ; i< em->Nt-1 ; i++) dvolt[i]    = ROSC*Cend*(em->volt[i+1] - em->volt[i])/em->At + em->volt[i] ;
	dvolt[em->Nt-1] = dvolt[em->Nt-2] ;
	
	#if SMOOTHING==1
	  Double_t tmin = em->time[0] ;
	  Double_t tmax = em->time[em->Nt-1] ;
	  Double_t At   = (tmax-tmin)/(em->Nt-1.0) ;
	  tmin = -0.5*At ;
	  tmax = tmin + em->Nt*At ;
	  TH1D *hvt=new TH1D( "hvt" , "deconv" , em->Nt , tmin , tmax ) ;
	  hvt->FillN(em->Nt, em->time , dvolt ) ;

	  hvt->Smooth( 4000 );
          cout << "Smoothing waveform " << iRead << flush <<"\r" ;

          for ( int i=0 ; i< em->Nt ; i++) em->volt[i] = hvt->GetBinContent(i) ;
          delete hvt;
          hvt=0;
        # else
	  for ( int i=0 ; i< em->Nt ; i++) em->volt[i] = dvolt[i] ;
 	#endif
	delete [] dvolt ; 
	dvolt=0;
      }
      
      //Estimate polarity
      if ( TMath::Abs(TMath::MaxElement(em->Nt,em->volt)) > TMath::Abs(TMath::MinElement(em->Nt,em->volt)) ) Polarity++ ;
      else Polarity-- ;

      //DATA-END-x
      in >> what ;
      if (DEBUG) cout <<" " << what<<endl ;

      em->event=iRead ;
      
      //Now postprocess this entry (find out baseline, rtime and so on
      TWaveform *wv = new TWaveform( em ) ;
	
      if (iRead==0) tree->Branch("proc" , &wv , 32000 , 1 );
       
      tree->Fill() ;
      iRead++   ;
      //if (iRead%1000 == 0) tree->AutoSave();

      if (iRead%100==0) cout << "Read " << iRead << flush <<"\r" ;
      delete wv ; 

      //Get the header of the event, or give an error, so the while will now stop
      in >> what >> dd >> c>> mm >> c >> yy>> hh >> c >> mn >>c >> ss ;
      iactual++ ;
      
      //Check if we want to read 1 out of SKIPLINES
      #if SKIPLINES!=0
	for (Int_t iskip=0 ; iskip<SKIPLINES ; iskip++ ) { 
          if ( iactual % Nzorig != 1 ) {   //We do not want to skip the first line of each Vbias scan
	    in.ignore(); getline(in,line);  //Skip one line, then read header again
	    //cout << "Skipping line "<< iactual << " ," <<hh << ":" << mn <<":"<< ss << endl;
	    in >> what >> dd >> c>> mm >> c >> yy>> hh >> c >> mn >>c >> ss ;
            iactual++ ;
          } 
	}
      #endif
    } 
    
    emh->Polarity = (Polarity>1) ? 1 : -1 ;
    tree->GetUserInfo()->Add( emh ) ;
    cout << endl ;
    cout << "Total read:" << iRead << endl ;
    em->Ntevent=iRead ; //It does not go into the tree, only in the class!
    
    in.close() ;
        
    //delete emh ;
    
    return(1);

}
/*----------------------------------------------------------*/
/**
 *
 * @param filename
 * @param em
 * @param tree
 * @param Cend
 * @return
 */
int parse_ifca( const char *filename, TMeas *em , TTree *tree , Double_t Cend ) {


    //Open file
    ifstream in( filename );
    if (!in || in.bad()) {
      cout << "Error opening " << filename << endl ;
      return(0);
    }
    
    cout << "I believe this is an IFCA measurement"   << endl ;    
    
    
    Double_t dt ;
    Int_t    Nt , NV ;

    //* t0      : 3.120000E-007 *
    string what , sval , unit , line ;
    getline(in,line) ; //Skip line

    //* dt      : 1.000000E-011 *
    getline(in,line) ;  stringstream myStream ; myStream<<line ;
    myStream  >> what >> what >> what >> dt >> what ;  
    em->At=dt*1.e9 ;
    
    //* Samples : 5000 *
    getline(in,line) ; myStream.str(""); myStream.clear() ; myStream << line ;
    cout << myStream.str() << endl ;
    myStream  >> what >> what >> what >> Nt >> what ;
    em->Nt = Nt ;

    em->volt = new Double_t [em->Nt] ;
    em->time = new Double_t [em->Nt] ;

    #if  XYTEMP_INFO==1
      //Read in header info related to motor steps. Caution, z coordinate=temp
      Double_t Ax, Ay ;
      UShort_t Nx, Ny ;
      getline(in,line) ; myStream.str(""); myStream.clear() ; myStream << line ;
      myStream  >> what >> Nx ;  
      getline(in,line) ; myStream.str(""); myStream.clear() ; myStream << line ;
      myStream  >> what >> Ax ;  
      getline(in,line) ; myStream.str(""); myStream.clear() ; myStream << line ;
      myStream  >> what >> Ny ;  
      getline(in,line) ; myStream.str(""); myStream.clear() ; myStream << line ;
      myStream  >> what >> Ay ;  
      getline(in,line) ; myStream.str(""); myStream.clear() ; myStream << line ;
    #endif
    #if  XYTEMP_INFO==2
      //Read in header info related to motor steps. Caution, z coordinate=temp
      Double_t Ax, Ay , Az ;
      UShort_t Nx, Ny , Nz ;
      getline(in,line) ; myStream.str(""); myStream.clear() ; myStream << line ;
      myStream  >> what >> Nx ;  
      getline(in,line) ; myStream.str(""); myStream.clear() ; myStream << line ;
      myStream  >> what >> Ax ;  
      getline(in,line) ; myStream.str(""); myStream.clear() ; myStream << line ;
      myStream  >> what >> Ny ;  
      getline(in,line) ; myStream.str(""); myStream.clear() ; myStream << line ;
      myStream  >> what >> Ay ;  
      getline(in,line) ; myStream.str(""); myStream.clear() ; myStream << line ;
      myStream  >> what >> Nz ;  
      getline(in,line) ; myStream.str(""); myStream.clear() ; myStream << line ;
      myStream  >> what >> Az ;  
    #endif
      
    //Count ner of voltages
    getline(in,line) ; myStream.str(""); myStream.clear(); myStream << line ;
    NV=0;
    Double_t val ;
    while ( myStream>>val ) { 
      NV++;
      if (myStream.peek() == '\n') break;  //Detects EOL
    }
    em->NV = NV ;

    TMeasHeader *emh=new TMeasHeader( em->NV ) ;
    emh->NV = em->NV;

    //Cend
    emh->Cend = Cend ;
    
    //Go back to the line with voltages
    Double_t *Vbias = new Double_t [em->NV] ;
    myStream.str(""); myStream.clear(); myStream << line ;
    Int_t iv = 0 ;
    while ( myStream>>val ) { 
       //Vbias[iv]=-1.*val ;if (iv==0) cout<<"Multiplying Vbias by -1 (negative signal)"<<endl ;
       Vbias[iv]=val ;
       emh->vVbias[iv] = Vbias[iv] ;
       iv++;
    }

    #if  XYTEMP_INFO==1

      emh->Nx = Nx ;
      emh->Ny = Ny ;
      emh->Nz = NV ;  
      emh->Ax = Ax ;
      emh->Ay = Ay ;

      //Read (x,y,temp) [I WILL PUT temperature in z coordinate!!!!!!!!]
      Char_t c ;
      Double_t valx, valy, valz ;
      Double_t *vx    = new Double_t [em->NV] ;
      Double_t *vy    = new Double_t [em->NV] ;
      Double_t *vz    = new Double_t [em->NV] ;   //temperature
      getline(in,line) ; myStream.str(""); myStream.clear(); myStream << line ;
      iv = 0 ;
      while ( myStream >> c >> valx >> c >> valy >> c >> valz >> c  ) { 
         vx[iv] = valx ;
         vy[iv] = valy ;
	 vz[iv] = valz ;
         iv++;
      }
    #endif   
    #if  XYTEMP_INFO==2

      emh->Nx = Nx ;
      emh->Ny = Ny ;
      emh->Nz = Nz ;  
      emh->Ax = Ax ;
      emh->Ay = Ay ;
      emh->Az = Az ;

      //Read (x,y,z,temp)
      Char_t c ;
      Double_t valx, valy, valz , valt ;
      Double_t *vx    = new Double_t [em->NV] ;
      Double_t *vy    = new Double_t [em->NV] ;
      Double_t *vz    = new Double_t [em->NV] ;   
      Double_t *vt    = new Double_t [em->NV] ;   //temperature
      getline(in,line) ; myStream.str(""); myStream.clear(); myStream << line ;
      iv = 0 ;
      while ( myStream >> c >> valx >> c >> valy >> c >> valz >> c >> valt >> c ) { 
         vx[iv] = valx ;
         vy[iv] = valy ;
	 vz[iv] = valz ;
	 vt[iv] = valt ;
         iv++;
      }
    #endif   

    /* NOW READ THE VOLTAGES IN COLUMNS!!!! */
    long FirstRow = in.tellg();
    
    int    iRead=0 , Polarity = 0 ;
    Double_t x0 , y0 , z0 ;
    for (Int_t CurrCol=0 ; CurrCol < NV ; CurrCol++) {
      
      em->Vbias  = Vbias[CurrCol] ;
      #if  XYTEMP_INFO==1
	em->x      = vx[CurrCol] ;
	em->y      = vy[CurrCol] ;
	em->z      = vz[CurrCol] ;   //temperature!!!
      #endif
      #if  XYTEMP_INFO==2
	em->x      = vx[CurrCol] ;
	em->y      = vy[CurrCol] ;
	em->z      = vz[CurrCol] ;   
	em->Temp   = vt[CurrCol] ;
      #endif
            
      //Calculate the bin slice along each coordinate
      if (CurrCol==0) { x0=em->x ;y0=em->y ;z0=em->z ; }
      em->ix = (em->Ax!=0.)? 1 + TMath::Nint((em->x - x0)/em->Ax):1 ;
      em->iy = (em->Ay!=0.)? 1 + TMath::Nint((em->y - y0)/em->Ay):1 ;
      em->iz = (em->Az!=0.)? 1 + TMath::Nint((em->z - z0)/em->Az):1 ;
      
      //Read one full column (~v/read)
      for ( Int_t iloop = 0 ; iloop<Nt; iloop++ ) {

         getline(in, line); myStream.str("") ; myStream.clear() ;
	 myStream << line ;

         //Go to the proper column
         Double_t val ;
         for (Int_t icol=0 ; icol < CurrCol ; icol++) myStream >> val ;

         //Read in the value
         myStream >> em->volt[iloop] ;
         em->time[iloop] = iloop*em->At ; 
                 
         //Discard the rest of the line. To mix >> with getline, one needs ignore after >>.
         //in.ignore() ; 
         //getline( in , what ) ;
         //cout << iloop <<" " << em->volt[iloop]<<endl;

      }
      
      //Check if we want deconvolution 
      if ( Cend ) {

	/* Deconvolute the pulse */
	Double_t *dvolt = new Double_t [em->Nt] ;
	for ( int i=0 ; i< em->Nt-1 ; i++) dvolt[i] = ROSC*Cend*(em->volt[i+1] - em->volt[i])/em->At + em->volt[i] ;
	dvolt[em->Nt-1] = dvolt[em->Nt-2] ;

	#if SMOOTHING==1
	  Double_t tmin = em->time[0] ;
	  Double_t tmax = em->time[em->Nt-1] ;
	  Double_t At   = (tmax-tmin)/(em->Nt-1.0) ;
	  tmin = -0.5*At ;
	  tmax = tmin + em->Nt*At ;
	  TH1D *hvt=new TH1D( "hvt" , "deconv" , em->Nt , tmin , tmax ) ;
	  hvt->FillN(em->Nt, em->time , dvolt ) ;

	  hvt->Smooth( 4000 );
          cout << "Smoothing waveform " << iRead << flush <<"\r" ;

          for ( int i=0 ; i< em->Nt ; i++) em->volt[i] = hvt->GetBinContent(i) ;
          delete hvt;
          hvt=0;
	# else
	  for ( int i=0 ; i< em->Nt ; i++) em->volt[i] = dvolt[i] ;
	#endif
	delete [] dvolt ; 
	dvolt=0;
      }
      
      em->event=iRead ;
            
      //Now postprocess this entry (find out baseline, rtime and so on
      TWaveform *wv = new TWaveform( em ) ;
      if (iRead==0) tree->Branch("proc" , &wv , 32000 , 1 );
      
      tree->Fill() ;
      iRead++ ;
      
      //Estimate polarity
      if ( TMath::Abs(TMath::MaxElement(em->Nt,em->volt)) > TMath::Abs(TMath::MinElement(em->Nt,em->volt)) ) Polarity++ ;
      else Polarity-- ;


      cout << "Read " << iRead << flush <<"\r" ;
      delete wv ; 

      if ( CurrCol<NV ) in.seekg(FirstRow,ios::beg);


    }
        
    emh->Polarity = (Polarity>1) ? 1 : -1 ;
    tree->GetUserInfo()->Add( emh ) ;
    cout << endl ;
    cout << "Total read:" << iRead << endl ;
    em->Ntevent=iRead ; //It does not go into the tree, only in the class!
    
    in.close() ;
    delete Vbias ;

    #if XYTEMP_INFO==1
      delete vx ;
      delete vy ;
      delete vz ;
    #endif
    #if XYTEMP_INFO==2
      delete vx ; delete vy ; delete vz ; delete vt ;
    #endif
        
    //delete emh ;
    
    return(1);

}

/*----------------------------------------------------------*/
/**
 *
 * @param filename
 * @param em
 * @param tree
 * @param Cend
 * @return
 */
int parse_tctplus( const char *filename, TMeas *em , TTree *tree , Double_t Cend ) {
    
    //Note that there is a difference in the coordinates. He calls X to my Y.
    //Open file
    ifstream in( filename );    

    string line, what ;
    for (Int_t iloop = 1 ; iloop<= 3 ; iloop ++ ) getline(in, line);
    stringstream myStream(line);
    Double_t version ; myStream >> what >> version ;

    //Maybe it would be better to read by keywords instead of by line position   
    Int_t offset = 0  ;
    if (version>=1.1) offset=1;
    
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

    //comment
    istart=iend+1; //7
    iend=istart+1; //8
    for (Int_t iloop = istart ; iloop<= iend ; iloop ++ ) getline(in, line);
    TString comment = TString( line ) ;
    
    //Wavelength
    istart=iend+1; //9
    iend=istart+1; //10
    for (Int_t iloop = istart ; iloop<= iend ; iloop ++ ) getline(in, line); 
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
    if ( direction.find("top")    ==0 ) Illum=1;
    if ( direction.find("edge")   ==0 ) Illum=0;
    if ( direction.find("bottom") ==0 ) Illum=-1;
    
    //Amp gain
    Double_t Gain ;
    getline(in, line);  
    myStream.str(""); myStream.clear() ; myStream << line ;
    myStream >> what >> Gain ;


    //Fluence and annealing
    Double_t Phi =0., iann =0. ;
    if (version <=1.2 ) getline(in, line); //read and discard biasteeresistance
    if (version ==1.3 ) {
      getline(in, line); //read and discard biasteeresistance
      
      //Fluence
      getline(in, line); 
      myStream.str(""); myStream.clear() ; myStream << line ;
      myStream >> what >> Phi ;
      
      //Annealing
      getline(in, line); 
      myStream.str(""); myStream.clear() ; myStream << line ;
      myStream >> what >> iann ;
      
    }
    
    //Frequency
    getline(in, line); 
    myStream.str(""); myStream.clear() ; myStream << line ;
    Int_t Freq ; myStream >> what >> Freq ;
    
    //Number of averages
    getline(in, line); 
    myStream.str(""); myStream.clear() ; myStream << line ;
    Int_t Nav ; myStream >> what >> Nav ;
    
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
    emh->comment=comment ;
    emh->Setup = 2 ;
    emh->Fluence = Phi ;
    emh->Nav  = Nav ;
    emh->Gain = Gain ;
    emh->iann = iann ;
    emh->Illum = Illum ;
   
    //Cend
    emh->Cend = Cend ;
    
    //Vbias vector
    getline(in, line); myStream.str(""); myStream.clear() ; myStream << line ;
    myStream >> what ; for ( int i=0 ; i< em->NV ; i++) myStream >> emh->vVbias[i] ;
    
    //Nominal power
    Int_t nPower ; Double_t Power ;
    for (Int_t iloop = 1 ; iloop<= 4 ; iloop ++ ) getline(in, line); 
    myStream.str(""); myStream.clear() ; myStream << line ;
    myStream >> what ; myStream >> nPower ;
    getline(in, line); 
    myStream.str(""); myStream.clear() ; myStream << line ;
    myStream >> what ; myStream >> Power ;
    emh->Power = Power ;


    //Ax
    istart = 1; 
    iend   = istart+8; 
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
    iend   = istart+3; 
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
    
    if ( emh->Illum == 0 ) cout << "TCT+: this is an eTCT measurement" << endl ;
    if ( emh->Illum != 0 ) cout << "TCT+: this is a normal TCT measurement" << endl ;
    
    istart = iend+1; 
    iend   = istart+3; 
    for (Int_t iloop = istart ; iloop<= iend ; iloop ++ ) getline(in, line); 
    
    int  iRead=0 , iactual = 1 ;
    int  Polarity = 0 ;
    Double_t ImA, x0, y0, z0;
    for (Int_t iloop = 0 ; iloop < NumOfScans ; iloop++ ) {
      
      getline(in, line); myStream.str(""); myStream.clear() ; myStream << line ;
      if ( in.eof() ) break;
      if (version==1.0)
        myStream >> tms >> Dval >> Dval >>Dval >>Dval >>em->Temp >> Dval >> em->Vbias >> ImA >> em->x >> em->y >> em->z >> Ival ;
      if (version==1.1)
        myStream >> tms >> Dval >> Dval >>Dval >>Dval >>em->Temp >> Dval >> Dval >> ImA >> em->Vbias >> em->x >> em->y >> em->z >> Ival ;    
      if (version==1.2 || version==1.3)
        myStream >> tms >> Dval >> Dval >>Dval >>Dval >>em->Temp >> Dval >> em->Vbias >> ImA >> Dval >> em->x >> em->y >> em->z >> Ival >> Dval >> Dval ;
      
      em->Itot = ImA * 1.e-3;

      //Calculate the bin slice along each coordinate
      if (iloop==0) { x0=em->x ;y0=em->y ;z0=em->z ; }
      em->ix = (em->Ax!=0.)? 1 + TMath::Nint((em->x - x0)/(em->Ax/1000.)):1 ;
      em->iy = (em->Ay!=0.)? 1 + TMath::Nint((em->y - y0)/(em->Ay/1000.)):1 ;
      em->iz = (em->Az!=0.)? 1 + TMath::Nint((em->z - z0)/(em->Az/1000.)):1 ;
      
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
      
      getline(in, line); myStream.str(""); myStream.clear() ; myStream << line ;
      myStream >> em->At >> em->Nt ; em->At = em->At*1.e9 ; emh->At = em->At ;
      #if READONLYPDWVF==1
        getline(in, line); myStream.str(""); myStream.clear() ; myStream << line ;
      #endif
      for ( int i=0 ; i< em->Nt ; i++) {
	myStream >> em->volt[i] ;
	em->time[i]=i*em->At ; 
      }

      //----------------- Photodiode -------------------------------
      //getline(in, line);    
      
      Double_t LPower = 0. , LNph = 0. ;
      Double_t Vpd , Qpd = 0. , Bline = 0. , Vpdmax = -999999.9 , Eph;

      //All photodiode's voltages
      #if READONLYPDWVF==0  //If we only study the PD, we already read this line
        getline(in, line); myStream.str(""); myStream.clear() ; myStream << line ;
      #endif

      //Baseline of photodiode
      Int_t Nbl = TMath::Nint(5.0/em->At) ;
      for ( int i=0 ; i< Nbl ; i++) {
        #if READONLYPDWVF==0  //If we only study the PD, we have already read this line
	  myStream >> Vpd ;
	  Bline+= Vpd ;  
	#else
	  Bline+= em->volt[i] ;  
	#endif
      }
      Bline  = Bline/Nbl ;

      //Baseline corrected
      for ( int i=Nbl ; i< em->Nt ; i++) {

        #if READONLYPDWVF==0  //If we only study the PD, we have already read this line
	  myStream >> Vpd ;
	  Vpd = Vpd - Bline ;
	#else
	  Vpd = em->volt[i] - Bline ;
	#endif
        if ( Vpd > Vpdmax ) Vpdmax = Vpd ;
	if ( Vpd > 0. ) Qpd+=Vpd ; //This is really only an approximation to the total charge !!!

      }

      //Power and number of photons
      LPower = (emh->Lambda==1064.0)? Vpdmax/(ROSC*AWIR) : Vpdmax/(ROSC*AWRED) ;
      if ( TMath::Nint(emh->Lambda)==660  ) LPower = 0.450*(10.0*LPower) ; //Laser->[10% Photodiode, 45% eTCT] ]
      if ( TMath::Nint(emh->Lambda)==1064 ) LPower = 0.225*(10.0*LPower) ; //Laser->[10% Photodiode, 90%[22.5% TCTup, 22.5% TCTdown, 45% eTCT] ]
      Eph = 1240.0/emh->Lambda*JeV ; //in Jules
      //LNph   = Qpd*em->At*1.e-9/(ROSC*RESPONS*Eph) ;
      LNph   = Qpd*em->At*1.e-9/(ROSC*QELEC) ;  // Qtot= Neh * q_e

      /* 
	If we provide an average power (3rd argument in the command line), this
	means we want to correct all measurements to this power.
      */
      //if ( Pavg !=0.) 
	//for ( int i=0 ; i< em->Nt-1 ; i++) em->volt[i]= (LPower!=0.)? Pavg/LPower*em->volt[i+1] : 0. ;
      
      
      //----------------- End of photodiode ------------------------

      //Estimate polarity
      if ( TMath::Abs(TMath::MaxElement(em->Nt,em->volt)) > TMath::Abs(TMath::MinElement(em->Nt,em->volt)) ) Polarity++ ;
      else Polarity-- ;

      em->event=iRead ;
      
      //Now postprocess this entry (find out baseline, rtime and so on
      TWaveform *wv = new TWaveform( em ) ;
	
      if (iRead==0) tree->Branch("proc" , &wv , 32000 , 1 );
       
      if ( LPower!=0. ) {
        wv->LPower = LPower ;
        //wv->LPower = Pavg ;
	wv->LNph   = LNph ;
      }
      
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

/*----------------------------------------------------------*/
/**
 *
 * @param filename
 * @param em
 * @param tree
 * @param Cend
 * @return
 */
int parse_TPA( const char *filename, TMeas *em , TTree *tree , Double_t Cend ) {
    
    //Note that there is a difference in the coordinates. He calls X to my Y.
    //Open file
    ifstream in( filename );    

    string line, what ;
    for (Int_t iloop = 1 ; iloop<= 4 ; iloop ++ ) getline(in, line);
    stringstream myStream(line);
    //Double_t version ; myStream >> what >> version ;
    
    Double_t version = 1.0;
    
    //StartTime
    for (Int_t iloop = 5 ; iloop<= 5 ; iloop ++ ) getline(in, line);
    TString sdate = TString( line ) ;
    UShort_t  dd, mm, yy , hh , mn, ss;
    dd=atoi( sdate(11,2).Data() );
    mm=atoi( sdate(14,2).Data() );
    yy=atoi( sdate(17,4).Data() );
    hh=atoi( sdate(22,2).Data() );
    mn=atoi( sdate(25,2).Data() );
    //ss=atoi( sdate(28,2).Data() );
    ss=0;
    TDatime date ;
    date.Set(yy, mm, dd, hh, mn, ss) ;
    em->utc = date ;

    //comment
    TString comment = TString("This is a TPA measurement") ;
    
    //Wavelength
    Double_t wvlength=1300.0 ;
    
    //Amp gain
    Double_t Gain = 100. ;

    //Fluence and annealing
    Double_t Phi =0., iann =0. ;
 
    //Frequency
    //Int_t Freq = 200 ;
    
    //Absolute (X,Y,Z)
    for (Int_t iloop = 6 ; iloop<= 9 ; iloop ++ ) getline(in, line);
    getline(in, line); myStream.str(""); myStream.clear() ; myStream << line ;
    Double_t Xabs ; myStream >> what >> what >> Xabs ;

    getline(in, line); myStream.str(""); myStream.clear() ; myStream << line ;
    Double_t Yabs ; myStream >> what >> what >> Yabs ;

    getline(in, line); myStream.str(""); myStream.clear() ; myStream << line ;
    Double_t Zabs ; myStream >> what >> what >> Zabs ;
        
    //Ax,Nsteps
    getline(in, line); myStream.str(""); myStream.clear() ; myStream << line ;
    myStream >> what >> what >> em->Ax >> em->Nx ;
    em->Ax  = em->Ax * 1000. ; 
 
    getline(in, line); myStream.str(""); myStream.clear() ; myStream << line ;
    myStream >> what >> what >> em->Ay >> em->Ny ;
    em->Ay  = em->Ay * 1000. ; 
 
    getline(in, line); myStream.str(""); myStream.clear() ; myStream << line ;
    myStream >> what >> what >> em->Az >> em->Nz ;
    em->Az  = em->Az * 1000. ; 
 
    //Total number of Scans
    Int_t eNx = (em->Nx==0)?1:em->Nx , eNy = (em->Ny==0)?1:em->Ny , eNz = (em->Nz==0)?1:em->Nz ;
    Int_t NumOfScans = eNx*eNy*eNz ;
        
    //Top/Bottom/edge
    //Try to decide on the type of scan this is 0=tct std, 1=eTCT (default), 2=TCT+ , 3=TPA
    Int_t Illum ;
    Illum = 3;   //Discarded variable, not used
    
    // NVoltages
    em->NV = 1 ; 

    //Scope time window [s] and Nt
    getline(in, line); 
    myStream.str(""); myStream.clear() ; myStream << line ;
    Double_t Tbase /*, At*/ ; myStream >> what >> what >> Tbase >> em->Nt ;
    em->At = Tbase / em->Nt ;
    
    //Change the time scale to [ns]
    em->At = em->At * 1.e9;
    
     //Number of averages
    getline(in, line); 
    myStream.str(""); myStream.clear() ; myStream << line ;
    Int_t samples, Nav ; myStream >> what >> what >> samples >> Nav ;
    
    //Create a TMeasHeader object to store info that does not depend on event
    TMeasHeader *emh=new TMeasHeader( em->NV ) ;
    emh->Lambda  = wvlength ;
    emh->NV      = em->NV;
    emh->comment = comment ;
    emh->Setup   = 6 ;
    emh->Fluence = Phi ;
    emh->Nav     = Nav ;
    emh->Gain    = Gain ;
    emh->iann    = iann ;
    emh->Illum   = Illum ;
    emh->Ax      = em->Ax;
    emh->Ay      = em->Ay;
    emh->Az      = em->Az;
    emh->Nx	 = em->Nx;
    emh->Ny	 = em->Ny;
    emh->Nz	 = em->Nz;
    emh->At      = em->At ;
   
    //Cend
    emh->Cend = Cend ;
    
    //Vbias vector
    //In the meantime ask the use
    cout << "What is the bias?" << endl;
    cin >> emh->vVbias[0] ;
    
    //emh->vVbias[0] = 80. ;
        
    cout << "This is a TPA measurement" << endl ;

    //Skip line contents
    for (Int_t iloop = 18 ; iloop<= 20 ; iloop ++ ) getline(in, line); 
        
    //READ WAVEFORMS
    int  iRead=0 , iactual = 1 ;
    int  Polarity = 0 ;
    Double_t /*ImA ,*/ LPower;
   // Double_t valx,valy,valz;
    
    Double_t x0 , y0 , z0 ;
    Double_t *V1=new Double_t [em->Nt] ; //CORRECTION OF SPA CONTAMINATION ON TPA 
   // Double_t I1 ;
    for (Int_t iloop = 0 ; iloop < NumOfScans ; iloop++ ) {
      
      getline(in, line); myStream.str(""); myStream.clear() ; myStream << line ;
      if ( in.eof() ) break;
      if (version==1.0) myStream >> em->x >> em->y >> em->z >> LPower ;  
      //if (version==1.0) myStream >> em->z >> em->y >> em->x >> LPower ;
      //if (version==1.0) myStream >> em->z >> em->y >> em->x >> LPower ;LPower = LPower * 0.885 ;
      //if (version==1.0) myStream >> em->z >> em->y >> em->x >> LPower ;LPower = LPower * 9.528 ;
      //if (version==1.0) myStream >> em->z >> em->y >> em->x >> LPower ;LPower = LPower * 1.031 ;
      em->x = Xabs + em->x ; em->y = Yabs + em->y ; em->z = Zabs + em->z ; 
      
      //Calculate the bin slice along each coordinate
      if (iloop==0) { x0=em->x ;y0=em->y ;z0=em->z ; }
      em->ix = (em->Ax>0.)? 1 + TMath::Nint((em->x - x0)/(em->Ax/1000.)) : 1 ;
      em->iy = (em->Ay>0.)? 1 + TMath::Nint((em->y - y0)/(em->Ay/1000.)) : 1 ;
      em->iz = (em->Az>0.)? 1 + TMath::Nint((em->z - z0)/(em->Az/1000.)) : 1 ;
      
      em->Temp = 20.0 ;
      em->Itot = 0 * 1.e-3;
      
      //em->Vbias = 80. ;
      em->Vbias = emh->vVbias[0] ;
      
      //TO BE UPDATED: TIME OF EACH SHOT should go here

      
      for ( int i=0 ; i< em->Nt ; i++) {
	myStream >> em->volt[i] ;
	em->time[i]=i*em->At ; 
      }
      
      #if SPACORR==1
        if (iRead==SPAEVENT)  {
	  for ( int i=0 ; i< em->Nt ; i++) V1[i]=em->volt[i] ; 
          I1 = LPower ; 
	  cout<<"Applying SPA correction"<<endl;
	}
	//Do not correct first waveform at all. Maybe good to avoid jumps in baseline, and other plots
        if (iRead!=SPAEVENT) for ( int i=0 ; i< em->Nt ; i++) em->volt[i] = ( em->volt[i] - LPower/I1 * V1[i] )*I1*I1 ; 
      #endif
      
      //TO BE UPDATED: CALCULATE NUMBER OF PHOTONS

      //Estimate polarity
      if ( TMath::Abs(TMath::MaxElement(em->Nt,em->volt)) > TMath::Abs(TMath::MinElement(em->Nt,em->volt)) ) Polarity++ ;
      else Polarity-- ;

      em->event=iRead ;
      
      //Now postprocess this entry (find out baseline, rtime and so on
      TWaveform *wv = new TWaveform( em ) ;
	
      if (iRead==0) tree->Branch("proc" , &wv , 32000 , 1 );
       
      if ( LPower!=0. ) {
        wv->LPower = LPower ;
        //wv->LPower = Pavg ;
	wv->LNph   = 1.e6 ;
      }
      
      if ( iRead==DISCARDEVENT ) {      //SKIP EVENT WITH A GLITCH
        delete wv ; 
        continue ;
      }
      
      #if SPACORR==0
        tree->Fill() ;
      #else
	if (iRead!=SPAEVENT) tree->Fill() ;
      #endif

      
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
    delete [] V1 ;
        
    //delete emh ;

    return(1) ;

}
/*----------------------------------------------------------*/
/**
 *
 * @param filename
 * @param em
 * @param tree
 * @param Cend
 * @return
 */
int parse_tctUHH( const char *filename, TMeas *em , TTree *tree , Double_t Cend ) {
    
    //Open file
    ifstream in( filename );    
    
    //Skip few lines
    string line ;
    for (Int_t iloop = 1 ; iloop<= 4 ; iloop ++ ) getline(in, line);

    //Get annealing time[min]@T[C]
    string what ; char c;
    getline(in, line); 
    stringstream myStream(line);
    Double_t tann, Tann ; myStream >> tann >> c >>  Tann ;
    
    //Skip
    for (Int_t iloop = 1 ; iloop<= 5 ; iloop ++ ) getline(in, line);

    //Frequency
    getline(in, line); myStream.str(""); myStream.clear() ; myStream << line ;
    Double_t freq ; myStream >> freq >> what ;

    //Skip
    for (Int_t iloop = 1 ; iloop<= 3 ; iloop ++ ) getline(in, line);
    
    //Lambda
    getline(in, line); myStream.str(""); myStream.clear() ; myStream << line ;
    Double_t wvlength ; myStream >> wvlength >> what ;
    
    //Illumination type
    getline(in, line); getline(in, line); myStream.str(""); myStream.clear() ; myStream << line ;
    myStream >> what  ;
    Int_t illum = -1 ;
    if (what.find("front") == 0 ) illum=1;
    
    //comment
    for (Int_t iloop = 1 ; iloop<= 3 ; iloop ++ ) getline(in, line);
    getline(in, line); 
    TString comment = TString( line ) ;
        
    //Number of steps

    //Nx
    getline(in, line); //Skip
    getline(in, line); myStream.str(""); myStream.clear() ; myStream << line ;
    Double_t Dval ;
    myStream >> Dval ; 
    em->Nx  = TMath::Nint(Dval) ;

    //Ax
    getline(in, line); myStream.str(""); myStream.clear() ; myStream << line ;
    myStream >> em->Ax ;
   
    //Ny
    getline(in, line); myStream.str(""); myStream.clear() ; myStream << line ;
    myStream >> Dval ; 
    em->Ny  = TMath::Nint(Dval);

    //Ay
    getline(in, line); myStream.str(""); myStream.clear() ; myStream << line ;
    myStream >> em->Ay ;
   
    //Nz
    getline(in, line); myStream.str(""); myStream.clear() ; myStream << line ;
    myStream >> Dval ; 
    em->Nz  = TMath::Nint(Dval);

    //Az
    getline(in, line); myStream.str(""); myStream.clear() ; myStream << line ;
    myStream >> em->Az ;
   
    //Electrical scan parameters
    getline(in, line);
    Double_t V0 , Vf , dV , NV ;
    getline(in, line); myStream.str(""); myStream.clear() ; myStream << line ;
    myStream >> V0  ;
    getline(in, line); myStream.str(""); myStream.clear() ; myStream << line ;
    myStream >> Vf  ;
    getline(in, line); myStream.str(""); myStream.clear() ; myStream << line ;
    myStream >> NV  ;
    getline(in, line); myStream.str(""); myStream.clear() ; myStream << line ;
    myStream >> dV  ;
    em->NV = TMath::Nint( NV+1 ) ;
   
    //Skip
    for (Int_t iloop = 1 ; iloop<= 6 ; iloop ++ ) getline(in, line);
    
    //Acquisition Time
    Double_t AcqTime;
    getline(in, line); myStream.str(""); myStream.clear() ; myStream << line ;
    myStream >> AcqTime ;  //ns


    //Nav
    for (Int_t iloop = 1 ; iloop<= 3 ; iloop ++ ) getline(in, line);
    getline(in, line); myStream.str(""); myStream.clear() ; myStream << line ;
    Int_t Nav ; myStream >> Nav ;  
    
    //Skip
    for (Int_t iloop = 1 ; iloop<= 23 ; iloop ++ ) getline(in, line);
    
    
    //Total number of Scans
    Int_t NumOfScans = em->NV * (em->Nx+1)* (em->Ny+1)* (em->Nz+1) ;
    
    
    //Create a TMeasHeader object to store info that does not depend on event
    TMeasHeader *emh=new TMeasHeader( em->NV ) ;
    emh->Lambda = wvlength ;
    emh->NV = em->NV;
    emh->comment=comment ;
    emh->tann = tann ;
    emh->TempAnn = Tann ;
    emh->Freq = freq ; 
    emh->Illum = illum ;
    emh->Nx = em->Nx ;
    emh->Ax = em->Ax;
    emh->Ny = em->Ny ;
    emh->Ax = em->Ay;
    emh->Nz = em->Nz ;
    emh->Ax = em->Az;
    emh->Nav = Nav ; 
    
    //Vbias vector
    for (Int_t i=0 ; i< em->NV ; i++ )  emh->vVbias[i] = V0 + i*dV ;

    //Cend
    emh->Cend = Cend ;
        
    
    //Try to decide on the type of scan this is 0=tct std, 1=eTCT (default), 2=TCT+
    emh->Setup = 3 ;
    cout << "Hbg measurement"   << endl ;
     
    //If bline at 0V needs to be removed, then skip all lines until Vbias==0, read it and store it before hand
    Double_t *V0bl = new Double_t[em->Nt] ;
    for ( int i=0 ; i< em->Nt ; i++)  V0bl[i] = 0. ; 
    #if BL0SUBS==1
       long WaveformsStart = in.tellg();
       for (Int_t iloop = 0 ; iloop < NumOfScans ; iloop++ ) {
	 getline(in, line); myStream.str(""); myStream.clear() ; myStream << line ;
	 if ( in.eof() ) break;
	 Double_t Vval ; myStream >> Dval >> Dval >> Dval >> Dval >> Vval >> Dval >> Dval >> Dval >> Dval >> Dval >> Dval >> Dval >> Dval >> Dval ;
	 if ( Vval != 0 ) continue ;   
	 //We have Vbias=0
	 for ( int i=0 ; i< em->Nt ; i++)  myStream >> V0bl[i] ; 
       }
       in.seekg( WaveformsStart );
    #endif
        
    int  iRead=0 , iactual = 1 ;
    int  Polarity = 0 ;
    Double_t x0 , y0 , z0 ;
    for (Int_t iloop = 0 ; iloop < NumOfScans ; iloop++ ) {
      
      getline(in, line); myStream.str(""); myStream.clear() ; myStream << line ;
      if ( in.eof() ) break;
      myStream >> Dval >> em->x >> em->y >> em->z >> em->Vbias >> Dval >> Dval >> em->Itot >> Dval >> em->Temp >> Dval >> Dval >> Dval >> em->At ;
      
      //Calculate the bin slice along each coordinate
      if (iloop==0) { x0=em->x ;y0=em->y ;z0=em->z ; }
      em->ix = (em->Ax!=0.)? 1 + TMath::Nint((em->x - x0)/em->Ax):1 ;
      em->iy = (em->Ay!=0.)? 1 + TMath::Nint((em->y - y0)/em->Ay):1 ;
      em->iz = (em->Az!=0.)? 1 + TMath::Nint((em->z - z0)/em->Az):1 ;
            
      em->At = em->At*1.e9 ;
      em->Nt = TMath::Nint(AcqTime/em->At) ;
      
      Double_t val ;
      for ( int i=0 ; i< em->Nt ; i++) {
	myStream >> val ; em->volt[i] = val - V0bl[i] ;
	em->time[i]=i*em->At ; 
      }
      
      //Estimate polarity
      if ( TMath::Abs(TMath::MaxElement(em->Nt,em->volt)) > TMath::Abs(TMath::MinElement(em->Nt,em->volt)) ) Polarity++ ;
      else Polarity-- ;

      em->event=iRead ;
      
      //Now postprocess this entry (find out baseline, rtime and so on
      TWaveform *wv = new TWaveform( em ) ;
	
      if (iRead==0) {
         emh->At = em->At ;
         tree->Branch("proc" , &wv , 32000 , 1 );
      }
       
      wv->LPower = 0. ;
      wv->LNph   = 0. ;
      
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
    delete V0bl ;
        
    //delete emh ;

    return(1) ;

}

/*----------------------------------------------------------*/
/**
 *
 * @param filename
 * @param em
 * @param tree
 * @param Cend
 * @return
 */
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
    Int_t Ival /*, t0s , tms*/ ;
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
    Double_t /*ImA,*/x0,y0,z0;
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

/*----------------------------------------------------------*/
/**
 *
 * @param filename
 * @param em
 * @param tree
 * @param Cend
 * @return
 */
int parse_IDLTS( const char *filename, TMeas *em , TTree *tree , Double_t Cend ) {
    
    //Note that there is a difference in the coordinates. He calls X to my Y.
    //Open file
    ifstream in( filename );    

    string line, what ;
    for (Int_t iloop = 1 ; iloop<= 3 ; iloop ++ ) getline(in, line);
    stringstream myStream(line);
    Double_t version ; myStream >> what >> version ;

    //Maybe it would be better to read by keywords instead of by line position   
    Int_t offset = 0  ;  
    if (version>=1.1) offset=1;
    
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

    //comment
    istart=iend+1; //7
    iend=istart+1; //8
    for (Int_t iloop = istart ; iloop<= iend ; iloop ++ ) getline(in, line);
    TString comment = TString( line ) ;
    
    //Wavelength
    istart=iend+1; //9
    iend=istart+1; //10
    for (Int_t iloop = istart ; iloop<= iend ; iloop ++ ) getline(in, line); 
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
    if ( direction.find("top")    ==0 ) Illum=1;
    if ( direction.find("edge")   ==0 ) Illum=0;
    if ( direction.find("bottom") ==0 ) Illum=-1;
    
    //Amp gain
    Double_t Gain ;
    getline(in, line);  
    myStream.str(""); myStream.clear() ; myStream << line ;
    myStream >> what >> Gain ;


    //Fluence and annealing
    Double_t Phi =0., iann =0. ;
    if (version ==1.3 ) {
      getline(in, line); //read and discard biasteeresistance
      
      //Fluence
      getline(in, line); 
      myStream.str(""); myStream.clear() ; myStream << line ;
      myStream >> what >> Phi ;
      
      //Annealing
      getline(in, line); 
      myStream.str(""); myStream.clear() ; myStream << line ;
      myStream >> what >> iann ;
      
    }
    
    //Frequency
    getline(in, line); 
    myStream.str(""); myStream.clear() ; myStream << line ;
    Int_t Freq ; myStream >> what >> Freq ;
    
    //Number of averages
    getline(in, line); 
    myStream.str(""); myStream.clear() ; myStream << line ;
    Int_t Nav ; myStream >> what >> Nav ;
    
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
    emh->NV     = em->NV;
    emh->comment=comment ;
    emh->Setup  = 6 ;
    emh->Fluence= Phi ;
    emh->Nav    = Nav ;
    emh->Gain   = Gain ;
    emh->iann   = iann ;
    emh->Illum  = Illum ;
   
    //Cend
    emh->Cend = Cend ;
    
    //Vbias vector
    getline(in, line); myStream.str(""); myStream.clear() ; myStream << line ;
    myStream >> what ; for ( int i=0 ; i< em->NV ; i++) myStream >> emh->vVbias[i] ;
    
    //Nominal power
    Int_t nPower ; Double_t Power ;
    for (Int_t iloop = 1 ; iloop<= 4 ; iloop ++ ) getline(in, line); 
    myStream.str(""); myStream.clear() ; myStream << line ;
    myStream >> what ; myStream >> nPower ;
    getline(in, line); 
    myStream.str(""); myStream.clear() ; myStream << line ;
    myStream >> what ; myStream >> Power ;
    emh->Power = Power ;


    //Ax
    istart = 1; 
    iend   = istart+8; 
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
    iend   = istart+3; 
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
    
    if ( emh->Illum == 0 ) cout << "TCT+: this is an eTCT measurement" << endl ;
    if ( emh->Illum != 0 ) cout << "TCT+: this is a normal TCT measurement" << endl ;
    
    istart = iend+1; 
    iend   = istart+3; 
    for (Int_t iloop = istart ; iloop<= iend ; iloop ++ ) getline(in, line); 
    
    int  iRead=0 , iactual = 1 ;
    int  Polarity = 0 ;
    Double_t ImA, x0, y0, z0;
    for (Int_t iloop = 0 ; iloop < NumOfScans ; iloop++ ) {
      
      getline(in, line); myStream.str(""); myStream.clear() ; myStream << line ;
      if ( in.eof() ) break;

      if (version==1.0)
        myStream >> tms >> em->Temp >> Dval >>Dval >>Dval >>Dval >> Dval >> em->Vbias >> ImA >> em->pWidth >> em->pAmpl >> Dval ;
      if (version==1.4)
        myStream >> tms >> em->Temp >> Dval >>Dval >>Dval >>Dval >> Dval >> em->Vbias >> ImA >> Dval >>em->pWidth >> em->pAmpl >> Dval ;
      
      em->Itot = ImA * 1.e-3;

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
      

      //Open the file with the waveforms
      TString Dfnm ;
      getline(in, line); myStream.str(""); myStream.clear() ; myStream << line ;
      myStream >> Dfnm ;
      ifstream Din( Dfnm );    
            
      //Parse filename
      Ssiz_t Lund = Dfnm.Last('_');
     // Ssiz_t Ldot = Dfnm.Last('.');
      
      //Calculate the bin slice along each coordinate
      TString tx = Dfnm(0,Lund) ; Lund=tx.Last('_'); tx = tx(Lund+1,tx.Length()-Lund);
      em->x=tx.Atof() ; em->y = 0. ; em->z = 0. ;
      if (iloop==0) { x0=em->x ;y0=em->y ;z0=em->z ; }
      em->ix = (em->Ax!=0.)? 1 + TMath::Nint((em->x - x0)/em->Ax):1 ;
      em->iy = (em->Ay!=0.)? 1 + TMath::Nint((em->y - y0)/em->Ay):1 ;
      em->iz = (em->Az!=0.)? 1 + TMath::Nint((em->z - z0)/em->Az):1 ;
      

      //Get the (Nt,At)
      Char_t c ;
      for (Int_t it=1;it<=5;it++) getline(Din, line);
      myStream.str(""); myStream.clear() ; myStream << line ;
      myStream >> what >> c >> em->Nt >> c >> Dval ;
      
      for (Int_t it=6;it<=9;it++) getline(Din, line);
      myStream.str(""); myStream.clear() ; myStream << line ;
      myStream >> what >> c >> em->At >> c >> Dval ;

      em->At = em->At*1.e9 ; emh->At = em->At ;
      for (Int_t it=10;it<=22;it++) getline(Din, line);
      
      Double_t LPower = 0. , LNph = 0. ;
      Double_t Vpd[em->Nt] , Qpd = 0. , Bline = 0. , Vpdmax = -999999.9 , Eph ;
      for ( int i=0 ; i< em->Nt ; i++) {
        getline(Din, line); myStream.str(""); myStream.clear() ; myStream << line ;
	myStream >> Dval>>c>>em->volt[i]>>c>>Vpd[i] ;
	em->time[i]=i*em->At ; 
      }
      Din.ignore() ;
      Din.close() ;


      //Baseline of photodiode
      Int_t Nbl = TMath::Nint(5.0/em->At) ;
      for ( int i=0 ; i< Nbl ; i++) Bline+= Vpd[i] ;  
      Bline  = Bline/Nbl ;

      //Baseline corrected
      for ( int i=Nbl ; i< em->Nt ; i++) {
	Vpd[i] = Vpd[i] - Bline ;
        if ( Vpd[i] > Vpdmax ) Vpdmax = Vpd[i] ;
	if ( Vpd[i] > 0. ) Qpd+=Vpd[i] ; //This is really only an approximation to the total charge !!!

      }

      //Power and number of photons
      LPower = (emh->Lambda==1064.0)? Vpdmax/(ROSC*AWIR) : Vpdmax/(ROSC*AWRED) ;
      if ( TMath::Nint(emh->Lambda)==660  ) LPower = 0.450*(10.0*LPower) ; //Laser->[10% Photodiode, 45% eTCT] ]
      if ( TMath::Nint(emh->Lambda)==1064 ) LPower = 0.225*(10.0*LPower) ; //Laser->[10% Photodiode, 90%[22.5% TCTup, 22.5% TCTdown, 45% eTCT] ]
      Eph = 1240.0/emh->Lambda*JeV ; //in Jules
      //LNph   = Qpd*em->At*1.e-9/(ROSC*RESPONS*Eph) ;
      LNph   = Qpd*em->At*1.e-9/(ROSC*QELEC) ;  // Qtot= Neh * q_e

      /* 
	If we provide an average power (3rd argument in the command line), this
	means we want to correct all measurements to this power.
      */
      //if ( Pavg !=0.) 
	//for ( int i=0 ; i< em->Nt-1 ; i++) em->volt[i]= (LPower!=0.)? Pavg/LPower*em->volt[i+1] : 0. ;
      
      
      //----------------- End of photodiode ------------------------

      //Estimate polarity
      if ( TMath::Abs(TMath::MaxElement(em->Nt,em->volt)) > TMath::Abs(TMath::MinElement(em->Nt,em->volt)) ) Polarity++ ;
      else Polarity-- ;

      em->event=iRead ;
      
      //Now postprocess this entry (find out baseline, rtime and so on
      TWaveform *wv = new TWaveform( em ) ;
	
      if (iRead==0) tree->Branch("proc" , &wv , 32000 , 1 );
       
      if ( LPower!=0. ) {
        wv->LPower = LPower ;
        //wv->LPower = Pavg ;
	wv->LNph   = LNph ;
      }
      
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
/*----------------------------------------------------------*/
/**
 *
 * @param filename
 * @return
 */
int ReadNerOfBins( const char *filename ) {
   
     Int_t Ner = 0;
  
    //Open file
    ifstream in( filename );
    if (!in || in.bad()) {
      cout << "Error opening " << filename << endl ;
      exit(1);
    }
        
    //Read header
    string what , line , Type , Dfnm ;
    getline(in,line) ;  
    stringstream myStream(line);
    Int_t Nskip ;
    if        ( line.find("================") == 0 ) {
       getline(in,line) ;  
       if        ( line.find("SSD simulation") == 0 ) { 
         Nskip=60;
       } else {

	 getline(in,line) ;  
	 Double_t version ; 
	 myStream.str(""); myStream.clear() ; myStream << line ;
	 myStream >> what >> version ;
	 if (version == 1.0 ) Nskip=55 ;
	 if (version == 1.1 || version == 1.2 ) Nskip=56 ;
	 if (version == 1.3 ) Nskip=58;
	 
	 //Only IDLTS
	 getline(in,line) ;getline(in,line) ;
	 myStream.str(""); myStream.clear() ; myStream << line ; 
	 myStream >> what >> Type ;
	 if ( Type.find("DLTS")==0 ) {
	   //Get filename where number of points is stored
	   Nskip=55;
	   for (int i=5 ; i<=Nskip ; i++) getline(in,line) ;myStream.str(""); 
	   myStream.clear() ; myStream << line ;  
	   myStream>>Dfnm;
           ifstream din( Dfnm );
	   if (!din || din.bad()) {
	     cout << "Error opening " << Dfnm << endl ;
	     exit(1);
	   }
	   for (int i=0 ; i<5 ; i++) getline(din,line) ; 
	   myStream.str(""); myStream.clear() ; myStream << line ;
	   char c ;
	   myStream >> what >> c>>Ner>>c>>Ner;
	   in.ignore() ;
	   in.close() ;
	   din.ignore() ;
	   din.close() ;
           return(Ner) ;
	 }
	 
         for (int i=4 ; i<Nskip ; i++) getline(in,line) ; 
       }
       
    } else if ( line.find("----------------")      == 0 ) {        //TPA
       Nskip=17;  //TPA
       for (int i=2 ; i<Nskip ; i++) getline(in,line) ; 

    } else if ( line.find("MTCT Header")      == 0 ) {             //UHH
       Nskip=69;  //UHH
       for (int i=2 ; i<Nskip ; i++) getline(in,line) ; 

    } else if ( line.find("* t0      :")      == 0 ) {             //IFCA
       getline(in,line) ;  getline(in,line) ;  
       myStream.str(""); myStream.clear() ; myStream << line ;
       char c ;
       myStream >> c >> what >> c >> Ner >> c ;
       cout << "Waveforms have " << Ner << " bins." <<endl ;

       in.ignore() ;
       in.close() ;

       return(Ner) ;
      
    } else {
      Nskip= 15;                                              //Standard eTCT
      for (int i=2 ; i<Nskip ; i++) getline(in,line) ; 
    }
    
       
    myStream.str(""); myStream.clear() ; myStream << line ;
    if ( Nskip==55 || Nskip==56 || Nskip==58 ) {
      Double_t val ;
      myStream  >> val >> Ner ;
    } else if (Nskip==15 ){
      // NTimePoints:    10003
      myStream  >> what >> Ner ;	     //14 NTimePoints:  10003
    } else if (Nskip==17 ){
      // timebase[s], datapoints	1.000000E-7	1000
      myStream  >> what >> what >> what >> Ner ;
    } else if (Nskip==51 ){
      // NTimePoints:    10003
      myStream  >> what >> Ner ;	     //14 NTimePoints:  10003
    } else if (Nskip==53 ){
      myStream  >> Ner ;	     //14 NTimePoints:  10003
    } else if (Nskip==60 ){
      for (int i=2 ; i<Nskip ; i++) getline(in,line) ; 
      myStream.str(""); myStream.clear() ; myStream << line ;
      myStream  >> Ner ;	     //TRACS
    } else {
      Ner = std::count( line.begin(), line.end(), '\t' ) - 14 + 1 ;
    }
    cout << "Waveforms have " << Ner << " bins." <<endl ;

    in.ignore() ;
    in.close() ;
    
    return(Ner) ;
}

/*----------------------------------------------------------*/
int main(int argc , char *argv[]) {

  //for (int i=1 ; i<argc ; i++)  tree (1,argv[i]);
  if ( argc ==2 ) tree (1,argv[1]);
  if ( argc ==3 ) tree (2,argv[1],argv[2]);   //argv[2]=Cend [pF]
  
}

