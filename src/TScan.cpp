/* 
Changelog:

8Feb2016: 

1) For TPA HVCMOS. I hardened the WeightedMean RMS filtering in TScan::RoughTRTL. The 2
consecutive RMS filters run with 0.8 and 0.6 times the RMS.
2) TScan::VsBias(), introducing case where NV=1

*/

#include "TScan.h"

extern Double_t        *Getz0t0( char *what , char *fnm , char *zlimit, TString ftype=TString("Erf") ) ;
extern vector<TString> plot1D( char *varexp , char *selection, int argc , ... ) ;
extern int             plot2D( char *varexp , char *selection, char *fnm ) ;
extern Int_t           plotvd( TString coord, char *zslection , Double_t t0, Double_t z0 , Double_t BLMean, char *selection, char *fnm , Double_t toffset=ATRISE , int i2corr=0 ) ;
extern Int_t           plotvd( TString coord, char *zslection , Double_t t0, Double_t z0 , Double_t Cend , Double_t At, Double_t BLMean, char *selection, char *fnm , Double_t toffset=ATRISE , int i2corr=0) ;
extern Double_t        *FitErf ( TH1D *h1 ) ;
extern Double_t        *LeftFitErf ( TH1D *h1 ) ;
Double_t               GaussBox(Double_t *x , Double_t *pars ) ;
Double_t               GaussAssyBox(Double_t *x , Double_t *pars ) ;
Double_t               GaussExpBox( Double_t *x, Double_t *pars ) ;
Double_t               GaussExpBoxBline( Double_t *x, Double_t *pars ) ;

ClassImp(TScan)

/* Default constructor, see ROOT manual page 270 */
TScan::TScan(  ) {

  Polarity = 0;
  Vbmax = Vbmin = Vmax = Vmin = t0 = z0 = 0.;
  BLMean= thicknessx= thicknessy =thicknessz = dz = dt = 0.;
  Fluence = iann = tann	= Etann	= TempAnn = Ldiffx = Ldiffy = Ldiffz = 0.;
  FwQx=FwQy=FwQz=FWHMQx = FWHMQy=FWHMQz=RMSQx=RMSQy=RMSQz=0;
  cce = Itot = 0. ;

}

TScan::TScan( char *fnm , Double_t Cend ) {
          
     
     /* OPEN RAW DATA FILES */
     
     this->fnm = fnm ;
     this->Cend = Cend ;

     char *pdot , *pbnm  ;
     pbnm  = basename( fnm ) ;
     
     //Get filename without path
     pdot  = strstr( pbnm , ".") ;
     Int_t nchar = pdot-pbnm ;
     bnm = new char [nchar+1];
     strncpy( bnm , pbnm , nchar) ;
     *(bnm+nchar) = '\0' ;
     
     //Find sensor name: either yyyymmddhhmmss_name.ext or yyyy-mm-dd_hh-mm-ss_name.ext
     string name(bnm);
     Int_t iund = name.rfind('_');
     name=name.substr(iund+1,name.length()-iund-1);
     
     //Get path
     nchar = pbnm-fnm ;
     dnm = new char [nchar+1];
     strncpy( dnm , fnm , nchar) ;
     *(dnm+nchar) = '\0' ;

     TString tfnm = TString( fnm )   ;

     TFile *f = new TFile( tfnm.Data() );
     if ( !f->IsOpen() ) exit(-1) ;

     TTree *tree = (TTree*) f->Get("edge");  

     TMeas   *em   = 0 ;
     TBranch *raw  = tree->GetBranch("raw") ;
     raw->SetAddress(&em) ;
     raw->GetEntry(0) ;  
     x0file=em->x;   
     y0file=em->y;   
     z0file=em->z;   
     
       
     TWaveform *wv = new TWaveform( em ) ;     
     TBranch *proc = tree->GetBranch("proc") ;
     proc->SetAddress(&wv) ;
     
     TMeasHeader *emh = (TMeasHeader *) tree->GetUserInfo()->At(0) ;
    
     tree->GetEntry(0) ;  //Commented 19Sept2012
          

     /* CREATE VECTORS  */
     
     dx = emh->Ax ; dy = emh->Ay ; dz = emh->Az ;
     Nx = emh->Nx ; Ny = emh->Ny ; Nz = emh->Nz ;
     ScanDir[0] =  (TMath::Abs(dx)>0.)? 1 : 0 ;
     ScanDir[1] =  (TMath::Abs(dy)>0.)? 1 : 0 ;
     ScanDir[2] =  (TMath::Abs(dz)>0.)? 1 : 0 ;  //edgeTCT
     ScanDir[3] =  (NV>0.)? 1 : 0 ;  
     
     dt     = emh->At ;
     NV     = emh->NV ;
     Vbhalf = emh->vVbias[NV/2] ;	
     vhsV_zi.reserve(NV)   ;  vlhsV_zi.reserve(NV) ;
     vhsV_Vbi.reserve(NV)  ;  vlhsV_Vbi.reserve(NV) ;
     vhvtz.reserve(NV)     ;  vhQ2D.reserve(NV) ; 
     vhvd2D.reserve(NV)    ;  vhvdEw2D.reserve(NV)    ;
     vhvtx.reserve(NV)     ;  vhtcoll2D.reserve(NV)  ;
     vhvty.reserve(NV)     ;
     CCE.reserve(NV)       ;  
     
     //Vectors that will go in the AutoPilot tree
     if ( Nx > 0 ) {
       Qx     = new Double_t [Nx];
       tcollx = new Double_t [Nx];
       x      = new Double_t [Nx];
       rQx    = new Double_t [Nx];
     }     
     if ( Ny > 0 ) {
       Qy     = new Double_t [Ny];
       tcolly = new Double_t [Ny];
       y      = new Double_t [Ny];
       rQy    = new Double_t [Ny];
     }     
     if ( Nz > 0 ) {
       Qz     = new Double_t [Nz];
       vdrz   = new Double_t [Nz];
       tcollz = new Double_t [Nz];
       z      = new Double_t [Nz];
       rQz    = new Double_t [Nz];
     }     
          
     Polarity = emh->GetFilePolarity() ;
     #if AUTO==0
     Polarity=POLARITY ;
     #endif
          
     //Vbmax  = ( Polarity == 1 ) ? emh->vVbias.Max() : emh->vVbias.Min() ;
     //Vbmin  = ( Polarity == 1 ) ? emh->vVbias.Min() : emh->vVbias.Max() ;

     //Robust method that works even with inverting amplifiers
     TVectorD AvV = emh->vVbias ;  AvV.Abs();
     Vbmax  = AvV.Max() ;
     Vbmin  = AvV.Min() ;
     
     Int_t sign = TMath::Nint(Vbhalf/AvV[NV/2]) ;
     Vbmax  = sign * Vbmax ;
     Vbmin  = sign * Vbmin ;
          
     /* GET VBIAS VECTOR */
     
     Vb         = new Double_t [NV] ;
     for (Int_t i=0 ; i<NV ; i++) Vb[i]=emh->vVbias[i] ;

     /* Vectors to dump averages */
     vTemp = new Double_t [NV] ;
     vItot = new Double_t [NV] ;
     
     /* Depleted thickness vs voltage */
     vFwQx    = new Double_t [NV] ; //Thickness calculated from fit to assumed function 
     vFwQy    = new Double_t [NV] ; //Thickness calculated from fit to assumed function 
     vFwQz    = new Double_t [NV] ; //Thickness calculated from fit to assumed function 
     
     vFwQx0   = new Double_t [NV] ; //Starting value for the fit		        
     vFwQy0   = new Double_t [NV] ; //Starting value for the fit		        
     vFwQz0   = new Double_t [NV] ; //Starting value for the fit		        
     
     vFWHMQx  = new Double_t [NV] ; //Qx FWHM					        
     vFWHMQy  = new Double_t [NV] ; //Qy FWHM					        
     vFWHMQz  = new Double_t [NV] ; //Qz FWHM					        
     
     vRMSQx   = new Double_t [NV] ; //Qx width measured as RMS					   
     vRMSQy   = new Double_t [NV] ; //Qy width measured as RMS					   
     vRMSQz   = new Double_t [NV] ; //Qz width measured as RMS
     
     vcQx     = new Double_t [NV] ; //Center of depletion width 		        
     vcQy     = new Double_t [NV] ; //Center of depletion width 		        
     vcQz     = new Double_t [NV] ; //Center of depletion width 		        
     vcQx0    = new Double_t [NV] ; //Starting value of the fit for vcQx	        
     vcQy0    = new Double_t [NV] ; //Starting value of the fit for vcQy	        
     vcQz0    = new Double_t [NV] ; //Starting value of the fit for vcQz	        
     
     vLdiffx  = new Double_t [NV] ; vLdiffy  = new Double_t [NV] ; vLdiffz  = new Double_t [NV] ;
     vGwidthx = new Double_t [NV] ; vGwidthy = new Double_t [NV] ; vGwidthz = new Double_t [NV] ;
     vFitResx = new Int_t [NV]    ; vFitResy = new Int_t [NV]	 ; vFitResz = new Int_t [NV]	;

     // Vector for number of e-h normalization
     N0.ResizeTo(NV) ;
          
     /* GET BlineMean AVERAGE FOR EACH Vbias */
     
     Double_t *vBlineMean = new Double_t [NV] ;
     TCanvas *cdum = new TCanvas("cdum","cdum",600,400);
     TH1D *h1;
     for (Int_t i=0 ; i<NV ; i++) {
       TString selection = Form( "Vbias==%f", Vb[i] ) ;
       //TString selection = Form( "Vbias==%f &&(z<4.94 || z>5.)", Vb[i] ) ;  //Sick cases where Bline is too short, this may help. 
       								            //In thoses cases, BlineMean changes with z
       tree->Draw("BlineMean",selection.Data()) ;
       h1 = (TH1D*) gPad->GetPrimitive("htemp"); 
       vBlineMean[i] = h1->GetMean(); 
     }
     BLMean = TMath::Mean(NV,vBlineMean);
     delete cdum ;
     delete [] vBlineMean ;

     
     /* 
     ESTIMATE <tright> FROM THE MIDDLE VOLTAGE 
     (it will be smaller than the max, so wider pulses) 
     */
     Double_t trtl[4] ;
     //TString vselection = Form("Vbias==%d",(Int_t) Vbhalf ) ;
     #ifdef ZLIMIT
       TString vselection = Form("Vbias==%d && %s",(Int_t) Vbhalf ,ZLIMIT) ;
     #else
       TString vselection = Form("Vbias==%d",(Int_t) Vbhalf ) ;
     #endif
     RoughTRTL( vselection , tree , wv , trtl );
     travg = (ScanDir[2])? trtl[1] : trtl[1] + 10.  ;
     

     /* ESTIMATE T0, <tleft> FROM THE HIGHEST VOLTAGE */
     #ifdef ZLIMIT
       vselection = Form("Vbias==%d && %s",(Int_t) Vbmax ,ZLIMIT) ;
     #else
       vselection = Form("Vbias==%d",(Int_t) Vbmax ) ;
     #endif
     
     #if DECONV==0
       RoughTRTL( vselection , tree , wv , trtl );
     #else
       RoughTRTL( vselection , tree , trtl , Cend );
     #endif
     tlavg = trtl[0] ;
     Atrl  = ( ATRL==0 ) ? travg-tlavg : ATRL ;
     Vmax  = trtl[2] ; 
     Vmin  = trtl[3] ; 
     
     /* Estimate the FWHM of the Q50 without any model assumption */

     Char_t coord ;
     for ( Int_t c=0 ; c<=2 ; c++) {

       if (c==0) coord='x' ; if (c==1) coord='y' ; if (c==2) coord='z' ;
       qfwhm[c] = 0. ;
       if ( ScanDir[c] ) {
         tree->Draw(Form("Q50:%c",coord),Form("Vbias==%f",Vbmax)) ;
	 TGraph *graph = (TGraph*)gPad->GetPrimitive("Graph"); // 2D
	 Int_t nentries = graph->GetN();
	 Double_t *zarr = tree->GetV2(); //do not delete these arrays, done by TTree
	 Double_t *qarr = tree->GetV1();
	 TH1D *hstat = new TH1D("hstat","hstat",nentries,zarr[0],zarr[nentries-1]);
	 for (Int_t il=0;il<nentries;il++) hstat->SetBinContent(il+1,qarr[il]);
	 qfwhm[c] = 2.355*hstat->GetRMS() ;
	 delete hstat ;  
       }

     }
      
      /* (T0,Z0) Refined calculations */
      t0 = tlavg ;
      Double_t tsmallest = tlavg ;
      if ( ScanDir[0] ) CalcWhereSignal( 0 );  if ( t0 < tsmallest ) tsmallest= t0;
      if ( ScanDir[1] ) CalcWhereSignal( 1 );  if ( t0 < tsmallest ) tsmallest= t0;
      if ( ScanDir[2] ) CalcWhereSignal( 2 );  if ( t0 < tsmallest ) tsmallest= t0;
      t0=tsmallest ;
      
      Print_t0( tree , t0 , Vbmax , Cend, dt, dnm , bnm );     
     
     /* Copy annealing info (if any) */
     
     iann    = emh->iann    ;   
     tann    = emh->tann    ;	
     Etann   = emh->Etann   ;  
     TempAnn = emh->TempAnn ;
     Fluence = emh->Fluence ;
     if (Fluence == 0.) Fluence=0.1; //Useful for logX plots
     
     /* Get some magnitudes vs bias voltage */
     VsBias( tree , em ) ;
     
     //Do a CCE plot
     PlotCCEt_vs_Vbias(tree);
     
     delete wv ;
     delete tree ;
     f->Close();   //It was commented up to April 13 2015

}

TScan::~TScan(  ) {
     
     delete [] Vb ;  N0.Clear()  ;
     vhsV_zi.clear() ;  vlhsV_zi.clear() ;
     vhsV_Vbi.clear() ;  vlhsV_Vbi.clear() ;
     vhvtz.clear() ;  vhQ2D.clear() ; 
     vhvd2D.clear()  ; vhvdEw2D.clear()  ; 
     vhvtx.clear() ;  vhtcoll2D.clear();
     vhvty.clear() ;
     CCE.clear()  ;
     
     delete [] vItot    ; delete [] vTemp    ; 
     delete [] Qx       ; delete [] Qy       ; delete [] Qz  ; 
     delete [] vFwQx    ; delete [] vFwQy    ; delete [] vFwQz ; 
     delete [] vFwQx0   ; delete [] vFwQy0   ; delete [] vFwQz0  ; 
     delete [] vFWHMQx  ; delete [] vFWHMQy  ; delete [] vFWHMQz ; 
     delete [] vRMSQx   ; delete [] vRMSQy   ; delete [] vRMSQz ; 
     delete [] vcQx     ; delete [] vcQy     ; delete [] vcQz  ; 
     delete [] vcQx0    ; delete [] vcQy0    ; delete [] vcQz0 ; 
     delete [] vLdiffx  ; delete [] vLdiffy  ; delete [] vLdiffz  ;
     delete [] vGwidthx ; delete [] vGwidthy ; delete [] vGwidthz ;
     delete [] vFitResx ; delete [] vFitResy ; delete [] vFitResz ;
}

//------------------------------------------------------------------------
void TScan::RoughTRTL( TString vselection , TTree *tree , TWaveform *wv , Double_t *trtl ) {

     tree->Draw( ">>myList" , (char*) vselection.Data() , "entrylistarray" ) ;
     TEntryListArray *list=(TEntryListArray*)gDirectory->Get("myList") ;
     tree->SetEntryList( list ) ; //Use tree->SetEventList(0) to switch off
                        	  
     Int_t  nev  = list->GetN() ; 
     if (nev==0) exit(-1) ;
     Double_t val , Vmax = -99999.9 , Vmin = 99999.9 ;
     
     vector<Double_t> tr , tl ;
     vector<Double_t>::iterator min , max ;
     for ( Int_t ii=0 ; ii < nev ; ii++ ) {

       Int_t iev = list->GetEntry(ii) ; 
       tree->GetEntry( iev );    
              
       tl.push_back( wv->GetTleft( ) ) ;
       tr.push_back( wv->GetTright( ) );
              
       val = wv->GetVmax() ;
       if (val > Vmax ) Vmax = val ; 
       val = wv->GetVmin() ;
       if (val < Vmin ) Vmin = val ; 

     }  
     
     Double_t Mean , RMS ;
     
     min = min_element( tl.begin() , tl.end() );
     max = max_element( tl.begin() , tl.end() );
     WeightedMean( tl , *min , *max , Mean , RMS )  ;
     
     /* TO CALCULATE TLEFT DO A WEIGHTED MEAN with RMS filtering */
     for ( Int_t i=0;i<5;i++ ) WeightedMean( tl , Mean , RMS ) ;
     Double_t tlavg = Mean ;

     min = std::min_element( tr.begin() , tr.end() );
     max = std::max_element( tr.begin() , tr.end() );
     WeightedMean( tr , *min , *max , Mean , RMS )  ;
     for ( Int_t i=0;i<5;i++ ) WeightedMean( tr , Mean , RMS ) ;
     Double_t travg = Mean ;
     if  (travg  <= tlavg )      travg=tlavg+8.0;
     if  (travg  - tlavg < 5.0 ) travg=tlavg+8.0;
     
     
     trtl[0]=tlavg;
     trtl[1]=travg;
     trtl[2]=Vmax;
     trtl[3]=Vmin;
     
     tr.clear();
     tl.clear();
     tree->SetEventList(0);
     
}

//------------------------------------------------------------------------
void TScan::RoughTRTL( TString vselection , TTree *tree , Double_t *trtl , Double_t Cend ) {

     //When we differenciate V(t;z) (so, we do dV(t;z)/dt), then t0 shifts
     

     TMeas   *em   = 0 ;
     TBranch *raw  = tree->GetBranch("raw") ;
     raw->SetAddress(&em) ;
     raw->GetEntry(0) ;  

     TWaveform *wv = new TWaveform( em ) ;     
     TBranch *proc  = tree->GetBranch("proc") ;
     proc->SetAddress(&wv) ;

     tree->Draw( ">>myList" , (char*) vselection.Data() , "entrylistarray" ) ;
     TEntryListArray *list=(TEntryListArray*)gDirectory->Get("myList") ;
     tree->SetEntryList( list ) ; //Use tree->SetEventList(0) to switch off
                       	  
     Int_t  nev  = list->GetN() ; 
     if (nev==0) exit(-1) ;
     Double_t val , Vmax = -99999.9 , Vmin = 99999.9 ;
     
     vector<Double_t> tr , tl ;
     vector<Double_t>::iterator min , max ;
     for ( Int_t ii=0 ; ii < nev ; ii++ ) {

       Int_t iev = list->GetEntry(ii) ; 
       tree->GetEntry( iev );    
       
       //We need to overwrite volt contents with the derivatieve: dvolt/dt
       Int_t Nt = wv->GetNbins() ;
       TMeas *dem = new TMeas( );
       dem->Nt = Nt ;
       dem->volt = new Double_t [ Nt ]; 
       dem->time = new Double_t [ Nt ]; 
       if (Cend!=0.) dem->Vbias = Vbmax ;
       for ( Int_t iv = 0 ; iv< Nt ; iv++ ) {        
	 dem->volt[iv] = ROSC*Cend*(em->volt[iv+1] - em->volt[iv])/dt + em->volt[iv] ;
	 dem->time[iv] = em->time[iv] ;
       }
       dem->event = iev ;
       TWaveform *dwv = new TWaveform( dem ) ;       
       
       tl.push_back( dwv->GetTleft( ) ) ;
       tr.push_back( dwv->GetTright( ) ) ;
       
       val = dwv->GetVmax() ;
       if (val > Vmax ) Vmax = val ; 
       val = dwv->GetVmin() ;
       if (val < Vmin ) Vmin = val ; 

       delete [] dem->volt ;
       delete [] dem->time ;
       delete dwv ;
     }  
     Double_t Mean  = TMath::Mean( tr.begin(), tr.end()) , RMS=TMath::RMS( tr.begin() , tr.end() ) ;
     Double_t travg = Mean+RMS ; 

     /* TO CALCULATE TLEFT DO A WEIGHTED MEAN with RMS filtering */
     Mean =  TMath::Mean( tl.begin() , tl.end() ) ; 
     RMS  =  TMath::RMS(  tl.begin() , tl.end() ) ;
     WeightedMean( tl , Mean , RMS ) ;
     WeightedMean( tl , Mean , RMS ) ;
     Double_t tlavg = Mean ;
          
     trtl[0]=tlavg;
     trtl[1]=travg;
     trtl[2]=Vmax;
     trtl[3]=Vmin;
     
     tr.clear();
     tl.clear();
     tree->SetEventList(0);
     delete wv ;
          
}

//------------------------------------------------------------------------

void TScan::CalcWhereSignal( Int_t coord ){

     /* USE ESTIMATED "T0" TO CALCULATE Efficiency AS Sum$(). THEN OBTAIN (TO,Z0) */

     Char_t xyz;
     if (coord==0) xyz = 'x' ; if (coord==1) xyz = 'y' ;  if (coord==2) xyz = 'z' ; 

     TString what ;
     #if DECONV==0
	 what = (ATRL==0)? Form("Sum$((volt-BlineMean)*(time>%f && time<%f))", tlavg, travg ) :
                           Form("Sum$((volt-BlineMean)*(time>%f && time<%f))", tlavg, tlavg + ATRL ) ;
	 if ( I2CORR==1 ) what=what+"/(LPower*LPower)";
	 what=what+Form(":%c", xyz);	  
     #else
	 what = (ATRL==0)? Form("Sum$(%5.2f*%7.5f*(volt[Iteration$+1]-volt[Iteration$])/%f+(volt-BlineMean)*(time>%f && time<%f)):%c", ROSC,Cend,dt,tlavg, travg , xyz) :
                           Form("Sum$(%5.2f*%7.5f*(volt[Iteration$+1]-volt[Iteration$])/%f+(volt-BlineMean)*(time>%f && time<%f)):%c", ROSC,Cend,dt,tlavg, tlavg+ATRL , xyz) ;	    
     #endif

     Double_t *z0t0;
     #ifdef ZLIMIT
       if (qfwhm[coord]>0.03) z0t0  = Getz0t0( (char *) what.Data() , fnm , ZLIMIT ) ;
       else                  z0t0  = Getz0t0( (char *) what.Data() , fnm , ZLIMIT , TString("gaus")) ; 
     #else
       TString empty="";  
       if ( coord ==0 && ScanDir[1] ) empty = empty+Form("TMath::Abs(y-%lf)<%lf",y0file+TMath::Nint(Ny/2.)*dy,0.5*dy);
       if ( coord ==0 && ScanDir[2] ) empty = empty+Form("TMath::Abs(z-%lf)<%lf",z0file+TMath::Nint(Nz/2.)*dz,0.5*dz);
       
       if ( coord ==1 && ScanDir[0] ) empty = empty+Form("TMath::Abs(x-%lf)<%lf",x0file+TMath::Nint(Nx/2.)*dx,0.5*dx);
       if ( coord ==1 && ScanDir[2] ) empty = empty+Form("TMath::Abs(z-%lf)<%lf",z0file+TMath::Nint(Nz/2.)*dz,0.5*dz);
       
       if ( coord ==2 && ScanDir[0] ) empty = empty+Form("TMath::Abs(x-%lf)<%lf",x0file+TMath::Nint(Nx/2.)*dx,0.5*dx);
       if ( coord ==2 && ScanDir[1] ) empty = empty+Form("TMath::Abs(y-%lf)<%lf",y0file+TMath::Nint(Ny/2.)*dy,0.5*dy);
       
       if (qfwhm[coord]>0.03) z0t0 = Getz0t0( (char *) what.Data() , fnm , (char *)empty.Data() ) ;
       else                  z0t0 = Getz0t0( (char *) what.Data() , fnm , (char *)empty.Data() , TString("gaus") ) ; 
       
     #endif

     if ( coord==0 ) { x0   = *z0t0     ;   xl   = *(z0t0+2) ;   xr = *(z0t0+3) ;  thicknessx = *(z0t0+4) ; }
     if ( coord==1 ) { y0   = *z0t0     ;   yl   = *(z0t0+2) ;   yr = *(z0t0+3) ;  thicknessy = *(z0t0+4) ; }
     if ( coord==2 ) { z0   = *z0t0     ;   zl   = *(z0t0+2) ;   zr = *(z0t0+3) ;  thicknessz = *(z0t0+4) ; }

     #if DECONV==0
       //If deconvolution, we keep the already calculated value. Getz0t0 loops on the values in the tree
       //and these are not differentiated
       t0  = *(z0t0+1) ; 
     #else
       t0  = tlavg ;
     #endif

     #if AUTO==0
       z0 = Z0 ; zr = ZR ;  zl = ZL ; t0=T0 ; thicknessz=THICKNESS;
     #endif

     TString cmd="mv Erf.pdf "+TString(dnm)+"plots/"+"Erf_"+TString(bnm)+".pdf" ;
     gSystem->Exec( cmd.Data() );
     
}

//------------------------------------------------------------------------
void  TScan::VoltVsTime_ci ( Int_t coord ) {
   
   TString what  = Form("volt-BlineMean:time-%f"    ,t0 ) ;
   #if I2CORR==1
     what  =  Form("(volt-BlineMean)/(LPower*LPower):time-%f"    ,t0 ) ;
   #endif
   #if DECONV==1
     what  = Form("%5.2f*%7.5f*(volt[Iteration$+1]-volt[Iteration$])/%5.3f+volt-BlineMean:time-%f" , ROSC , Cend , dt, t0 ) ;
   #endif
   
   
   //Find the measured z value falling closest to z0
   Char_t xyz ;
   Double_t dc , dcmm , c0 , c0file ;
   if (coord==0) { xyz='x' ; dc = dx ; dcmm = dx ; c0 = x0 ; c0file = x0file ;}
   if (coord==1) { xyz='y' ; dc = dy ; dcmm = dy ; c0 = y0 ; c0file = y0file ;}
   if (coord==2) { xyz='z' ; dc = dz ; dcmm = dz ; c0 = z0 ; c0file = z0file ;}
   Int_t Ncreal=TMath::Nint((c0-c0file)/dcmm);
   Double_t c0R = c0file + Ncreal * dcmm ;

   Double_t DeltaC = TMath::Nint( (qfwhm[coord]/5)/dcmm )* dcmm;
   Double_t  d1 = c0R + DeltaC, d2 = d1+DeltaC, d3 = d2+DeltaC, d4 = d3+DeltaC , d5=d4+DeltaC ;
   

   TString pthnm = TString(dnm)+"plots/";
   TString pdfnm = pthnm + Form("Vt%ci_",xyz) + TString(bnm) +"_deconv.pdf" , pdf0=pdfnm+"[", pdff=pdfnm+"]";   
   #if DECONV==0
     pdfnm = pthnm + Form("Vt%ci_",xyz) + TString(bnm) +".pdf" , pdf0=pdfnm+"[", pdff=pdfnm+"]";
   #endif
   #if I2CORR==1
     pdfnm = pthnm + Form("Vt%ci_",xyz) + TString(bnm) +"_I2corr.pdf" , pdf0=pdfnm+"[", pdff=pdfnm+"]";
   #endif
   
   
   vector<TString> tit ;
   TString ccut  ;
   ccut= Form("(TMath::Abs(%c-%lf)<%lf",xyz,d1,0.5*dc);
   Int_t tval = TMath::Nint(1000*(d1-c0R)) ; tit.push_back( Form("%d #mum",tval) );
   
   ccut  =  ccut + Form(" || TMath::Abs(%c-%lf)<%lf",xyz,d2,0.5*dc)   ;
   tval = TMath::Nint(1000*(d2-c0R)) ; tit.push_back(Form("%d #mum",tval ));
   
   ccut  =  ccut + Form(" || TMath::Abs(%c-%lf)<%lf",xyz,d3,0.5*dc)   ;
   tval = TMath::Nint(1000*(d3-c0R)) ; tit.push_back(Form("%d #mum",tval));
   
   ccut  =  ccut + Form(" || TMath::Abs(%c-%lf)<%lf",xyz,d4,0.5*dc)   ;
   tval = TMath::Nint(1000*(d4-c0R)) ; tit.push_back(Form("%d #mum",tval));
   
   ccut  =  ccut + Form(" || TMath::Abs(%c-%lf)<%lf",xyz,d5,0.5*dc)   ;
   tval = TMath::Nint(1000*(d5-c0R)) ; tit.push_back(Form("%d #mum",tval));
   
   
   ccut=ccut+TString(")");
   
   //In cases where we have a 2D scan, we need to cut in the other variable   
   
   if ( coord==0 && ScanDir[1] ) ccut=ccut+Form("&& TMath::Abs(y-%lf)<%lf",y0file+TMath::Nint(Ny/2.)*dy,0.5*dy);
   if ( coord==0 && ScanDir[2] ) ccut=ccut+Form("&& TMath::Abs(z-%lf)<%lf",z0file+TMath::Nint(Nz/2.)*dz,0.5*dz);
   
   if ( coord==1 && ScanDir[0] ) ccut=ccut+Form("&& TMath::Abs(x-%lf)<%lf",x0file+TMath::Nint(Nx/2.)*dx,0.5*dx);
   if ( coord==1 && ScanDir[2] ) ccut=ccut+Form("&& TMath::Abs(z-%lf)<%lf",z0file+TMath::Nint(Nz/2.)*dz,0.5*dz);
   
   if ( coord==2 && ScanDir[0] ) ccut=ccut+Form("&& TMath::Abs(x-%lf)<%lf",x0file+TMath::Nint(Nx/2.)*dx,0.5*dx);
   if ( coord==2 && ScanDir[1] ) ccut=ccut+Form("&& TMath::Abs(y-%lf)<%lf",y0file+TMath::Nint(Ny/2.)*dy,0.5*dy);
   
   TCanvas *c1 ;
   TLegend *leg ;
   THStack *hs ;
   TH1D *h1;
   Double_t tVmax=-1000. , tVmin=1000. ;
   for (Int_t iv=0 ; iv<NV ;iv++) {
     TString selection  = ccut + Form(" && Vbias==%d",(Int_t) Vb[iv]) ;
     cout << "Executing plot1D("<<what.Data()<<","<<selection.Data()<<endl ;
     plot1D( (char*) what.Data() , (char*) selection.Data() , 1 , fnm ) ;
     hs = (THStack *) gPad->GetPrimitive("hs") ;
     
     //Protection against cases with step width not matching multiples of 10
     Int_t dval=20.;
     while (hs==0) {
       return;
       TString swhat=Form("==%d ",dval) ;
       dval--;
       TString sby=Form("==%d ",dval) ;
       selection.ReplaceAll( swhat.Data() , sby.Data() );
       std::cout<<"Lowering the zcut threshold zcut="<<dval<<endl;
       plot1D( (char*) what.Data() , (char*) selection.Data() , 1 , fnm ) ;
       hs = (THStack *) gPad->GetPrimitive("hs") ;
     }
     
     vhsV_zi[iv] = hs ;
     TList *hls = (TList *) hs->GetHists()  ;
     for ( Int_t ih = 0 ; ih < hls->GetSize(); ih++ ) {
       h1 = (TH1D * ) hls->At(ih) ;
       cout<<ih<<" " <<tit[ih].Data()<<endl;
       h1->SetTitle(tit[ih]);
     }
     
     //Needed for deconvolution plots, cause Vmax is for non-deconvoluted
     if ( hs->GetMaximum("nostack") > tVmax)  tVmax = hs->GetMaximum("nostack") ;
     if ( hs->GetMinimum("nostack") < tVmin)  tVmin = hs->GetMinimum("nostack") ;
     gPad->Clear();
   }
   #if DECONV==1
     Vmax = tVmax ;
     Vmin = tVmin ;
   #endif
   
   for (Int_t iv=0 ; iv<NV ;iv++) {
     hs = vhsV_zi[iv]  ;
     TString tVbias = Form( "Vbias=%d V" , (Int_t) Vb[iv] ) ;
     c1  = new TCanvas();
     hs->GetXaxis()->SetRangeUser(-2.0,Atrl);
     hs->SetMaximum( tVmax );
     hs->SetMinimum( tVmin );
     hs->Draw("nostack");
     c1->SetGridx(); c1->SetGridy(); 
     leg = c1->BuildLegend( 0.75 , 0.35 , 0.9 , 0.65 , tVbias.Data() );
     leg->SetFillColor( kWhite ) ;
     leg->SetTextSize( 0.033 ) ;
     vlhsV_zi[iv] = leg ;
	
     leg->Draw();
     hs->GetYaxis()->SetTitle( "Signal [V]" ) ;
     hs->GetXaxis()->SetTitle( "Time [ns]"  ) ;
     if (iv==0)    c1->Print( pdf0.Data() ); //It opens the pdf but not writes to it
     c1->Print( pdfnm.Data() ) ;
   }
   c1->Print( pdff.Data() ) ;
   delete c1 ;
     
   
}
//------------------------------------------------------------------------
void  TScan::PlotVar ( Int_t iwhat , Int_t coord ) {
   
   //Plots a single tree variable (tright-tleft,...) versus one of the moving coordinates
   Char_t xyz;
   Double_t c0;
   if (coord==0) { xyz='x' ; c0 = x0  ;}
   if (coord==1) { xyz='y' ; c0 = y0  ;}
   if (coord==2) { xyz='z' ; c0 = z0  ;}

   TString what ;
   if ( iwhat==1) what  = Form("tright-%f:%c-%f" ,t0, xyz , c0 ) ;   
      

   TString pthnm = TString(dnm)+"plots/", pdfnm ;
   if ( iwhat==1) pdfnm = pthnm + Form("tcoll%c_",xyz) + TString(bnm) ;      
      
   //In cases where we have a 2D scan, we need to cut in the other variable   
   TString ccut ;
   if ( coord==0 && ScanDir[1] ) ccut=Form("TMath::Abs(y-%lf)<%lf",y0file+TMath::Nint(Ny/2.)*dy,0.5*dy);
   if ( coord==0 && ScanDir[2] ) ccut=Form("TMath::Abs(z-%lf)<%lf",z0file+TMath::Nint(Nz/2.)*dz,0.5*dz);
   
   if ( coord==1 && ScanDir[0] ) ccut=Form("TMath::Abs(x-%lf)<%lf",x0file+TMath::Nint(Nx/2.)*dx,0.5*dx);
   if ( coord==1 && ScanDir[2] ) ccut=Form("TMath::Abs(z-%lf)<%lf",z0file+TMath::Nint(Nz/2.)*dz,0.5*dz);
   
   if ( coord==2 && ScanDir[0] ) ccut=Form("TMath::Abs(x-%lf)<%lf",x0file+TMath::Nint(Nx/2.)*dx,0.5*dx);
   if ( coord==2 && ScanDir[1] ) ccut=Form("TMath::Abs(y-%lf)<%lf",y0file+TMath::Nint(Ny/2.)*dy,0.5*dy);
   
   cout << "Executing plot1D("<<what.Data()<<","<<ccut.Data()<<endl ;
   plot1D( (char*) what.Data() , (char*) ccut.Data() , 1 , fnm ) ;
   TCanvas *c1 ;
   TLegend *leg ;
   THStack *hs = (THStack *) gPad->GetPrimitive("hs") ;
          
   if (coord==0) vhstcollx.push_back( hs ) ;
   if (coord==1) vhstcolly.push_back( hs ) ;
   if (coord==2) vhstcollz.push_back( hs ) ;
   
   c1  = new TCanvas();
   hs->SetMaximum( TCOLL );
   hs->SetMinimum( 0 );
   hs->Draw("nostack");
   c1->SetGridx(); c1->SetGridy(); 
   leg = c1->BuildLegend( 0.75 , 0.35 , 0.9 , 0.65 );
   leg->SetFillColor( kWhite ) ;
   leg->SetTextSize( 0.033 ) ;
     
   c1->Print( pdfnm+".pdf" );
}
//------------------------------------------------------------------------
void  TScan::VoltVsTime ( ) {
  
   TString what  = Form("volt-BlineMean:time-%f" ,t0 ) ;
   #if I2CORR==1
     what  = Form("(volt-BlineMean)/(LPower*LPower):time-%f" ,t0 ) ;
   #endif

   TString pthnm = TString(dnm)+"plots/";
   TString pdfnm = pthnm + "Vt_"+TString(bnm) +".pdf" , pdf0=pdfnm+"[", pdff=pdfnm+"]";
   #if I2CORR==1
     pdfnm = pthnm + "Vt_"+TString(bnm) +"_I2corr.pdf" , pdf0=pdfnm+"[", pdff=pdfnm+"]";
   #endif
   
   
   TCanvas *c1 ;
   TLegend *leg ;
   THStack *hs ;
   TH1D *h1;
   Double_t tVmax=-2000. , tVmin=2000. ;

   //TString selection=TString("Vbias>-400 && Vbias%20==0");
   TString selection=TString("");
   //TString selection=TString("Vbias%40==0");
   plot1D( (char*) what.Data() , (char*)selection.Data() , 1 , fnm ) ;
   hs = (THStack *) gPad->GetPrimitive("hs") ;
   vhsV_zi[0] = hs  ;
   TList *hls = (TList *) hs->GetHists()  ;
   for ( Int_t ih = 0 ; ih < hls->GetSize(); ih++ ) {
     h1 = (TH1D * ) hls->At(ih) ;
     TString tit = h1->GetTitle();
     Int_t ipos = tit.Index("V=");
     TString stit( tit(ipos+2,7) );
     h1->SetTitle( stit );
   }


   //Needed for deconvolution plots, cause Vmax is for non-deconvoluted
   if ( hs->GetMaximum("nostack") > tVmax)  tVmax = hs->GetMaximum("nostack") ;
   if ( hs->GetMinimum("nostack") < tVmin)  tVmin = hs->GetMinimum("nostack") ;
   
   hs = vhsV_zi[0]  ;
   c1  = new TCanvas();
   
   hs->GetXaxis()->SetRangeUser(-2.0,Atrl);
   //hs->GetXaxis()->SetRangeUser(-2.0,0.75*(hs->GetXaxis()->GetXmax()));
   
   hs->SetMaximum( Vmax );
   hs->SetMinimum( Vmin-0.05*(Vmax-Vmin) );
   hs->Draw("nostack");
   c1->SetGridx(); c1->SetGridy(); 
   leg = c1->BuildLegend( 0.9 , 0.1 , 0.99 , 0.78  );
   leg->SetFillColor( kWhite ) ;
   leg->SetTextSize( 0.033 ) ;
   vlhsV_zi[0] = leg ;

   leg->Draw();
   hs->GetYaxis()->SetTitle( "Signal [V]" ) ;
   #if I2CORR==1
     hs->GetYaxis()->SetTitle( "Signal [I^{2} corr, a.u.]" ) ;
   #endif
  
   hs->GetXaxis()->SetTitle( "Time [ns]"  ) ;
   c1->Print( pdfnm.Data() ) ;
   delete c1 ;
   
}

//------------------------------------------------------------------------
void TScan::TimeCoordVolt2D( Int_t coord ) {

   TH2D *h2 ;
   TCanvas *c2 ;
   vector<TText *> tv ; 
   TText *text ;
   TString pthnm = TString(dnm)+"plots/";
   TString pdfnm = pthnm+"V2D_" + TString(bnm)  ;
   #if DECONV==1
     pdfnm = pdfnm+TString("_deconv");
   #endif

   #if I2CORR==1
     pdfnm = pdfnm+TString("_I2corr");
   #endif
   
   if (HIGHQPDF == 1) pdfnm = pdfnm + ".pdf" ; 
   else               pdfnm = pdfnm + ".png" ;  
   TString pdf0=pdfnm+"[" ; 
   TString pdff=pdfnm+"]";
   
   gStyle->SetPadRightMargin(0.25);   //Distance between plot and pad
   gStyle->SetPadBottomMargin(0.09);   //Distance between plot and pad
   
   //Get the histograms
   Char_t xyz ; 
   Double_t val0;
   if (coord==0) { xyz='x'; val0=x0 ;} 
   if (coord==1) { xyz='y'; val0=y0 ;} 
   if (coord==2) { xyz='z'; val0=z0 ;} 
   for ( Int_t iv = 0 ; iv < NV ; iv++ ) {
   
     TString selection  = Form("Vbias==%d",(Int_t) Vb[iv]) ;
 
     //If this is a 2D scan, then go to the center of the other coordinate
     if ( coord==0 && ScanDir[1] ) selection=selection+Form("&& TMath::Abs(y-%lf)<%lf",y0file+TMath::Nint(Ny/2.)*dy,0.5*dy);
     if ( coord==0 && ScanDir[2] ) selection=selection+Form("&& TMath::Abs(z-%lf)<%lf",z0file+TMath::Nint(Nz/2.)*dz,0.5*dz);
     
     if ( coord==1 && ScanDir[0] ) selection=selection+Form("&& TMath::Abs(x-%lf)<%lf",x0file+TMath::Nint(Nx/2.)*dx,0.5*dx); 
     if ( coord==1 && ScanDir[2] ) selection=selection+Form("&& TMath::Abs(z-%lf)<%lf",z0file+TMath::Nint(Nz/2.)*dz,0.5*dz); 
     
     if ( coord==2 && ScanDir[0] ) selection=selection+Form("&& TMath::Abs(x-%lf)<%lf",x0file+TMath::Nint(Nx/2.)*dx,0.5*dx); 
     if ( coord==2 && ScanDir[1] ) selection=selection+Form("&& TMath::Abs(y-%lf)<%lf",y0file+TMath::Nint(Ny/2.)*dy,0.5*dy); 
     
     TString what  = Form("time-%f:%c-%f:(volt-BlineMean)", t0,xyz,val0) ;
     #if I2CORR==1
       what = Form("time-%f:%c-%f:(volt-BlineMean)/(LPower*LPower)", t0,xyz,val0) ;
     #endif
     
     #if DECONV==1
       what  = Form("time-%f:z-%f:%5.2f*%f*(volt[Iteration$+1]-volt[Iteration$])/%f+volt-BlineMean" , t0,z0, ROSC , Cend , dt ) ;
     #endif
     plot2D( (char*) what.Data() , (char*) selection.Data() , fnm ) ;
     
     what  = Form("(volt-BlineMean)=f(time-%f,%c-%f)",t0, xyz, val0) ;
     #if DECONV==1
       what  = Form("%5.2f*%f*(volt[Iteration$+1]-volt[Iteration$])/%f+volt-BlineMean=f(time-%f,z-%f)", ROSC , Cend , dt,t0, z0) ;
     #endif
     #if I2CORR==1
       what  = Form("(volt-BlineMean)/(LPower*LPower)=f(time-%f,%c-%f)",t0, xyz,val0) ;
     #endif
     
     //h2 = (TH2D *) gDirectory->Get(Form("volt-BlineMean=f(time-40.059797,z)") ;
     h2 = (TH2D *) gPad->GetPrimitive( what.Data() ) ;
     if (Atrl+2<30.) h2->GetXaxis()->SetRangeUser(-2.0,Atrl+2.0) ;
     else            h2->GetXaxis()->SetRangeUser(-2.0,30.) ;
     Double_t thickness ;
     if (coord==0) thickness = thicknessx ;
     if (coord==1) thickness = thicknessy ;
     if (coord==2) thickness = thicknessz ;
     h2->GetYaxis()->SetRangeUser(val0-0.01,val0+thickness+0.01) ;
     h2->GetXaxis()->SetTitle("Time [ns]") ;
     h2->GetYaxis()->SetTitle( Form("%c [mm]",xyz) ) ;
     TPaletteAxis *palette = (TPaletteAxis*)h2->GetListOfFunctions()->FindObject("palette");
     palette->GetAxis()->SetTitleOffset(1.1); palette->GetAxis()->SetTitleSize(0.045);
     palette->GetAxis()->SetTitle( "Signal [V]" );
     #if I2CORR==1
       palette->GetAxis()->SetTitleOffset(1.25); palette->GetAxis()->SetTitleSize(0.045);
       palette->GetAxis()->SetTitle( "Signal [I^{2} corr., a.u.]" );
     #endif
     if (coord==0) vhvtx.push_back(h2) ;     
     if (coord==1) vhvty.push_back(h2) ;     
     if (coord==2) vhvtz.push_back(h2) ;     

     selection=Form("%d V",(Int_t) Vb[iv]) ;
     if (Atrl+2<30.) text=new TText( 0.8*Atrl, 0.8*thickness , selection );
     else            text=new TText( 0.8*30. , 0.8*thickness , selection );
     tv.push_back(text) ;
     cout << endl; cout << endl;

   }
   
   //Plot them (needed in case we wanted to have subpads
   if (PORTRAIT==1)  c2 = new TCanvas("c2","c2",400,600);
   else              c2 = new TCanvas("c2","c2",600,400);
   
   c2->Divide(NPX,NPY) ;
   
   c2->Print( pdf0.Data() );
   UInt_t iDraw = 0 ;
   for ( Int_t iv = 0 ; iv < NV ; iv++ ) {
     iDraw++ ;
     if (coord==0) h2 = (TH2D*) vhvtx[iv] ; 
     if (coord==1) h2 = (TH2D*) vhvty[iv] ; 
     if (coord==2) h2 = (TH2D*) vhvtz[iv] ; 
     c2->cd(iDraw);
     
     if ( NPX*NPY != 1) {
       if ( iDraw <= NPX*(NPY-1) ) gPad->SetBottomMargin(1.e-5);
       if ( iDraw > NPX )          gPad->SetTopMargin(1.e-5);
       if ( iDraw > NPX*(NPY-1) )  gPad->SetBottomMargin(0.1);
     }     
     
     h2->Draw("colz") ;
     text = (TText *) tv[iv];
     text->Draw();
     if (iDraw==NPX*NPY) { 
       c2->Print( pdfnm.Data() );
       iDraw=0;
       c2->Clear();
       c2->Divide(NPX,NPY) ;
     }
     
   }  
    
   if ( iDraw!=0 ) c2->Print( pdfnm.Data() );
   
   c2->Print( pdff.Data() );
   
   //Convert png to pdf
   if  (HIGHQPDF != 1) {
     TString cmd="convert -trim -quality 100 -rotate 90 "+pdfnm+" kk.pdf"; 
     gSystem->Exec( cmd.Data() );
     cmd ="mv kk.pdf "+TString(dnm)+"plots/"+"V2D_"+TString(bnm)+".pdf" ;
     gSystem->Exec( cmd.Data() );
     cmd ="rm -f "+TString(pdfnm);
     gSystem->Exec( cmd.Data() );
   }
      
   tv.clear();
   delete c2 ;  
   delete h2 ; 
}

//------------------------------------------------------------------------
void TScan::Maps2D( Int_t iwhat ) {

   //Plots 2D maps of variables (Q, tcoll...). 
   
   TString Plot;
   if ( iwhat == 0) Plot=TString("Q2D");
   if ( iwhat == 1) Plot=TString("tcoll2D");
   if ( iwhat == 2) Plot=TString("vd2D");
   cout <<  "Plotting " << Plot.Data()<<" maps" <<endl ;
   TString myfnm = TString(dnm)+"histos/"+Plot+"_" + TString(bnm)+".hroot"  ;
   
   TCanvas *c2 ;
   vector<TString> tv ; 
   TString pthnm = TString(dnm)+"plots/";
   TString pdfnm = pthnm+ Plot +"_" + TString(bnm)  ;
   #if DECONV==1
     pdfnm = pdfnm + TString("_deconv");
   #endif
   #if I2CORR==1
     if ( Plot.EqualTo("Q2D") ) pdfnm = pdfnm + TString("_I2corr");
   #endif
   
   if (HIGHQPDF == 1) pdfnm = pdfnm + ".pdf" ; 
   else               pdfnm = pdfnm + ".png" ;  
   TString pdf0=pdfnm+"[" ; 
   TString pdff=pdfnm+"]";
   
   gStyle->SetPadRightMargin(0.25);   //Distance between plot and pad
   gStyle->SetPadBottomMargin(0.09);   //Distance between plot and pad
   
   /* Get the histograms */
   //Build the variable we want to plot

   //Q2D and vd2D are only different by the integration time
   Double_t TheAtrl ;
   if ( Plot.EqualTo("Q2D") )  TheAtrl=Atrl;
   if ( Plot.EqualTo("vd2D") ) TheAtrl=0.4;

   TString what , twhat ;
   if ( Plot.EqualTo("Q2D") || Plot.EqualTo("vd2D")    ) what= Form( "Sum$((volt-BlineMean)*(time>%f &&time<%f))",t0,t0+TheAtrl );
   if ( Plot.EqualTo("tcoll2D") ) what= TString( Form("tright-%f",tlavg ) );

   if ( ScanDir[0] && ScanDir[1] ) what  = Form("x-%f:y-%f:", x0, y0 ) + what ;
   if ( ScanDir[0] && ScanDir[2] ) what  = Form("x-%f:z-%f:", x0, z0 ) + what ;
   if ( ScanDir[1] && ScanDir[2] ) what  = Form("y-%f:z-%f:", y0, z0 ) + what ;
   #if I2CORR==1
     if ( Plot.EqualTo("Q2D") || Plot.EqualTo("vd2D") )      { 
       what  = Form( "Sum$((volt-BlineMean)*(time>%f &&time<%f))/(LPower*LPower)",t0,t0+TheAtrl ) ;
       if ( ScanDir[0] && ScanDir[1] ) what  = Form("x-%f:y-%f:", x0, y0 ) + what ;
       if ( ScanDir[0] && ScanDir[2] ) what  = Form("x-%f:z-%f:", x0, z0 ) + what ;
       if ( ScanDir[1] && ScanDir[2] ) what  = Form("y-%f:z-%f:", y0, z0 ) + what ;
     }
   #endif


   //Titles of the histograms that are going to be created
   if ( Plot.EqualTo("Q2D") || Plot.EqualTo("vd2D") ) {
     if ( ScanDir[0] && ScanDir[1]) twhat  = Form("Sum$((volt-BlineMean)*(time>%f &&time<%f))=f(x-%f,y-%f)",t0, t0+TheAtrl,x0,y0) ;
     if ( ScanDir[0] && ScanDir[2]) twhat  = Form("Sum$((volt-BlineMean)*(time>%f &&time<%f))=f(x-%f,z-%f)",t0, t0+TheAtrl,x0,z0) ;
     if ( ScanDir[1] && ScanDir[2]) twhat  = Form("Sum$((volt-BlineMean)*(time>%f &&time<%f))=f(y-%f,z-%f)",t0, t0+TheAtrl,y0,z0) ;
     #if I2CORR==1
       if ( ScanDir[0] && ScanDir[1]) twhat  = Form("Sum$((volt-BlineMean)*(time>%f &&time<%f))/(LPower*LPower)=f(x-%f,y-%f)",t0, t0+TheAtrl,x0,y0) ;
       if ( ScanDir[0] && ScanDir[2]) twhat  = Form("Sum$((volt-BlineMean)*(time>%f &&time<%f))/(LPower*LPower)=f(x-%f,z-%f)",t0, t0+TheAtrl,x0,z0) ;
       if ( ScanDir[1] && ScanDir[2]) twhat  = Form("Sum$((volt-BlineMean)*(time>%f &&time<%f))/(LPower*LPower)=f(y-%f,z-%f)",t0, t0+TheAtrl,y0,z0) ;
     #endif
   } else if ( Plot.EqualTo("tcoll2D" ) ) {
     if ( ScanDir[0] && ScanDir[1]) twhat  = Form("tright-%f=f(x-%f,y-%f)",tlavg,x0,y0) ;
     if ( ScanDir[0] && ScanDir[2]) twhat  = Form("tright-%f=f(x-%f,z-%f)",tlavg,x0,z0) ;
     if ( ScanDir[1] && ScanDir[2]) twhat  = Form("tright-%f=f(y-%f,z-%f)",tlavg,y0,z0) ;     
   }


   TH2D *h2 ;
   for ( Int_t iv = 0 ; iv < NV ; iv++ ) {
   
     TString selection  = Form("Vbias==%d",(Int_t) Vb[iv]) ;
     
     plot2D( (char*) what.Data() , (char*) selection.Data() , fnm ) ;
          
     h2 = (TH2D *) gPad->GetPrimitive( twhat.Data() ) ; 
     

     if ( ScanDir[0] && ScanDir[1] ) { h2->GetXaxis()->SetTitle("X [mm]") ; h2->GetYaxis()->SetTitle("Y [mm]") ; }
     if ( ScanDir[0] && ScanDir[2] ) { h2->GetXaxis()->SetTitle("X [mm]") ; h2->GetYaxis()->SetTitle("Z [mm]") ; }
     if ( ScanDir[1] && ScanDir[2] ) { h2->GetXaxis()->SetTitle("Y [mm]") ; h2->GetYaxis()->SetTitle("Z [mm]") ; }
     TPaletteAxis *palette = (TPaletteAxis*)h2->GetListOfFunctions()->FindObject("palette");
     palette->GetAxis()->SetTitleOffset(1.1); palette->GetAxis()->SetTitleSize(0.045);
     
     if ( Plot.EqualTo("Q2D" ) || Plot.EqualTo("vd2D") ) { 
       if ( Plot.EqualTo("Q2D") ) {
         palette->GetAxis()->SetTitle( Form("Charge in %d ns [a.u.]",TMath::Nint(TheAtrl)) );
	 #if I2CORR==1
	   palette->GetAxis()->SetTitle( Form("Charge in %d ns [I^{2} corr., a.u.]",TMath::Nint(TheAtrl)) );
	 #endif
       } 
       if ( Plot.EqualTo("vd2D") ) {
         palette->GetAxis()->SetTitle( "Drift velocity [a.u.]" );
	 #if I2CORR==1
	   palette->GetAxis()->SetTitle( "Drift velocity [I^{2} corr., a.u.]" );
	 #endif
       }
       
       //Rename the histogram
       TString oldname= h2->GetName();
       if ( Plot.EqualTo("Q2D" ) ) oldname.ReplaceAll("Sum$((volt-BlineMean)*","Q");
       if ( Plot.EqualTo("vd2D") ) oldname.ReplaceAll("Sum$((volt-BlineMean)*","vd");
       oldname.ReplaceAll("))/",")/");
       h2->SetNameTitle(oldname,oldname);
              
       if ( Plot.EqualTo("Q2D" ) )  vhQ2D[iv]=h2 ;    
       if ( Plot.EqualTo("vd2D") ) vhvd2D[iv]=h2 ;    
     }
     
     if ( Plot.EqualTo("tcoll2D" ) ) {
       palette->GetAxis()->SetTitle( Form("Collection Time [ns]") );
       vhtcoll2D[iv]=h2 ;    
     }

     selection=Form("%d V",(Int_t) Vb[iv]) ;
     
     TLatex text = TLatex();
     text.SetNDC();
     text.DrawLatex(0.17, 0.8, selection.Data() );

     tv.push_back(selection) ;
     cout << endl; cout << endl;

   }
   
   //MULTIPLY DRIFT VELOCITY BY WEIGHTING FIELD
   #if WEIGHTINGFIELD>0
     if ( Plot.EqualTo("vd2D" ) ) vd_times_Ew() ;
   #endif
   
      
   //Write histos to file
   TFile *f2Dout=new TFile(myfnm.Data(),"RECREATE");
   f2Dout->cd();
   for ( Int_t iv = 0 ; iv < NV ; iv++ ) {
     if ( Plot.EqualTo("Q2D" ) )     h2 = (TH2D*) vhQ2D[iv] ; 
     if ( Plot.EqualTo("vd2D" ) )    {
       h2   = (TH2D*) vhvd2D[iv] ; 
       #if WEIGHTINGFIELD>0
         TH2D *h2Ew = (TH2D*) vhvdEw2D[iv] ; h2Ew->Write();
       #endif
     }
     if ( Plot.EqualTo("tcoll2D" ) ) h2 = (TH2D*) vhtcoll2D[iv] ; 
     h2->Write();
   }
   f2Dout->Write();
   f2Dout->Close();
   
     
   //Plot them (needed in case we wanted to have subpads
   if (PORTRAIT==1)  c2 = new TCanvas("c2","c2",400,600);
   else              c2 = new TCanvas("c2","c2",600,400);
   
   if (NV>1) c2->Divide(NPX,NPY) ;
   
   
   c2->Print( pdf0.Data() );
   UInt_t iDraw = 0 ;
   for ( Int_t iv = 0 ; iv < NV ; iv++ ) {
     iDraw++ ;
     if ( Plot.EqualTo("Q2D" ) )     h2 = (TH2D*) vhQ2D[iv] ; 
     if ( Plot.EqualTo("vd2D" ) )    h2 = (TH2D*) vhvd2D[iv] ; 
     if ( Plot.EqualTo("tcoll2D" ) ) { 
       h2 = (TH2D*) vhtcoll2D[iv] ; 
       h2->SetMinimum(0); h2->SetMaximum(TCOLL);
     }
     c2->cd(iDraw);
     
     if ( NPX*NPY != 1) {
       if ( iDraw <= NPX*(NPY-1) ) gPad->SetBottomMargin(1.e-5);
       if ( iDraw > NPX )          gPad->SetTopMargin(1.e-5);
       if ( iDraw > NPX*(NPY-1) )  gPad->SetBottomMargin(0.1);
     }     
     
     h2->Draw("colz") ;
     
     TString selval = tv[iv];
     TLatex text = TLatex() ;  text.SetNDC() ; text.DrawLatex(0.17, 0.8, selval.Data() );
     if (iDraw==NPX*NPY) { 
       c2->Print( pdfnm.Data() );
       iDraw=0;
       c2->Clear();
       c2->Divide(NPX,NPY) ;
     }
     
   }  
    
   if ( iDraw!=0 ) c2->Print( pdfnm.Data() );
   
   c2->Print( pdff.Data() );
   
   //Convert png to pdf
   if  (HIGHQPDF != 1) {
     TString cmd="convert -trim -quality 100 -rotate 90 "+pdfnm+" kk.pdf"; 
     gSystem->Exec( cmd.Data() );
     cmd ="mv kk.pdf "+TString(dnm)+"plots/"+"Q2D_"+TString(bnm)+".pdf" ;
     gSystem->Exec( cmd.Data() );
     cmd ="rm -f "+TString(pdfnm);
     gSystem->Exec( cmd.Data() );
   }
      
   tv.clear();

   delete c2 ; 
   delete h2 ; 
}
//------------------------------------------------------------------------
void TScan::vd_times_Ew() {   
     
     //For a 2D scan, multiplies drift velocity (calculated withtout Ew)
     //by an external mapping of Ew.
     
     //Open the Weighting Field (careful, this histogram comes in um)
     TFile *fEw=new TFile("/home/ssdstorage/TPA/data/Feb2016_HVCMOSv3_unirrad/HVCMOS_TRACS_WeightingField_w100um_d150um.root");
     TH2D *hEw = (TH2D*) fEw->Get("hwf");

     Double_t XPN=0.063 , YPN = 0.0188 ; // Measured point for center of DNW implant(XPN,YPN) [Simulation center=(250,0)]
          
     for ( Int_t iv = 0 ; iv < NV ; iv++ ) {

        //Create a new measured histogram hvEw that starts at Y=YPN, 
	//cause below that we do not have info on Ew
        TH2D* h2 = (TH2D*) vhvd2D[iv] ;
	
	//Find the bin where YPN happens, then create a new histogram
	TAxis *xaxm = h2->GetXaxis() ; 
	TAxis *yaxm = h2->GetYaxis() ; 
	Double_t iby = yaxm->FindBin( YPN );
	
	Int_t Ny       = h2->GetNbinsY() , Nyp = Ny - iby + 1 ;
	Double_t Ay    = yaxm->GetBinCenter(2) - yaxm->GetBinCenter(1) ;
	Double_t Yminp = yaxm->GetBinCenter(iby) - 0.5*Ay ;
	Double_t Ymax  = h2->GetYaxis()->GetXmax();

	Int_t    Nx   = h2->GetNbinsX() ;
	Double_t Xmin = h2->GetXaxis()->GetXmin() , Xmax = h2->GetXaxis()->GetXmax()  ;
	
	TString newtit = TString("Ew ") + h2->GetTitle() ;
	TH2D *hvEw=new TH2D( newtit,newtit,Nx,Xmin,Xmax,Nyp,Yminp,Ymax );
	

	for ( Int_t ix=1 ; ix<=Nx ; ix++ ) {
	  
	  for ( Int_t iy=1 ; iy<=Ny ; iy++ ) {

	     //Measured drift velocity
	     Int_t gbin = h2->GetBin(ix,iy) ;
	     Double_t vdmeas = h2->GetBinContent( gbin ) ;
	     Double_t xmeas = xaxm->GetBinCenter( ix );
	     Double_t ymeas = yaxm->GetBinCenter( iy );
             
	     if ( ymeas >= YPN ) {
	       //Simulated Weighting Field
	       Double_t xsim = 250  + 1000.*(xmeas-XPN) ;
	       Double_t ysim =   0. + 1000.*(ymeas-YPN) ;
	       Double_t   Ew = hEw->GetBinContent( hEw->FindBin( xsim , ysim ) ) ;
	       std::cout <<  "(x,y,Ew)=(" << xsim <<"," <<ysim << ", "<<Ew << ")"<<std::endl ;
	       
	       //New drift velocity
	       if (Ew!=0.) hvEw->Fill( xmeas , ymeas , vdmeas / Ew );
	     }
	     
	  } //Loop on Y
	  
	}   //Loop on X
	
	//Overwrite current drift velocity
	vhvdEw2D[iv] = hvEw ;

     }      //Loop on voltages
     
     fEw->Close();
     
}
//------------------------------------------------------------------------
void TScan::Q2D_IntegrationTime( ) {

   //Plots 2D maps of variables (Q, tcoll...). In case of charge maps, 
   //do several plots, increasing the collection time. Ex: Q(1ns),Q(2ns),...Q(Atrl)
   
   TString Plot=TString("Q2D");
   Int_t iwhat = 0 ; 
   cout <<  "Plotting " << Plot.Data()<<" maps" <<endl ;
   TString myfnm = TString(dnm)+"histos/"+Plot+"_" + TString(bnm)+".hroot"  ;
   
   TCanvas *c2 ;
   TString pthnm = TString(dnm)+"plots/";
   TString pdfnm = pthnm+ Plot +"_" + TString(bnm)  ;
   #if DECONV==1
     pdfnm = pdfnm + TString("_deconv");
   #endif
   #if I2CORR==1
     if ( Plot.EqualTo("Q2D") ) pdfnm = pdfnm + TString("_I2corr");
   #endif
   
   if (HIGHQPDF == 1) pdfnm = pdfnm + ".pdf" ; 
   else               pdfnm = pdfnm + ".png" ;  
   TString pdf0=pdfnm+"[" ; 
   TString pdff=pdfnm+"]";
   
   gStyle->SetPadRightMargin(0.25);   //Distance between plot and pad
   gStyle->SetPadBottomMargin(0.09);   //Distance between plot and pad
   
   /* Get the histograms */
   //Build the variable we want to plot

   TFile *f2Dout=new TFile(myfnm.Data(),"RECREATE");
   
   if (PORTRAIT==1)  c2 = new TCanvas("c2","c2",400,600);
   else              c2 = new TCanvas("c2","c2",600,400);
   c2->cd();
   c2->Print( pdf0.Data() );

   for (Int_t it=1;it<=Atrl;it++) {

     TString what , twhat ;
     if ( Plot.EqualTo("Q2D")     ) what= Form( "Sum$((volt-BlineMean)*(time>%f &&time<%f))",t0,t0+it );

     if ( ScanDir[0] && ScanDir[1] ) what  = Form("x-%f:y-%f:", x0, y0 ) + what ;
     if ( ScanDir[0] && ScanDir[2] ) what  = Form("x-%f:z-%f:", x0, z0 ) + what ;
     if ( ScanDir[1] && ScanDir[2] ) what  = Form("y-%f:z-%f:", y0, z0 ) + what ;
     #if I2CORR==1
       if ( Plot.EqualTo("Q2D") )      { 
	 what  = Form( "Sum$((volt-BlineMean)*(time>%f &&time<%f))/(LPower*LPower)",t0,t0+it ) ;
	 if ( ScanDir[0] && ScanDir[1] ) what  = Form("x-%f:y-%f:", x0, y0 ) + what ;
	 if ( ScanDir[0] && ScanDir[2] ) what  = Form("x-%f:z-%f:", x0, z0 ) + what ;
	 if ( ScanDir[1] && ScanDir[2] ) what  = Form("y-%f:z-%f:", y0, z0 ) + what ;
       }
     #endif


     //Titles of the histograms that are going to be created
     if ( Plot.EqualTo("Q2D") ) {
       if ( ScanDir[0] && ScanDir[1]) twhat  = Form("Sum$((volt-BlineMean)*(time>%f &&time<%f))=f(x-%f,y-%f)",t0, t0+it,x0,y0) ;
       if ( ScanDir[0] && ScanDir[2]) twhat  = Form("Sum$((volt-BlineMean)*(time>%f &&time<%f))=f(x-%f,z-%f)",t0, t0+it,x0,z0) ;
       if ( ScanDir[1] && ScanDir[2]) twhat  = Form("Sum$((volt-BlineMean)*(time>%f &&time<%f))=f(y-%f,z-%f)",t0, t0+it,y0,z0) ;
       #if I2CORR==1
	 if ( ScanDir[0] && ScanDir[1]) twhat  = Form("Sum$((volt-BlineMean)*(time>%f &&time<%f))/(LPower*LPower)=f(x-%f,y-%f)",t0, t0+it,x0,y0) ;
	 if ( ScanDir[0] && ScanDir[2]) twhat  = Form("Sum$((volt-BlineMean)*(time>%f &&time<%f))/(LPower*LPower)=f(x-%f,z-%f)",t0, t0+it,x0,z0) ;
	 if ( ScanDir[1] && ScanDir[2]) twhat  = Form("Sum$((volt-BlineMean)*(time>%f &&time<%f))/(LPower*LPower)=f(y-%f,z-%f)",t0, t0+it,y0,z0) ;
       #endif
     }


     TH2D *h2 ;
     for ( Int_t iv = 0 ; iv < NV ; iv++ ) {

       TString selection  = Form("Vbias==%d",(Int_t) Vb[iv]) ;

       plot2D( (char*) what.Data() , (char*) selection.Data() , fnm ) ;

       h2 = (TH2D *) gPad->GetPrimitive( twhat.Data() ) ; 


       if ( ScanDir[0] && ScanDir[1] ) { h2->GetXaxis()->SetTitle("X [mm]") ; h2->GetYaxis()->SetTitle("Y [mm]") ; }
       if ( ScanDir[0] && ScanDir[2] ) { h2->GetXaxis()->SetTitle("X [mm]") ; h2->GetYaxis()->SetTitle("Z [mm]") ; }
       if ( ScanDir[1] && ScanDir[2] ) { h2->GetXaxis()->SetTitle("Y [mm]") ; h2->GetYaxis()->SetTitle("Z [mm]") ; }
       TPaletteAxis *palette = (TPaletteAxis*)h2->GetListOfFunctions()->FindObject("palette");
       palette->GetAxis()->SetTitleOffset(1.1); palette->GetAxis()->SetTitleSize(0.045);

       if ( Plot.EqualTo("Q2D" ) ) { 
	 palette->GetAxis()->SetTitle( Form("Charge in %d ns [a.u.]",it) );
	 #if I2CORR==1
	   palette->GetAxis()->SetTitle( Form("Charge in %d ns [I^{2} corr., a.u.]",it) );
	 #endif

	 //Rename the histogram
	 TString oldname= h2->GetName();
	 oldname.ReplaceAll("Sum$((volt-BlineMean)*","Q");
	 oldname.ReplaceAll("))/",")/");
	 h2->SetNameTitle(oldname,oldname);
	 selection=Form("%d V",(Int_t) Vb[iv]) ;


	 TLatex text = TLatex();
	 text.SetNDC();
	 text.DrawLatex(0.17, 0.8, selection.Data() );
         
	 c2->cd();
	 h2->Draw("colz");
         c2->Print( pdfnm.Data() );
        
	 f2Dout->cd() ; h2->Write() ; 
	 
	 if (it==TMath::Nint(Atrl)) vhQ2D.push_back(h2) ;    
	 
       }


   
     }
     delete h2;

   }

   c2->Print( pdff.Data() );

   delete c2 ; 

}
//------------------------------------------------------------------------
void TScan::TimeVbiasVolt2D() {

   TH2D *h2 ;
   TCanvas *c2 ;
   TString pthnm = TString(dnm)+"plots/";
   TString pdfnm = pthnm+"V2D_" + TString(bnm)  ;
   #if I2CORR==1
     pdfnm = pdfnm + "_I2corr" ; 
   #endif
   
   if (HIGHQPDF == 1) pdfnm = pdfnm + ".pdf" ; 
   else               pdfnm = pdfnm + ".png" ;  
   TString pdf0=pdfnm+"[" ; 
   TString pdff=pdfnm+"]";
   
   gStyle->SetPadRightMargin(0.25);   //Distance between plot and pad
   gStyle->SetPadBottomMargin(0.09);   //Distance between plot and pad
   
   //Get the histograms
   
   TString what  = Form("time-%f:Vbias:volt-BlineMean", t0) ;
   #if I2CORR==1
     what  = Form("time-%f:Vbias:(volt-BlineMean)/(LPower*LPower)", t0) ;
   #endif
   
   //TString selection=TString("Vbias>-300");
   TString selection=TString("");
   plot2D( (char*) what.Data() , (char*) selection.Data() , fnm ) ;
   what  = Form("volt-BlineMean=f(time-%f,Vbias)",t0) ;
   #if I2CORR==1
     what  = Form("(volt-BlineMean)/(LPower*LPower)=f(time-%f,Vbias)",t0) ;
   #endif

   h2 = (TH2D *) gPad->GetPrimitive( what.Data() ) ;
   if (Atrl+2<30.) h2->GetXaxis()->SetRangeUser(-2.0,Atrl+2.0) ;
   else            h2->GetXaxis()->SetRangeUser(-2.0,30.) ;
   h2->GetXaxis()->SetTitle("Time [ns]") ;
   h2->GetYaxis()->SetTitle("Bias voltage [V]") ;
   TPaletteAxis *palette = (TPaletteAxis*)h2->GetListOfFunctions()->FindObject("palette");
   palette->GetAxis()->SetTitleOffset(1.1); palette->GetAxis()->SetTitleSize(0.045);
   palette->GetAxis()->SetTitle( "Signal [V]" );
   #if I2CORR==1
     palette->GetAxis()->SetTitle( "Signal [I^{2} corr., a.u.]" );
   #endif

   cout << endl; cout << endl;

   
   //Plot them (needed in case we wanted to have subpads
   if (PORTRAIT==1)  c2 = new TCanvas("c2","c2",400,600);
   else              c2 = new TCanvas("c2","c2",600,400);
   
   
     
   h2->Draw("colz") ;
   c2->Print( pdfnm.Data() );
      
   //Convert png to pdf
   if  (HIGHQPDF != 1) {
     TString cmd="convert -trim -quality 100 "+pdfnm+" kk.pdf"; 
     gSystem->Exec( cmd.Data() );
     cmd ="mv kk.pdf "+TString(dnm)+"plots/"+"V2D_"+TString(bnm)+".pdf" ;
     gSystem->Exec( cmd.Data() );
     cmd ="rm -f "+TString(pdfnm);
     gSystem->Exec( cmd.Data() );
   }
      
   delete c2 ;  
   delete h2 ;
   
}
//------------------------------------------------------------------------
void TScan::Vdt( Double_t toffset , char *cond , Int_t c ) {
   
   TString xyz, Cond ;
   Double_t cl , cr , c0 ;
   if ( c==0 ) { xyz=TString("x") ; cl=xl ; cr=xr; c0=x0; }
   if ( c==1 ) { xyz=TString("y") ; cl=yl ; cr=yr; c0=y0; }
   if ( c==2 ) { xyz=TString("z") ; cl=zl ; cr=zr; c0=z0; }
   
   //TString zsel = Form("%s>%f && %s<%f", xyz.Data(), cl , xyz.Data(), cr ) ;
   TString zsel = "" ;
   if ( xyz =="x" && ScanDir[1] ) Cond = Form("TMath::Abs(y-%lf)<%lf",y0file+TMath::Nint(Ny/2.)*dy,0.5*dy);
   if ( xyz =="x" && ScanDir[2] ) Cond = Form("TMath::Abs(z-%lf)<%lf",z0file+TMath::Nint(Nz/2.)*dz,0.5*dz);

   if ( xyz =="y" && ScanDir[0] ) Cond = Form("TMath::Abs(x-%lf)<%lf",x0file+TMath::Nint(Nx/2.)*dx,0.5*dx);
   if ( xyz =="y" && ScanDir[2] ) Cond = Form("TMath::Abs(z-%lf)<%lf",z0file+TMath::Nint(Nz/2.)*dz,0.5*dz);

   if ( xyz =="z" && ScanDir[0] ) Cond = Form("TMath::Abs(x-%lf)<%lf",x0file+TMath::Nint(Nx/2.)*dx,0.5*dx);
   if ( xyz =="z" && ScanDir[1] ) Cond = Form("TMath::Abs(y-%lf)<%lf",y0file+TMath::Nint(Ny/2.)*dy,0.5*dy);

   //Calculate Vdt or efficiency (depends on toffset)
   #if DECONV == 0
     if (toffset == 0.) {
       cout <<  "Calculating DRIFT VELOCITY" <<endl ;
       plotvd( xyz ,(char *) zsel.Data() , t0 , c0 , BLMean, (char *)Cond.Data() , fnm , ATRISE, I2CORR);
     } else {
       cout <<  "Calculating EFFICIENCY" <<endl ;
       plotvd( xyz ,(char *) zsel.Data() , t0 , c0 , BLMean, (char *)Cond.Data() , fnm , toffset , I2CORR);
     }
   #else
     if (toffset == 0.) {
       cout <<  "Calculating DRIFT VELOCITY" <<endl ;
       plotvd( xyz ,(char *) zsel.Data() , t0 , c0 , Cend , dt, BLMean, (char *)Cond.Data() , fnm , ATRISE, I2CORR);
     } else {
       cout <<  "Calculating EFFICIENCY" <<endl ;
       //plotvd( (char *) zsel.Data() , t0 , c0 , Cend , dt, BLMean, cond , fnm , toffset , I2CORR);  //Q(z) using deconv pulses
       plotvd( xyz ,(char *) zsel.Data() , t0 , c0 , BLMean, (char *)Cond.Data() , fnm , toffset , I2CORR);               //Q(z) using std pulses
     }

   #endif
   
   
   //Legends should only display voltage values
   THStack *hs = (THStack *) gPad->GetPrimitive("hs") ;
   RenameHistosStack( hs );
      
   #if SCALEVD>0
      //Save root file with histograms maybe for later
      #if DECONV==0
	TString rfnm = (toffset==0.)? TString(dnm) + "histos/"+ bnm + TString(".vdn")+xyz+ Form("%d",SCALEVD)  : 
                                      TString(dnm) + "histos/"+ bnm + Form("_%dns",TMath::Nint(Atrl)) + TString(".effn")+xyz+ Form("%d",SCALEVD) ;
      #else
	TString rfnm = (toffset==0.)? TString(dnm) + "histos/"+ bnm + TString("_deconv.vdn")+xyz+ Form("%d",SCALEVD)  : 
                                      TString(dnm) + "histos/"+ bnm + Form("_%dns",TMath::Nint(Atrl)) + TString("_deconv.effn")+xyz+ Form("%d",SCALEVD) ;
      #endif
      ScaleTHStack( hs , rfnm , toffset , c );
   #else
      #if DECONV==0
	TString rfnm = (toffset==0.)? TString(dnm) + "histos/"+ bnm + TString(".vd")+xyz  : 
                                      TString(dnm) + "histos/"+ bnm + Form("_%dns",TMath::Nint(Atrl)) + TString(".eff")+xyz ;
      #else
	TString rfnm = (toffset==0.)? TString(dnm) + "histos/"+ bnm + TString("_deconv.vd")+xyz  : 
                                      TString(dnm) + "histos/"+ bnm + Form("_%dns",TMath::Nint(Atrl)) + TString("_deconv.eff")+xyz ;
      #endif
      TString cmd="mv 1Dhistos.root "+ rfnm ; 
      gSystem->Exec( cmd.Data() );
      
      
   #endif
   
   //Save annealing info   
   TFile *fout = new TFile( rfnm.Data() ,"update" );
   TVectorD AnnInfo(4) ;
   AnnInfo(0) = iann  ;
   AnnInfo(1) = tann  ;
   AnnInfo(2) = TempAnn ; 
   AnnInfo(3) = Etann ;
   AnnInfo.Write("annealing");
   fout->Close();   

   //Save stacks for later (we maybe overlay them)
   if   ( toffset == 0 ) {
     if (c==0) vhvdx.push_back(hs)  ; if (c==1) vhvdy.push_back(hs)  ; if (c==2) vhvdz.push_back(hs) ;
   } else  {
     if (c==0) vhsQx.push_back(hs)  ; if (c==1) vhsQy.push_back(hs)  ; if (c==2) vhsQz.push_back(hs) ;
   }
   
   //Calculate FWHM of charge profiles
   if   ( toffset != 0 ) {

     TList *hls = (TList *) hs->GetHists()  ;
     for ( Int_t ih = 0 ; ih < hls->GetSize(); ih++ ) {
       TH1D *h1 = (TH1D * ) hls->At(ih) ;
       if ( TMath::Abs(h1->GetMinimum()) > TMath::Abs(h1->GetMaximum()) ) h1->Scale(-1);      // Gaussian fits with positive polarity don't fai
       int    bin1 = h1->FindFirstBinAbove(h1->GetMaximum()/2);
       int    bin2 = h1->FindLastBinAbove(h1->GetMaximum()/2);
       double fwhm = h1->GetBinCenter(bin2) - h1->GetBinCenter(bin1);
       double rms  = h1->GetRMS( ) ;  
       if ( c == 0 ) { vFWHMQx[ih] = fwhm ; vRMSQx[ih] = rms ; }
       if ( c == 1 ) { vFWHMQy[ih] = fwhm ; vRMSQy[ih] = rms ; }
       if ( c == 2 ) { vFWHMQz[ih] = fwhm ; vRMSQz[ih] = rms ; }
     }
     
   }
   
   //Collection time stacks
   if ( toffset==0. ) GetCollectionTimeStack( xyz ) ;
   
 
}
//------------------------------------------------------------------------
void TScan::PlotVdt( Double_t toffset , vector<THStack *> vhvd , Int_t c ) {

   /* 
      The number of plots produced is the dimension of the vector vhvd.
      Each vhvd[i] is a stack of several histograms
   
   */
   TString xyz ;
   if ( c==0 )  xyz=TString("x") ; 
   if ( c==1 )  xyz=TString("y") ; 
   if ( c==2 )  xyz=TString("z") ; 
   
   TString pthnm  = TString(dnm)+"plots/";
   TString dpdfnm = pthnm + "vd" + xyz + "_"  +TString(bnm)  ;
   TString epdfnm = pthnm + "eff" + xyz + "_"  +TString(bnm) + Form("_%dns",TMath::Nint(Atrl)) ;
   #if SCALEVD>0
      #if DECONV==0
	dpdfnm = dpdfnm + "_norm" +Form("%d",SCALEVD)+".pdf"        ; epdfnm = epdfnm + "_norm" +Form("%d",SCALEVD)+".pdf" ;
      #else
	dpdfnm = dpdfnm + "_norm" +Form("%d",SCALEVD)+"_deconv.pdf" ; epdfnm = epdfnm + "_norm" +Form("%d",SCALEVD)+"_deconv.pdf" ;
      #endif
   #else
      #if DECONV==0
	dpdfnm = dpdfnm + ".pdf"      ; epdfnm = epdfnm + ".pdf" ;
      #else
	dpdfnm = dpdfnm + ".pdf"      ; epdfnm = epdfnm + "_deconv.pdf" ;
      #endif
   #endif
   TString dpdf0  = dpdfnm+"[" , dpdff = dpdfnm+"]" , epdf0  = epdfnm+"[", epdff=epdfnm+"]"  ;
   
   Int_t NStacks = vhvd.size() ;
   if ( NStacks == 1 ) {

     THStack *hs = (THStack *) vhvd[0];
     TCanvas *c3  = new TCanvas();
     c3->SetGridx(); c3->SetGridy(); 
     hs->Draw("nostack");
     #if SCALEVD==2
       if (toffset ==0) { 
         hs->GetYaxis()->SetTitle("Drift velocity [cm/s]") ;
	 #if I2CORR==1
	   hs->GetYaxis()->SetTitle("Drift velocity [I^{2} corr, cm/s]") ;
	 #endif
       }
       if (toffset !=0) {
         hs->GetYaxis()->SetTitle("Efficiency") ;
	 #if I2CORR==1
           hs->GetYaxis()->SetTitle("Efficiency [I^{2} corr]") ;
	 #endif
	 //hs->SetMaximum(1.3); //NEW on 28th May 2013 to show annealed sensors in the same scale
       }
     #endif
     TLegend *leg = c3->BuildLegend( 0.85 , 0.11 , 0.97 , 0.71 );
     leg->SetFillColor( kWhite ) ;
     leg->Draw();
     if (toffset ==0)    c3->Print( dpdfnm.Data() );
     else                c3->Print( epdfnm.Data() );
     
     delete c3;

   } else {
     
     StackSetMaxMin( vhvd );  
     if   (toffset==0.) PlotSplit( vhvd  ,  dpdfnm ) ;
     else               PlotSplit( vhvd  ,  epdfnm ) ;
     
   }
   
}

//------------------------------------------------------------------------
void TScan::GetCollectionTimeStack( TString xyz ) {  //UNSAFE !!!
   
   //Note that in the case where we choose deconv this function is still giving the tcoll
   //from non-deconvoluted
   TString what=Form("(%s>%f && %s<%f)*(tright-tleft):%s-%f",xyz.Data(),z0-0.016,xyz.Data(),z0+thicknessz+0.016,xyz.Data(),z0) , selection ;

   //Case we want to limit the E(z) fit later on to V>Vdep
   //if (Vbhalf>0.) selection=Form("Vbias>%d",(int)Vbhalf);
   //if (Vbhalf<0.) selection=Form("Vbias<%d",(int)Vbhalf);
   if (Vbhalf>0.) selection=Form("Vbias>%d",(int)0.);
   if (Vbhalf<0.) selection=Form("Vbias<%d",(int)0.);
   plot1D((char*)what.Data(),(char*)selection.Data(),1,fnm) ;
   TString pthnm = TString(dnm)+"plots/";
   TString pdfnm = pthnm+"tcoll_" + TString(bnm)+".pdf"  ;
   gPad->Print( pdfnm.Data() ) ;
   
   TString hfnm=TString(dnm)+"histos/"+TString(bnm)+".tcoll";
   TString cmd="mv 1Dhistos.root "+ hfnm ;
   gSystem->Exec( cmd.Data() );
   
   //z0 and thickness needed by the fitting program
   TFile *fout = new TFile( hfnm ,"UPDATE" ) ;
     TVectorD z0d(2) ;
     z0d[0]=z0 ;
     z0d[1]=thicknessz;
     z0d.Write("z0d") ;
   //fout->Close(); //Deleting the file already closes it
   delete fout;
   
   //THStack *hs = (THStack *) gPad->GetPrimitive("hs") ;
   //hstcoll.push_back( hs ) ;
   //TString rfnm = pthnm+ TString(bnm)+".tcoll"  ;
   //TFile *tfile = new TFile( rfnm.Data() ,"RECREATE");
   //hs->Write("nostack");
   //tfile->Close();
   //delete tfile;
   
}
//------------------------------------------------------------------------
void TScan::ScaleTHStack( THStack *hs , TString rfnm , Double_t toffset , Int_t c ) {

   TFile *fout = new TFile( rfnm.Data() ,"recreate" );
   TH1D *h1 ;
   TList *hls = (TList *) hs->GetHists()  ;
   Double_t norm ;
   Int_t scalevd = SCALEVD ;
   Int_t WNorm   = (scalevd == 1 ) ? 1 : 0 ;
   Int_t TNorm   = (scalevd == 3 ) ? 1 : 0 ;

   if ( WNorm ) {
     /* Normalize both eff and vdrift */
     for ( Int_t ih = 0 ; ih < hls->GetSize(); ih++ ) {

       h1 = (TH1D * ) hls->At(ih) ;
       norm = h1->GetSumOfWeights();
       if (norm) h1->Scale(1.0/norm);
       h1->Write();

     }
     fout->Close();
     return;

   }
   
   /* For corrections other than normalization by weight, 
      we want efficiencies normalized to 1. 
      1st) Calculate average number of created e-h pairs. 
      2nd) Calculate actual normalization "norm" values 
   */
   if ( toffset != 0 ) { //Efficiency
     for ( Int_t ih = 0 ; ih < hls->GetSize(); ih++ ) {
     
       h1 = (TH1D * ) hls->At(ih) ;
       N0[ih] = GetN0AtVbias( h1 ) ;       
     
     }
     
     // Fit N0 vs V, fit it to Erf, get top plateau
     if (N0.Sum()<0) N0*=-1;
     N0avg = GetPlateauN0vsV( c ) ; 
     //N0avg = 35000.0;
     N0.Write("N0") ;
    
     norm = dt*1.e-9/( N0avg*AMP*Q0*ROSC ) ;   //Q(z) was not multiplied by dt, neither divided by 50 Ohm
   
   } else {              //vdrift

     //Normalization factor
     Int_t NSum   = TMath::Nint(ATRISE/dt) ;
     norm = 0.1*thicknessz / ( N0avg*AMP*Q0*NSum*ROSC ) ;  //cm/C --> vdrift [cm/s]

     //if ( TNorm ) norm = 1.0 ;
     
     //Save z0 and thickness info in vdrift file
     TVectorD z0d(2) ;
     z0d[0]=z0 ;
     z0d[1]=thicknessz;
     z0d.Write("z0d") ;
     
   } 
      
   //Save eff to disk always. Save vdrift to disk except if we ask for tcoll normalization
   for ( Int_t ih = 0 ; ih < hls->GetSize(); ih++ ) {

     h1   = (TH1D * ) hls->At(ih) ;
     h1->Scale(norm) ;

     if ( toffset!=0 ) h1->Write();                //Save normalized eff to disk
     if ( toffset == 0 && !TNorm ) h1->Write() ;   //Save vdrift to disk, TNorm is on


     //Update stack with normalized version, remove raw version
     hls->RemoveAt(ih) ;       
     hls->AddAt( h1 , ih);
   } 
   
   //Normalize vdrift by collection time
   if ( toffset ==0 && TNorm ) TCollNormalization( hls , fout ) ; 
          
   fout->Close();   
   
}   

//------------------------------------------------------------------------

Double_t TScan::GetN0AtVbias( TH1D * h1 ) {

   // Get Noavg from charge normalization, 
   Double_t SumQ   = h1->Integral( h1->FindBin(0),h1->FindBin(thicknessz) )*dt*1.e-9 ; //A*s (=C)
   Int_t    NSum   = h1->FindBin(thicknessz)-h1->FindBin(0)+1 ;
   Double_t N0V    = (SumQ/NSum)/(ROSC*AMP*Q0) ;     

   return( N0V );
    
}

//------------------------------------------------------------------------

Double_t TScan::GetPlateauN0vsV( Int_t c ) {
    
   //Build a double-erf function out of N0, then use the double Erf fit
   //Update: fill some bins before the actual N0 with 0's, to provide a smooth cte+Erf(z) fit
   Int_t offset = 10 , dim = 2*(NV+offset) ;
   TVectorD N0Erf( dim ), VErf( dim );
   
   if (N0[NV-1]> N0[0]) {
     for (Int_t i=0  ; i<2*(NV+offset) ; i++) { N0Erf[i]=N0[0]; VErf[i]=Vb[0];  }
     for (Int_t i=0  ; i<NV ; i++) {
       N0Erf[i+offset]    = N0[i]      ; VErf[i+offset]    = Vb[i]	  ;
       N0Erf[i+offset+NV] = N0[NV-1-i] ; VErf[i+offset+NV] = Vb[NV-1-i] ;
     }
   } else {
     for (Int_t i=0  ; i<2*(NV+offset) ; i++) { N0Erf[i]=N0[NV-1]; VErf[i]=Vb[NV-1];  }
     for (Int_t i=0  ; i<NV ; i++) {
       N0Erf[i+offset]    = N0[NV-1-i] ; VErf[i+offset]    = Vb[NV-1-i] ;
       N0Erf[i+offset+NV] = N0[i]      ; VErf[i+offset+NV] = Vb[i]	  ;   
     }
   }

   //Fit double-Erf to N0 vs index, then calculate the top plateau
   TH1D *hN0 = new TH1D( N0Erf );
   //hN0->Write("N0.h1");
   Double_t *pars = new Double_t [7];
   if (qfwhm[c]>0.03) pars  =  FitErf(hN0) ;
   else           FitGaussian(hN0 , pars );
   Double_t N0avg  = *(pars+6) ;
   //Double_t N0avg  = (N0.Sum()>0.)? N0.Max():N0.Min() ;
   Double_t xmean  = *pars ;
   Double_t sigmal = *(pars+4) ; 
   Int_t ileft     = hN0->FindBin(xmean+sigmal); 
   
   //NEW: We can estimate Vdep from the plot of N0 vs V, cause N0 has to reach a plateau
   //Later on, one can precisely get Vdep from CCE
   Vdep = VErf[ileft] ;
      
   TString cmd="mv Erf.pdf "+TString(dnm)+"plots/"+"N0_"+TString(bnm)+".pdf" ;
   gSystem->Exec( cmd.Data() );

   delete [] pars ;
   
   delete hN0 ; 
   
   return N0avg ; 

}

//------------------------------------------------------------------------
void TScan::TCollNormalization( TList *hls,  TFile *fout ) {
   
   fout->cd();
   
   
   /* Loop over vdrift histograms. If above depletion, then get corresponding tcoll histogram
      Apply tcoll correction (either one constant for each voltage or bin-wise).
      h1     = vdrift histogram (vs z)
      hiv    = inverse of drift velocity (vs z)
      htcoll = collection time (vs z)
   */
   for ( Int_t ih = 0 ; ih < hls->GetSize(); ih++ ) {

     if ( TMath::Abs(Vb[ih]) >= TMath::Abs(Vdep) ) {
       TH1D *h1   = (TH1D * ) hls->At(ih) ;
       TH1D *hiv  = (TH1D *) h1->Clone();
       hiv->Reset();

       for (Int_t i=1;i< h1->GetNbinsX();i++)  
	 hiv->SetBinContent( i , (h1->GetBinContent(i)!=0.)? 1./h1->GetBinContent(i) : 0. ) ;
      
       #if TCollStrat==1
         // In principle vdrift has no units, they will be forced by the constant we use to scale
	 TCollVoltageWise( ih , hiv , h1 );
       #endif
       
       #if TCollStrat==2
	 TCollBinByBin( ih  , hiv , h1 );
       #endif
      
       /*Crosscheck: calculate collection time. It has to be very similar to measured on.
         Note that vdrift came out of the normalization fulfilling Int(1/vd*dz)=tcoll [s]. 
         Since dz[=]mm, we convert to cm
        */
       for (Int_t i=1;i< h1->GetNbinsX();i++)  
	 hiv->SetBinContent( i , (h1->GetBinContent(i)!=0.)? 1./h1->GetBinContent(i) : 0. ) ;

       cout << "t_collection(V=" << Vb[ih] <<")=" << hiv->Integral( hiv->FindBin( 0. ), hiv->FindBin( thicknessz ),"width" )*0.1 <<endl ;


       h1->Write();  //CAREFUL, if you remove this, you have to change h1->Write() above

     }
     
   } 
      
}

//------------------------------------------------------------------------
void TScan::TCollVoltageWise( Int_t iVbias , TH1D* hiv , TH1D *h1 ) {

    vector<Double_t> vnorm;
    
    //Get histogram of tright-tleft for current voltage
    THStack  *hs = (THStack *) hstcoll[0] ;
    TList  *htls = (TList *)   hs->GetHists()  ;
    TH1D *htcoll = (TH1D * )   htls->At(iVbias) ;
    htcoll->GetXaxis()->SetRange(thicknessz/4.,thicknessz-0.02);
    Double_t tcollavg = htcoll->GetMean() ;
    Double_t tInt     = hiv->Integral( h1->FindBin(0.) , h1->FindBin(thicknessz/4.) ,"width" ) +
                        hiv->Integral( h1->FindBin(0.) , h1->FindBin(thicknessz-0.02) ,"width" ) ;
    Double_t norm     = tInt*0.1 / (tcollavg*1.e-9) ;
    h1->Scale( norm ) ;        
    
}

//------------------------------------------------------------------------
void TScan::TCollBinByBin( Int_t iVbias ,  TH1D *hiv , TH1D *hvd ) {

    /*
        iVbias  index of the current Vbias element
	hiv     the inverse drift velocity vs z  
    */

    //Here we use the collection time from waveforms, bin by bin, and we scale vdrift
    vector<Double_t> vnorm;
    THStack  *hs = (THStack *) hstcoll[0]       ; //stack with collection time vs z, at different Vbias
    TList  *htls =   (TList *) hs->GetHists()   ; 
    TH1D *htcoll =   (TH1D * ) htls->At(iVbias) ;
    
    //Polarity=1 for N-bulk ; -1 for P-bulk
    Int_t i0 = ( Polarity==1 ) ? htcoll->FindBin( thicknessz/4. )    : htcoll->FindBin( 0.020 ) ;
    Int_t i1 = ( Polarity==1 ) ? htcoll->FindBin( thicknessz-0.025 ) : htcoll->FindBin( 3.*thicknessz/4. )  ;


    
    for ( Int_t i=i0 ; i< i1 ; i++ ) { 
      Double_t tcollwv   = htcoll->GetBinContent(i) ;  //ns
      Double_t tcollint  = hiv->Integral( 1 , i ) * 0.1 ;   //hiv comes in no units
      Double_t norm	 = tcollint / (tcollwv*1.e-9) ;
      vnorm.push_back(norm);
      cout <<"z="<<htcoll->GetBinCenter(i)<< " vdr=" << hvd->GetBinContent(i) <<" ,tcoll="<<tcollwv<<" ns, ";
      cout <<"Int(dz/v)"<< tcollint/1.e-9<<" ns ,norm="<<norm<<endl ;
    }
    
    

#ifdef DELME
    TH1D *thiv = (TH1D *) hiv->Clone();
    thiv->Reset();
    for ( Int_t i=i0 ; i<i1 ; i++ )                  thiv->SetBinContent(i,hiv->GetBinContent(i)  ) ;
    for ( Int_t i=1  ; i<i0 ; i++ )                  thiv->SetBinContent(i,hiv->GetBinContent(i0) ) ;
    for ( Int_t i=i1 ; i<=thiv->GetNbinsX() ; i++ )  thiv->SetBinContent(i,hiv->GetBinContent(i1) ) ;

    for ( Int_t i=i0 ; i< i1 ; i++ ) { 
      Double_t tcollwv   = htcoll->GetBinContent(i) ;  //ns
      Double_t tcollint  = thiv->Integral( 1 , i ) * 0.1 ;   //hiv comes in no units
      Double_t norm	 = tcollint / (tcollwv*1.e-9) ;
      vnorm.push_back(norm);
    }
#endif

    for ( Int_t i=i0 ; i< i1 ; i++ ) {
      Double_t vdval = hvd->GetBinContent( i ) ;
      hvd->SetBinContent( i , vdval * vnorm[i-i0] ) ;
      cout << "vdrift("<<hvd->GetBinCenter(i)<<")="<<vdval * vnorm[i-i0]<<endl;
    }
    
    Double_t Mean  = TMath::Mean( vnorm.begin(), vnorm.end()) , RMS=TMath::RMS( vnorm.begin() , vnorm.end() ) ;
    WeightedMean( vnorm , Mean , RMS ) ;
    for ( Int_t i=1  ;       i< i0         ; i++ )  hvd->SetBinContent( i , hvd->GetBinContent(i) * vnorm[0] ) ;
    for ( Int_t i=i1 ; i< hiv->GetNbinsX() ; i++ )  hvd->SetBinContent( i , hvd->GetBinContent(i) * vnorm[vnorm.size()-1] ) ;
    
    vnorm.clear(); 

    
}

//------------------------------------------------------------------------
Double_t  TScan::FindErfStart( TTree *tree , TString what , TString selection ) {

       tree->Draw( what, selection,"goff" );
       Int_t nent = tree->GetSelectedRows(); Double_t *qxy  = tree->GetV1(); Double_t *yarr = tree->GetV2();
       TH1D *hqxy = new TH1D( "hqxy" ,"hqxy" , nent , yarr[0] , yarr[nent-1] );
       for ( Int_t il=0 ; il<nent ;il++ ) hqxy->SetBinContent( il , qxy[il-1]);
       Double_t *pars = FitErf( hqxy ) ;
       Double_t val = *pars ;  
       delete hqxy;
       
       return val;

}
//------------------------------------------------------------------------
void TScan::PlotCCExyz_vs_Vbias( Int_t c ) {

     //Calculate CCE
     TH1 *h1 ;
     vector<Double_t> sVb ;
   
     Int_t NStacks ; 
     if (c==0 ) NStacks = vhsQx.size() ; if (c==1 ) NStacks = vhsQy.size() ; if (c==2 ) NStacks = vhsQz.size() ; 
   
     Char_t xyz;
     if (c==0) xyz = 'x' ; if (c==1) xyz = 'y' ;  if (c==2) xyz = 'z' ; 
   
     Double_t val;
     for ( Int_t il = 0 ; il < NStacks ; il++ ) {

       THStack *hs ;
       if (c==0 ) hs = vhsQx[il] ;
       if (c==1 ) hs = vhsQy[il] ;
       if (c==2 ) hs = vhsQz[il] ;
       TList *hls = (TList *) hs->GetHists()  ;
       for ( Int_t ih = 0 ; ih < hls->GetSize(); ih++ ) {
	 h1 = (TH1D * ) hls->At(ih) ;
	 //Double_t val = 1./thickness * h1->Integral( h1->FindBin(0),h1->FindBin(thickness), "width" ); //New "width" by 16Oct2014
	 #if FIXTHICKNESS!=0
	   if ( h1->FindBin(FIXTHICKNESS/1000.) < h1->GetNbinsX() )
	      val = h1->Integral( h1->FindBin(0),h1->FindBin(FIXTHICKNESS/1000.), "width" ); //HVCMOSv3
	   else
	      val = h1->Integral( h1->FindBin(0),h1->GetNbinsX(), "width" ); //HVCMOSv3
	 #else
	   Double_t thickness ;
	   if ( c==0 ) thickness=thicknessx  ;if ( c==1 ) thickness=thicknessy ; if ( c==2 ) thickness=thicknessz ;
	   val = h1->Integral( h1->FindBin(0),h1->FindBin(thickness), "width" ); 
	 #endif
	 CCE[ih]= val  ;  //Space was already booked in the constructor
	 TString tit = h1->GetTitle() ;
	 Int_t iend = tit.Index(" V") ;
	 tit=tit(0,iend);
	 sVb.push_back( tit.Atof() ) ;
       }
     }
     
     
     //Plot CCE
     TCanvas *c4  = new TCanvas() ;
     gStyle->SetOptStat(0) ;
     gStyle->SetMarkerStyle(20) ;
     gStyle->SetMarkerSize(1) ;
     c4->SetLeftMargin(0.16); c4->SetGridx(); c4->SetGridy(); 
     TGraph *geff = new TGraph( NV , &sVb[0] , &CCE[0] ); geff->SetTitle("");
     geff->GetXaxis()->SetTitle("Bias voltage [V]");
     //geff->GetYaxis()->SetTitle(Form("CCE=#frac{1}{d}#int_{0}^{d}Q(z, 0-%d ns) #bf{dz} [a.u.]",TMath::Nint(Atrl))) ;
     
     #if FIXTHICKNESS!=0
       geff->GetYaxis()->SetTitle(Form("CC=#int_{0}^{%d} Q(%c, %d ns) #bf{dz} [a.u.]",FIXTHICKNESS,xyz,TMath::Nint(Atrl))) ;
     #else
       Double_t thickness ;
       if ( c==0 ) thickness=thicknessx  ;if ( c==1 ) thickness=thicknessy ; if ( c==2 ) thickness=thicknessz ;
       geff->GetYaxis()->SetTitle(Form("CC=#int_{0}^{%f} Q(%c, %d ns) #bf{dz} [a.u.]",thickness,xyz,TMath::Nint(Atrl))) ;
     #endif
     
     geff->GetYaxis()->SetTitleSize(0.045) ;geff->GetYaxis()->SetTitleOffset(1.4) ;
     geff->Draw( "awlp" ) ;
     #if SCALEVD>0
       TString cpdfnm = TString(dnm)+"plots/"+"CCE"+xyz+"_"+TString(bnm)+Form("_%dns",TMath::Nint(Atrl))+"_norm" +Form("%d",SCALEVD);
     #else
       TString cpdfnm = TString(dnm)+"plots/"+"CCE"+xyz+"_"+TString(bnm)+Form("_%dns",TMath::Nint(Atrl)) ;
     #endif
     #if DECONV==1
       cpdfnm = cpdfnm + "_deconv" ;
     #endif
     cpdfnm=cpdfnm+".pdf" ;
     c4->Print( cpdfnm.Data() );
     
     //Save graph to file
     #if SCALEVD>0
       TString rfnm = TString(dnm)+"histos/"+"CCE"+xyz+"_"+TString(bnm)+Form("_%dns",TMath::Nint(Atrl))+".ccen" +Form("%d",SCALEVD) ;
       #if DECONV==1
	 rfnm = TString(dnm)+"histos/"+"CCE"+xyz+"_"+TString(bnm)+Form("_%dns",TMath::Nint(Atrl))+"_deconv.ccen" +Form("%d",SCALEVD) ;
       #endif
     #else
       TString rfnm = TString(dnm)+"histos/"+"CCE"+xyz+"_"+TString(bnm)+Form("_%dns",TMath::Nint(Atrl))+".cce" ;
       #if DECONV==1
	 rfnm = TString(dnm)+"histos/"+"CCE"+xyz+"_"+TString(bnm)+Form("_%dns",TMath::Nint(Atrl))+"_deconv.cce" ;
       #endif
     #endif
     TFile *fout  = new TFile( rfnm.Data() ,"recreate" ) ;
     geff->Write() ;
     TVectorD AnnInfo(4) ;
     AnnInfo[0] = iann  ;
     AnnInfo[1] = tann  ;
     AnnInfo[2] = TempAnn ; 
     AnnInfo[3] = Etann ;
     AnnInfo.Write("annealing");
     fout->Close();
     
     delete c4;
     
     //Calculate Vdep (not needed anymore)
//      TH1D *hCCE = new TH1D( "hCCE" , "CCE vs voltage" , sVb.size() , &sVb[0] );
//      for ( Int_t i = 0 ; i < sVb.size(); i++ ) hCCE->SetBinContent( i+1 , CCE[i]  ) ;
//      Double_t *pars  =  LeftFitErf( hCCE ) ;
//      Double_t xmean  = *pars ;
//      Double_t sigmal = *(pars+2) ;
//      Vdep   = hCCE->FindBin(xmean+1.0*sigmal);   //Class Member
//      Now go and look for Capacitance file and look up the capacitance
//      Int_t    icend  = hCCE->FindBin(1.2*Vdep);
//      Cend   = (icend>sVb.size()) ? CCE[sVb.size()-1] : CCE[icend] ;  //Class Member
//      delete hCCE ;     
//      delete pars ;
     
}

//------------------------------------------------------------------------
void TScan::PlotCCEt_vs_Vbias( TTree *tree ) {

   TString what  = Form("Sum$((volt-BlineMean)*(time>%f &&time<%f)):Vbias" ,t0,t0+Atrl ) ;
   #if I2CORR==1
     what  = Form("Sum$((volt-BlineMean)*(time>%f &&time<%f))/(LPower*LPower):Vbias" ,t0,t0+Atrl ) ;
   #endif   
   
   //TString what  = Form("Sum$((volt-BlineMean)*(time>%f &&time<%f)):Vbias" ,t0,t0+50.0 ) ;
   TString pthnm = TString(dnm);
   TString pdfnm = pthnm + "plots/"+"CCE_"+TString(bnm) +".pdf" , pdf0=pdfnm+"[", pdff=pdfnm+"]";
   //TString selection = (Polarity==1) ? Form( "Vbias<%f" ,Vbmax-20) : Form( "Vbias>%f" ,Vbmax+50.) ; 
   TString selection = ""  ; 
   
   TCanvas *c4 = new TCanvas("c4","CC",600,400);
   Int_t Nval = tree->Draw(what.Data(),selection.Data(),"l") ;
   TH2D *h2D = (TH2D*) gPad->GetPrimitive("htemp");  
   h2D->GetXaxis()->SetTitle("Bias voltage [V]");

   h2D->GetYaxis()->SetTitle(Form("Q=#int_{0}^{%d} I(t) #bf{dt} [a.u.]",TMath::Nint(Atrl))) ;
   #if I2CORR==1 
     h2D->GetYaxis()->SetTitle(Form("Q=#int_{0}^{%d} I(t) #bf{dt} [I^{2}corr, a.u.]",TMath::Nint(Atrl)))  ;
   #endif
    
   c4->SetGridx(); c4->SetGridy(); 
   c4->Print(pdfnm);
   
   Double_t *Vv    = tree->GetV2() ;
   Double_t *CCEv  = tree->GetV1() ;


   //CCE[NV]. If we are in a XYZ scan, then the plot has Nx*Ny*Nz*NV points!
   if ( (ScanDir[0]+ScanDir[1]+ScanDir[2]) == 0 ) for (Int_t il=0;il<Nval ;il++) CCE[il]=CCEv[il];

   TGraph *g=new TGraph( Nval , Vv , CCEv );
   g->SetMarkerStyle(20); gStyle->SetOptTitle(0);
   g->Draw();
   g->GetXaxis()->SetTitle("Bias voltage [V]");
   g->GetYaxis()->SetTitle("CC [a.u.]");

   TString rfnm = TString(dnm)+"histos/" + bnm + TString(".cce")  ; 
   TFile *fout  = new TFile( rfnm.Data() ,"recreate" ) ;
   g->Write() ;
   fout->Close();

   delete c4 ;
     
}

//------------------------------------------------------------------------
void TScan::FitChargeProfile( Int_t coord ) {

     //File to dump the fitted histograms
     Char_t xyz;
     if (coord==0) xyz = 'x' ; if (coord==1) xyz = 'y' ;  if (coord==2) xyz = 'z' ; 
     TString hfnm = TString(dnm)+"histos/"+"FitEff"+xyz+"_"+TString(bnm)+Form("_%dns.hroot",TMath::Nint(Atrl));
     TFile *fout = new TFile( hfnm , "RECREATE");
     fout->cd();
     
     //Do the fits
     TH1D *h1 ;
     vector<Double_t> sVb ;
     Int_t NStacks ;
     if ( coord == 0 ) NStacks = vhsQx.size() ;
     if ( coord == 1 ) NStacks = vhsQy.size() ;
     if ( coord == 2 ) NStacks = vhsQz.size() ;
     vector <TH1D *> vh1;
     TString pdffit = TString(dnm)+"plots/"+"FitEff1by"+xyz+"_"+TString(bnm)+Form("_%dns",TMath::Nint(Atrl))+".pdf" ;
     TCanvas *c41  = new TCanvas("c41","Effz fits",600,400) ;
     c41->cd();
     gStyle->SetOptFit();
     Int_t Convergence;
     for ( Int_t il = 0 ; il < NStacks ; il++ ) {

       THStack *hs ;
       if ( coord == 0 ) hs = (THStack *) vhsQx[il] ;
       if ( coord == 1 ) hs = (THStack *) vhsQy[il] ;
       if ( coord == 2 ) hs = (THStack *) vhsQz[il] ;
       TList *hls = (TList *) hs->GetHists()  ;
       Double_t *pars =new Double_t [7];
       for ( Int_t ih = 0 ; ih < hls->GetSize(); ih++ ) {
	 h1 = (TH1D * ) hls->At(ih) ;
         gStyle->SetOptFit();gROOT->ForceStyle();
	 if ( TMath::Abs(h1->GetMinimum()) > TMath::Abs(h1->GetMaximum()) ) h1->Scale(-1);      // Gaussian fits with positive polarity don't fail
	 
	 #if   FITQPROF==1
	  Convergence = FitGaussian( h1 , pars );  
	 #elif FITQPROF==2
	  Convergence = FitGaussianBox( h1 , pars );  
	 #elif FITQPROF==3
	  Convergence = FitGaussianAssyBox( h1 , pars );  
	 #elif FITQPROF==4
	  Convergence = FitGaussianPol1( h1 , pars );  
	 #elif FITQPROF==5 
	  Convergence = FitGaussianAssyBox( h1 , pars );  
	 #elif FITQPROF==6 
	  Convergence = FitGaussianAssyBox( h1 , pars );  
	 #elif FITQPROF==7
	  Convergence = FitGaussianAssyBox( h1 , pars );  
	 #elif FITQPROF==8 
	  Convergence = FitGaussianAssyBox( h1 , pars );  
	 #endif
	 
	 TPaveText *doc = new TPaveText( 0.12,0.8,0.3,0.86,"NDC" );
	 doc->SetTextSize(0.029);
	 if ( Convergence==0 ) doc->AddText("Converged") ; else doc->AddText("Failed") ;
	 h1->GetListOfFunctions()->Add(doc); 
	 h1->Draw();
	 gPad->UseCurrentStyle();
	 gPad->Update();
	 h1->Write() ;
	 if (ih==0) c41->Print(pdffit+"[");
	 gPad->Print(pdffit);
	 vh1.push_back(h1) ; 
	 
	 //Check the following
         if (coord==0) vcQx[ih]     = pars[0]     ; if (coord==1)  vcQy[ih]     = pars[0]     ; if (coord==2) vcQz[ih]	   = pars[0] ;
	 if (coord==0) vFwQx[ih]    = pars[1]     ; if (coord==1)  vFwQy[ih]    = pars[1]     ; if (coord==2) vFwQz[ih]	   = pars[1] ;
         if (coord==0) vLdiffx[ih]  = pars[2]     ; if (coord==1)  vLdiffy[ih]  = pars[2]     ; if (coord==2) vLdiffz[ih]  = pars[2] ;
         if (coord==0) vGwidthx[ih] = pars[3]     ; if (coord==1)  vGwidthy[ih] = pars[3]     ; if (coord==2) vGwidthz[ih] = pars[3] ;
	 if (coord==0) vFwQx0[ih]   = pars[4]     ; if (coord==1)  vFwQy0[ih]	= pars[4]     ; if (coord==2) vFwQz0[ih]   = pars[4] ;
         if (coord==0) vcQx0[ih]    = pars[5]     ; if (coord==1)  vcQy0[ih]    = pars[5]     ; if (coord==2) vcQz0[ih]    = pars[5] ;
	 if (coord==0) vFitResx[ih] = Convergence ; if (coord==1) vFitResy[ih]  = Convergence ; if (coord==2) vFitResz[ih] = Convergence ;
       }
       delete [] pars;
       c41->Print(pdffit+"]");
   }
     
     fout->Close();
     delete c41;
     
     //Plot all fits
     TLegend *leg = new TLegend( 0.85 , 0.05 , 0.97 , 0.65 );
     TCanvas *c4  = new TCanvas("c4","Effz fits",600,400) ;
     THStack *hsfit = new THStack("hsfit","Q(z) fits");
     
     for ( UInt_t ih = 0 ; ih < vh1.size(); ih++ ) {
       hsfit->Add(vh1[ih]) ;
       TString tit = vh1[ih]->GetTitle() ;
       Int_t iend  = tit.Index(" V") ;
       tit = tit(0,iend);
       sVb.push_back( tit.Atof() ) ;
       if (coord==0) std::cout << "d(V=" <<tit.Atof()<<")="<<vFwQx[ih]<<std::endl;
       if (coord==1) std::cout << "d(V=" <<tit.Atof()<<")="<<vFwQy[ih]<<std::endl;
       if (coord==2) std::cout << "d(V=" <<tit.Atof()<<")="<<vFwQz[ih]<<std::endl;
       leg->AddEntry(vh1[ih],tit.Data(),"l");
       std::cout<<ih<<" done"<<std::endl;
     }
     
     hsfit->Draw("nostack");
     leg->SetFillColor( kWhite ) ;
     leg->Draw();
     
     hsfit->GetXaxis()->SetTitle("Bias voltage [V]");
     hsfit->GetYaxis()->SetTitle(Form("CC=#int_{0}^{%d}Q(z, 0-%d ns) #bf{dz} [a.u.]",FIXTHICKNESS,TMath::Nint(Atrl))) ;
     hsfit->GetYaxis()->SetTitleSize(0.025) ;hsfit->GetYaxis()->SetTitleOffset(1.1) ;
     #if SCALEVD>0
       TString cpdfnm = TString(dnm)+"plots/"+"FitEff_"+TString(bnm)+Form("_%dns",TMath::Nint(Atrl))+"_norm" +Form("%d",SCALEVD);
     #else
       TString cpdfnm = TString(dnm)+"plots/"+"FitEff_"+TString(bnm)+Form("_%dns",TMath::Nint(Atrl)) ;
     #endif
     #if DECONV==1
       cpdfnm = cpdfnm + "_deconv" ;
     #endif
     cpdfnm=cpdfnm+".pdf" ;
     c4->Print( cpdfnm.Data() );
          
     //delete hsfit ;
     delete c4;
}
     
//------------------------------------------------------------------------
void TScan::DumpToTree (  ) {

  //Extend this for the cases of surface scans and no Z
  
  TString fout = TString(dnm) + TString("AP_") + TString(bnm)+TString(".root");
  TFile froot( fout , "RECREATE" );

  // Create a ROOT Tree with 1 branch
  TTree *tap = new TTree("auto","eTCT measurement");
  tap->Branch("tscan", this,32000,1);
  
  /* Dump information */
  
  for ( Int_t iv = 0 ; iv<NV;iv++ ) { 

    Vbias = Vb[iv] ;

    //CCE, temperature, leakage current
    if ( XOrY()==0 ) cce   = CCE[iv];  //CC defined only for Vscans or ZVscans
    Temp     = vTemp[iv]  ;
    Itot     = vItot[iv]  ;
    
    //Width of depletion region out of fits
    FwQx  = vFwQx[iv]   ; FwQy  = vFwQy[iv]  ; FwQz  = vFwQz[iv]   ; 
    FwQx0 = vFwQx0[iv]  ; FwQy0 = vFwQy0[iv] ; FwQz0 = vFwQz0[iv]  ; 
    
    //Starting value for fits
    FWHMQx   = vFWHMQx[iv] ;
    FWHMQy   = vFWHMQy[iv] ;
    FWHMQz   = vFWHMQz[iv] ;
    
    RMSQx    = vRMSQx[iv] ;
    RMSQy    = vRMSQy[iv] ;
    RMSQz    = vRMSQz[iv] ;
    
    //Calculated error (adhoc) on depletion depth for 10 Ohm.cm and 10 um beam width
    sFWHMQx  = TMath::Exp(-0.2484-169.1*FWHMQx) ;
    sFWHMQy  = TMath::Exp(-0.2484-169.1*FWHMQy) ;
    sFWHMQz  = TMath::Exp(-0.2484-169.1*FWHMQz) ;
    
    //Position of the center of depletion witdth, calculated from fit...
    cQx    = vcQx[iv]   ; cQy    = vcQy[iv]   ; cQz    = vcQz[iv] ;

    //...and its initial value
    cQx0   = vcQx0[iv]  ; cQy0   = vcQy0[iv]  ; cQz0   = vcQz0[iv]  ;
    
    
    Ldiffx     = vLdiffx[iv]  ; Ldiffy	   = vLdiffy[iv]   ; Ldiffz     = vLdiffz[iv]  ;
    sigmax     = vGwidthx[iv] ; sigmay	   = vGwidthy[iv]  ; sigmaz     = vGwidthz[iv] ;
    FitResultx = vFitResx[iv] ; FitResulty = vFitResy[iv] ; FitResultz = vFitResz[iv] ;
    
    FitFunc  = FITQPROF        ; 
    filename = TString( bnm ); 

    //Decide what was the requested temperature
    if (Temp<-18.0)           Tset=-20.;
    if (Temp>-2.0 && Temp<2.) Tset=  0.;
    if (Temp>18.0)            Tset= 20.;
        
    //efficiency(z) profiles 
    if ( ScanDir[0] ) {
               StackToArray( vhsQx[0] , iv , Qx);
       StackAbscissaToArray( vhsQx[0] , iv , x) ;
      StackRunningQcToArray( vhsQx[0] , iv , 0) ;
               StackToArray( vhstcollx[0] , iv , tcollx);
    }
    
    if ( ScanDir[1] ) {
               StackToArray( vhsQy[0] , iv , Qy);
       StackAbscissaToArray( vhsQy[0] , iv , y) ;
      StackRunningQcToArray( vhsQy[0] , iv , 1) ;
               StackToArray( vhstcolly[0] , iv , tcolly);
    }
    
    if ( ScanDir[2]) {
               StackToArray( vhsQz[0] , iv , Qz)   ;
               StackToArray( vhvdz[0] , iv , vdrz) ;
       StackAbscissaToArray( vhsQz[0] , iv , z)   ;
      StackRunningQcToArray( vhsQz[0] , iv , 2 ) ;
      //StackToArray( hstcoll[0] , iv , tcoll) ;
               StackToArray( vhstcollz[0] , iv , tcollz);
    }
            
    tap->Fill(); 

  }

  froot.Write();
  delete tap ;
  froot.Close();
  
}
//------------------------------------------------------------------------
void TScan::StackToArray ( THStack *hs , Int_t iVb , Double_t *arr ) {

   TH1D *h1 ;
   TList *hls = (TList *) hs->GetHists()  ;
   h1 = (TH1D * ) hls->At(iVb) ;
   for ( Int_t il=1;il<=h1->GetNbinsX();il++) arr[il-1]=h1->GetBinContent(il);
   
}

//------------------------------------------------------------------------
void TScan::StackAbscissaToArray ( THStack *hs , Int_t iVb , Double_t *arr ) {

   TH1D *h1 ;
   TList *hls = (TList *) hs->GetHists()  ;
   h1 = (TH1D * ) hls->At(iVb) ;
   for ( Int_t il=1;il<=h1->GetNbinsX();il++) arr[il-1]=h1->GetBinCenter(il);
   
}

//------------------------------------------------------------------------
void TScan::StackRunningQcToArray ( THStack *hs , Int_t iVb , Int_t coord) {

   //Running charge along the z coordinate (not the time coordinate)
   TH1D *h1 ;
   TList *hls = (TList *) hs->GetHists()  ;
   h1 = (TH1D * ) hls->At(iVb) ;
   if (coord==0) for ( Int_t il=1;il<=h1->GetNbinsX();il++) rQx[il-1]=h1->Integral( 1 , il , "width" );
   if (coord==1) for ( Int_t il=1;il<=h1->GetNbinsX();il++) rQy[il-1]=h1->Integral( 1 , il , "width" );
   if (coord==2) for ( Int_t il=1;il<=h1->GetNbinsX();il++) rQz[il-1]=h1->Integral( 1 , il , "width" );
   
}

//------------------------------------------------------------------------
void TScan::Print_t0 ( TTree *tree , Double_t t0, Double_t Vbmax , Double_t Cend, Double_t dt, TString dnm , TString bnm ) {
     
     TCanvas *ccont = new TCanvas("ccont","ccont",600,400);
     TString selection = Form( "Vbias==%f && (time>%f && time<%f) && (event%%10==1)", Vbmax, t0-2.0, t0+10.0 ) ;

     #if DECONV == 0
       tree->Draw("volt-BlineMean:time",selection,"l") ;
     #else
       TString what = Form("%5.2f*%7.5f*(volt[Iteration$+1]-volt[Iteration$])/%5.3f+volt-BlineMean:time" , ROSC , Cend , dt);    
       tree->Draw(what.Data(),selection.Data(),"l") ; 
     #endif
     TH1D *htemp = (TH1D*) gPad->GetPrimitive("htemp");  
     if (htemp==0) return;
     Double_t Vmin = htemp->GetYaxis()->GetXmin() ;
     Double_t Vmax = htemp->GetYaxis()->GetXmax() ;
     TLine *t0l=new TLine( t0, Vmin , t0, Vmax );
     t0l->Draw();
     
     TString vfnm = TString(dnm)+"plots/" + TString("vt_") +bnm + TString(".pdf") ;
     ccont->Print( vfnm.Data() ) ;
     
     delete ccont ; 
     delete t0l ;
}
//------------------------------------------------------------------------
 
void TScan::WeightedMean( vector<Double_t> tl , Double_t &Mean , Double_t &RMS) {

     TH1D *h1=new TH1D("h1","h1",100,Mean-RMS,Mean+RMS);
     //TH1D *h1=new TH1D("h1","h1",100,Mean-6*RMS,Mean+6*RMS);
     h1->Sumw2();
     for ( UInt_t ii=0 ; ii < tl.size() ; ii++ ) h1->Fill( tl[ii] ) ;
     //If the underflow and overflow of the histogram is negligible, then proceed with Mean,RMS evaluation
     if ( (h1->GetBinContent(0) + h1->GetBinContent( h1->GetNbinsX()+1)) != tl.size() ) {
       Mean = h1->GetMean() ;
       RMS  = h1->GetRMS()  ;
     }
          
     delete h1 ; 

}
//------------------------------------------------------------------------
 
void TScan::WeightedMean( vector<Double_t> tl , Double_t min , Double_t max , Double_t &Mean , Double_t &RMS) {

     TH1D *h1=new TH1D("h1","h1",100,min,max);
     h1->Sumw2();
     for ( UInt_t ii=0 ; ii < tl.size() ; ii++ ) h1->Fill( tl[ii] ) ;
     //If the underflow and overflow of the histogram is negligible, then proceed with Mean,RMS evaluation
     if ( (h1->GetBinContent(0) + h1->GetBinContent( h1->GetNbinsX()+1)) != tl.size() ) {
       Mean = h1->GetMean() ;
       RMS  = h1->GetRMS()  ;
     }
          
     delete h1 ; 

}
 
/*
//------------------------------------------------------------------------
void TScan::WeightedMean( vector<Double_t> tl , Double_t &Mean , Double_t &RMS) {

     //TH1D *h1=new TH1D("h1","h1",100,Mean-0.5*RMS,Mean+0.5*RMS);
     Double_t tmin = TMath::MinElement(tl.size(),&tl[0]) , 
     	      tmax = TMath::MaxElement(tl.size(),&tl[0]) ;
	      
     TH1D *h1=new TH1D("h1","h1",100,tmin,tmax);
     h1->Sumw2();
     for ( UInt_t ii=0 ; ii < tl.size() ; ii++ ) h1->Fill( tl[ii] , 1.0 ) ;
     //If the underflow and overflow of the histogram is negligible, then proceed with Mean,RMS evaluation
     if ( (h1->GetBinContent(0) + h1->GetBinContent( h1->GetNbinsX()+1)) != tl.size() ) {
       Mean = h1->GetMean() ;
       RMS  = h1->GetRMS()  ;
     }
          
     delete h1 ; 

}
 */
 
//------------------------------------------------------------------------
void TScan::RenameHistosStack( THStack *hs ) {

   TH1D *h1 ;
   TList *hls = (TList *) hs->GetHists()  ;
   TString tit ; 
   for ( Int_t ih = 0 ; ih < hls->GetSize(); ih++ ) {
     h1 = (TH1D * ) hls->At(ih) ;
     tit = h1->GetTitle() ;
     Int_t iV = tit.Index("V=");
     Int_t icom = tit.Index(", ");
     tit=tit(iV+2,icom-(iV+2));
     h1->SetTitle( tit ) ;
   }
      
}

//------------------------------------------------------------------------
void TScan::StackSetMaxMin( vector<THStack *> hsv ) {
    
     Double_t Max=-9.e12, Min=9.e12 ;
     for ( Int_t ih = 0 ; ih<hsv.size() ; ih++ ) {
       Double_t Mx0 = hsv[ih]->GetMaximum("nostack"), Mn0 = hsv[ih]->GetMinimum("nostack"); 
       if ( Mx0 > Max ) Max = Mx0 ;
       if ( Mn0 < Min ) Min = Mn0 ;
     }

     for ( Int_t ih = 0 ; ih<hsv.size() ; ih++ ) {
       hsv[ih]->SetMinimum(Min) ; hsv[ih]->SetMaximum(Max) ;
     }
    
}

//------------------------------------------------------------------------
void TScan::PlotSplit( vector<THStack *> hsv , TString fnm ) {

     TLegend *leg ; 
     TCanvas *c3 = new TCanvas();
     hsv[0]->Draw("nostack")     ; leg = c3->BuildLegend( 0.85 , 0.05 , 0.97 , 0.65 )  ; 
     leg->SetFillColor( kWhite ) ; leg->Draw();
     TString open=fnm+"[" ;
     c3->Print( open.Data() )    ; c3->Print( fnm.Data() ) ; 
     
     hsv[1]->Draw("nostack")     ;  leg = c3->BuildLegend( 0.85 , 0.05 , 0.97 , 0.65 )  ; 
     leg->SetFillColor( kWhite ) ; leg->Draw();
     TString close=fnm+"]" ;
     c3->Print( fnm )            ; c3->Print( close.Data() ) ; 

}
//------------------------------------------------------------------------
string TScan::GetName(  ) {
   return name ;
}

//------------------------------------------------------------------------
Int_t TScan::XOrY(  ) {
   
   Int_t val = ( (ScanDir[0]+ScanDir[1])>=1 ) ? 1 : 0 ;
   return val ;
   
}
//------------------------------------------------------------------------
void TScan::VsBias( TTree *tree , TMeas *em  ) {
    
    TProfile *hIp , *hTp ;
    if (NV>1) {

      //Some magnitudes vs Vbias
      //We need variable bin sizes. Therefore add one lower-edge
      Double_t *arr = new Double_t [NV+1];
      for ( Int_t il=0 ; il<NV;il++) arr[il]=Vb[il];

      //Sort in increasing order, just in case the bias vector comes in descending order
      Int_t index[NV] ;
      TMath::Sort( NV , arr , index , kFALSE );
      Double_t *acpy = new Double_t [NV+1];
      for ( Int_t il=0 ; il<NV;il++) acpy[il] = arr[index[il]]  ;
      for ( Int_t il=0 ; il<NV;il++) arr[il]  = acpy[il] ;


      arr[NV]=arr[NV-1]+(arr[NV-1]-arr[NV-2]);
      hIp = new TProfile("hIp","Itot vs Vbias",NV,arr,"s") ;
      hTp = new TProfile("hTp","Temp vs Vbias",NV,arr,"s") ;
      delete [] arr ; delete [] acpy ;
      
    } else {
      
      hIp = new TProfile("hIp","Itot vs Vbias",1,Vb[0]-0.1*TMath::Abs(Vb[0]),Vb[0]+0.1*TMath::Abs(Vb[0]),"s") ;
      hTp = new TProfile("hTp","Temp vs Vbias",1,Vb[0]-0.1*TMath::Abs(Vb[0]),Vb[0]+0.1*TMath::Abs(Vb[0]),"s") ;
    }
    
    Long_t nentries = tree->GetEntries() ; 
    for ( Int_t ii=0 ; ii < nentries ; ii++ ) {

      tree->GetEntry(ii);
      hIp->Fill( em->Vbias, em->Itot );
      hTp->Fill( em->Vbias, em->Temp );

    }

    for ( Int_t ii=0 ; ii < NV ; ii++ ) {
      vItot[ii] = hIp->GetBinContent(ii+1); 
      vTemp[ii] = hTp->GetBinContent(ii+1); 
    }
    
    delete hIp; delete hTp; 
    
}

//------------------------------------------------------------------------
Int_t TScan::FitGaussian( TH1D* h1 , Double_t *pars ) {

        TCanvas *c1=new TCanvas("c1","Gaussian fits",600,400);
        Int_t Convergence = h1->Fit("gaus");
        c1->Update();
        c1->Print("Erf.pdf");
        TFormula *fit = (TFormula *) h1->GetFunction("gaus");

        Double_t Mean=fit->GetParameter(1);
        Double_t Sigma=fit->GetParameter(2);
        std::cout<<"Sigma="<<Sigma<<std::endl;
        Double_t Thickness = 4.*Sigma ;


	pars[0]=Mean ;
	pars[1]=Thickness ;
	pars[2]=0. ;
	pars[3]=Sigma ;
	pars[4]=0.;
	pars[5]=0.;
	pars[6]=0.;
   
        std::cout << "Thickness="<<Thickness<<std::endl;
	
        delete c1;
	
	return Convergence ;
	
}

//------------------------------------------------------------------------
Int_t TScan::FitGaussianPol1( TH1D* h1 , Double_t *pars ) {

        //Guessing starting parameters
	h1->Fit("gaus");
	TFormula *fit = (TFormula *) h1->GetFunction("gaus");
	Double_t norm=fit->GetParameter(0);
	Double_t mu=fit->GetParameter(1);
	Double_t sigma=fit->GetParameter(2);

        TCanvas *c1=new TCanvas("c1","Gaussian+pol1",600,400);
	TF1 *gausp1 = new TF1( "gausp1" , "gaus(0)+pol1(3)" , h1->GetXaxis()->GetXmin() , h1->GetXaxis()->GetXmax() );
        gausp1->SetParameters(norm,mu,sigma,0.,0.6);
        Int_t Convergence = h1->Fit( "gausp1","","",mu-3.*sigma,mu+5.*sigma );
	gStyle->SetOptFit();
        c1->Update();
        c1->Print("Erf.pdf");
        Double_t Mean=gausp1->GetParameter(1);
        Double_t Sigma=gausp1->GetParameter(2);
        Double_t Thickness = 4.*Sigma ;


	pars[0] = Mean;
	pars[1] = Thickness;
	pars[2] = 0.;
	pars[3] = Sigma ;
	pars[4] = sigma;
	pars[5] = mu ;
	pars[6] = 0.;
   
        std::cout << "Thickness="<<Thickness<<std::endl;
	
        delete c1;
	
	return Convergence ;
	
}

//------------------------------------------------------------------------
Int_t TScan::FitGaussianBox( TH1D* h1 , Double_t *pars ) {

        //Guessing starting parameters
	h1->Fit("gaus");
	TFormula *fit = (TFormula *) h1->GetFunction("gaus");
	Double_t norm=fit->GetParameter(0);
	Double_t mu=fit->GetParameter(1);
	Double_t sigma=fit->GetParameter(2);

        TCanvas *c1=new TCanvas("c1","Gaussian conv. Box",600,400);
	TF1 *gausbox = new TF1( "gausbox" , GaussBox , h1->GetXaxis()->GetXmin() , h1->GetXaxis()->GetXmax(),4 );
        gausbox->SetParameters(mu,norm,2.*sigma,0.01);
        gausbox->SetParNames("Box center","Norm","Box Length","sigma");
	gausbox->FixParameter(3,0.01); //Fix beam width
        Int_t Convergence = h1->Fit( "gausbox","WW","",0.,0.09 );
	gStyle->SetOptFit();
        c1->Update();
        c1->Print("Erf.pdf");
        Double_t Mean      = gausbox->GetParameter(0);
        Double_t BoxLength = gausbox->GetParameter(2);
	Double_t Gwidth    = TMath::Abs(gausbox->GetParameter(3));


	pars[0]=Mean;
	pars[1]=BoxLength;
	pars[2]=0.;
	pars[3]=Gwidth;
	pars[4]=2*sigma ;
	pars[5]=mu;
	pars[6]=0.;
   
        std::cout << "Thickness="<<BoxLength<<std::endl;
	
        delete c1;
	
	return Convergence ;
	
}

//------------------------------------------------------------------------
Int_t TScan::FitGaussianAssyBox( TH1D* h1 , Double_t *pars ) {

//Suggestions: pack error of 2 um. We are fitting the RMS of mean values
//but each of those has a 2 um error!
//https://root.cern.ch/phpBB3/viewtopic.php?t=9320

        //Guessing starting parameters: gaussian
	Double_t mean = h1->GetMean() , rms = h1->GetRMS();
	h1->Fit("gaus","WWR","",mean-rms,mean+rms);
	gPad->Print("gauss.pdf");
	TFormula *fit = (TFormula *) h1->GetFunction("gaus");
	Double_t mu=fit->GetParameter(1);
	Double_t sigma=fit->GetParameter(2);
        Double_t eLdiff, bexpo ;
	
	#if FITQPROF==5 || FITQPROF==6 || FITQPROF==8
          //Guessing starting parmeters: exponential
	  h1->Draw();h1->Fit("expo","WW","R,sames",mu+0.66*sigma,h1->GetXaxis()->GetXmax());
	  cout<<"[]="<<mu+0.2*sigma<<" "<<h1->GetXaxis()->GetXmax()<<endl;
	  gPad->Print("expo.pdf");
	  fit = (TFormula *) h1->GetFunction("expo");
	  bexpo = fit->GetParameter(1);
	  eLdiff = 1./TMath::Abs(bexpo) ;
	#elif FITQPROF==3 || FITQPROF==7
          //Guessing starting parmeters: exponential
	  h1->Draw();h1->Fit("pol1","WW","R,sames",mu+0.66*sigma,h1->GetXaxis()->GetXmax());
	  cout<<"[]="<<mu+0.2*sigma<<" "<<h1->GetXaxis()->GetXmax()<<endl;
	  gPad->Print("expo.pdf");
	  fit = (TFormula *) h1->GetFunction("pol1");
	  bexpo = fit->GetParameter(1);	
	  eLdiff =  -fit->GetParameter(0)/fit->GetParameter(1) - mu;
	#endif
	
	//Estimate depleted thickness
	int bin1 = h1->FindFirstBinAbove(h1->GetMaximum()/2);
	int bin2 = h1->FindLastBinAbove(h1->GetMaximum()/2);
	double fwhm = h1->GetBinCenter(bin2) - h1->GetBinCenter(bin1);
	//Double_t edepl = fwhm - 0.005 - eLdiff ; //Removing beam width (left), and diffusion length (right)
	Double_t edepl = TMath::Sqrt(fwhm*fwhm - 0.005*0.005)  ; //Removing beam width (left), and diffusion length (right)
	if (edepl<0) edepl=0.;

        TString tit = h1->GetTitle() ;
        Int_t iend  = tit.Index(" V") ;
        tit = tit(0,iend);	
	Double_t eVb = tit.Atof();
//	Double_t edepl = 0.3e-3*TMath::Sqrt(OHMCM*TMath::Abs(eVb)) ;
			   
        //TMinuit *gMinuit = new TMinuit(5);  //initialize TMinuit with a maximum of 5 params

        TCanvas *c1=new TCanvas("c1","Gaussian conv. Assymetric Box",600,400);
	Double_t norm=h1->Integral(h1->FindBin(-0.02),h1->FindBin( 0.09),"width");
	TF1 *gausabox ;
	#if FITQPROF==3
	  gausabox = new TF1( "gausabox" , GaussAssyBox , h1->GetXaxis()->GetXmin() , h1->GetXaxis()->GetXmax() , 5 );
	  gausabox->SetParameters(mu,norm,edepl,0.01,0.06);
          gausabox->SetParNames("#mu_{box}","Norm","#sigma_{box}","#sigma_{laser}","d_{diff}");
	#elif FITQPROF==5
	  gausabox = new TF1( "gausabox" , GaussExpBox , -0.02 , 0.09, 5 );
	  gausabox->SetParameters(mu,norm,edepl,0.01,bexpo);
	  cout<<mu<<" "<<norm<<" "<<edepl<<" "<<bexpo<<endl;
          gausabox->SetParNames("#mu_{box}","Norm","#sigma_{box}","#sigma_{laser}","b_{diff}");
	#elif FITQPROF==6
	  Double_t qbl   = h1->Integral(1,h1->FindBin(-0.01))/(1+h1->FindBin(-0.01)) ;
	  Double_t slope = (h1->GetBinContent(h1->GetNbinsX())-h1->GetBinContent(1))/(h1->GetXaxis()->GetXmax() - h1->GetXaxis()->GetXmin()) ;
	  //gausabox = new TF1( "gausabox" , GaussExpBoxBline , h1->GetXaxis()->GetXmin() , h1->GetXaxis()->GetXmax(), 7 );
	  //gausabox->SetParameters(mu,norm,edepl,0.01,bexpo,qbl ,slope );
	  gausabox = new TF1( "gausabox" , GaussExpBoxBline , h1->GetXaxis()->GetXmin() , h1->GetXaxis()->GetXmax(), 6 );
	  gausabox->SetParameters(mu,norm,edepl,0.01,bexpo,qbl ,slope );
          gausabox->SetParNames("#mu_{box}","Norm","#sigma_{box}","#sigma_{laser}","b_{diff}", "Bline");
	#elif FITQPROF==7
	  Double_t qbl   = h1->Integral(1,h1->FindBin(-0.01))/(1+h1->FindBin(-0.01)) ;
	  gausabox = new TF1( "gausabox" , GaussAssyBox , h1->GetXaxis()->GetXmin() , h1->GetXaxis()->GetXmax() , 6 );
	  //gausabox->SetParameters(mu,norm,edepl,0.01,0.06);
	  gausabox->SetParameters(0.5*edepl,norm,edepl,0.01,eLdiff,qbl);
          gausabox->SetParNames("#mu_{box}","Norm","#sigma_{box}","#sigma_{laser}","d_{diff}", "Bline");
	#elif FITQPROF==8
	  Double_t qbl   = h1->Integral(1,h1->FindBin(-0.01))/(1+h1->FindBin(-0.01)) ;
	  Double_t slope = (h1->GetBinContent(h1->GetNbinsX())-h1->GetBinContent(1))/(h1->GetXaxis()->GetXmax() - h1->GetXaxis()->GetXmin()) ;
	  gausabox = new TF1( "gausabox" , GaussExpBoxBline , h1->GetXaxis()->GetXmin() , h1->GetXaxis()->GetXmax(), 6 );
	  gausabox->SetParameters(0.5*edepl,norm,edepl,0.01,bexpo,qbl );
          gausabox->SetParNames("#mu_{box}","Norm","#sigma_{box}","#sigma_{laser}","b_{diff}", "Bline");
        #endif

	//Double_t sbn=0.3*edepl; if (sbn<0) sbn=0.; //Seems ok for 2e16
	//Double_t sbn = 0.1e-3*TMath::Sqrt(OHMCM*TMath::Abs(eVb)) ;
	gausabox->FixParameter(2,edepl);   //Fix box width
	gausabox->FixParameter(3,0.01);   //Fix beam width
	#if FITQPROF!=3 && FITQPROF!=7
	  gausabox->FixParameter(4,bexpo);   //Fix exponential decay constant
	#endif

	//Double_t mean=h1->GetMean() , rms=h1->GetRMS(); 
	//gMinuit->Command("SET STR 0");
	//gMinuit->Command("SIM");
	Double_t RangeMax =  h1->GetXaxis()->GetXmax();
	if ( RangeMax > 0.09 )  RangeMax=0.09;
	cout<<"Start..."<<h1->GetName()<<endl;
	std::cout<<"Init parameters: fwhm="<<fwhm<<" ["<<"center="<<0.5*edepl<<" width="<<edepl<<" Lwidth="<<1./TMath::Abs(bexpo)<<std::endl;
	
	Int_t Convergence ;
	
        #if FITQPROF==6 || FITQPROF==7 || FITQPROF==8 || FITQPROF==3 
	  Convergence = h1->Fit( "gausabox","WW","",h1->GetXaxis()->GetXmin(), h1->GetXaxis()->GetXmax()); //Set all weights =1
        #else  
	  Convergence = h1->Fit( "gausabox","WW","",-0.02,0.09 ); //Set all weights =1
	#endif

        gausabox->ReleaseParameter(2);
	gausabox->SetParLimits(2,0.,0.5); //Limit box dimension
        gausabox->ReleaseParameter(3);
        gausabox->ReleaseParameter(4);
	#if FITQPROF==6 || FITQPROF==7 || FITQPROF==8 || FITQPROF==3 
	  Convergence = h1->Fit( "gausabox","WW","",h1->GetXaxis()->GetXmin(), h1->GetXaxis()->GetXmax()); //Set all weights =1
	#else  
	  Convergence = h1->Fit( "gausabox","WW","",-0.02,0.09 ); //Set all weights =1
        #endif
	
	cout<<"End of "<<h1->GetName()<<endl;
	gStyle->SetOptFit();
        c1->Update();
        c1->Print("Erf.pdf");
        Double_t Mean       = gausabox->GetParameter(0);
        Double_t BoxLength  = TMath::Abs(gausabox->GetParameter(2));
	Double_t Gwidth     = TMath::Abs(gausabox->GetParameter(3));
	Double_t FallLength ;
	#if FITQPROF==3 || FITQPROF==7
          FallLength = gausabox->GetParameter(4);
	#elif FITQPROF==5 || FITQPROF==6  || FITQPROF==8
	  Double_t b=gausabox->GetParameter(4);
          FallLength = TMath::Abs(1./b);
        #endif
        Double_t Thickness = BoxLength ;


	pars[0]=Mean       ;
	pars[1]=Thickness  ;
	pars[2]=FallLength ;
	pars[3]=Gwidth     ;
	pars[4]=edepl ;   //Starting value for wV
	pars[5]=mu ;     //Starting value for cdepl
	pars[6]=0. ;
   
        //std::cout << "Thickness="<<Thickness<<std::endl;
	
        delete c1;
	
	return Convergence ;
	
}

//------------------------------------------------------------------------

Double_t GaussBox( Double_t *x, Double_t *par ) {

//par[0]=box center, par[1]=Normalization, par[2]=box width  , par[3]=Gaussian width  

   Double_t bw        = 0.001  ;  //Bin width
   Double_t BoxLength = par[2] ;
   Double_t Norm      = par[1] ;
   Double_t halfBW    = 0.5*BoxLength ;
   Double_t sigma     = par[3] ;

   Double_t Xmin      = par[0]-1.5*BoxLength ;
   Double_t Xmax      = par[0]+1.5*BoxLength ;
   Xmin      = -0.2 ;
   Xmax      = 0.2 ;
   
   Int_t  Nbins = 1+TMath::Nint((Xmax-Xmin)/bw) ; 
         
   //Gaussian is always infinitely wide, so in principle there is always overlap between the box and gaussian
   Double_t sum = 0. ;
   for ( Int_t j=1;j<=Nbins;j++  ) {
     Double_t xit = Xmin+(j-1)*bw ;
     Double_t box = ( (par[0]-halfBW) <= xit && xit<= (par[0]+halfBW) ) ? 1./BoxLength : 0. ; //Normalizing the box so area =1
     sum+= Norm * TMath::Gaus(xit,x[0],sigma,kTRUE) * box * bw ;
   }

   
   return sum ;

}

//------------------------------------------------------------------------

Double_t GaussAssyBox( Double_t *x, Double_t *par ) {

//par[0]=box center, par[1]=Normalization, par[2]=box width  , par[3]=Gaussian width  , par[4]=Diffussion Region

   Double_t bw         = 0.002  ;  //Bin width
   Double_t BoxCenter  = par[0] ;
   Double_t BoxLength  = par[2] ;
   Double_t DiffLength = par[4] ;
   Double_t Bline      = par[5] ;
   Double_t Xmin       = -0.2 ;
   Double_t Xmax       = 0.2 ;
   Double_t Norm       = par[1] ;
   Double_t sigma      = par[3] ;
   
   Int_t  Nbins = 1+TMath::Nint((Xmax-Xmin)/bw) ; 
         
   Double_t sum = 0. , box=0. ;
   for ( Int_t j=1;j<=Nbins;j++  ) {
     Double_t xit = Xmin+(j-1)*bw ;
     Double_t x1 = BoxCenter-0.5*BoxLength , x2 = BoxCenter+0.5*BoxLength , x3=x2+DiffLength;
     Double_t norm = 1./(1.0*BoxLength+0.5*DiffLength*1.0) ;
     #if   FITQPROF==3
       box = norm*( ( xit<x1) *0. +
                           ( x1 <=xit && xit<=x2 ) *  1.0 + 
                           ( xit>x2  && xit<=x3 ) * (1./DiffLength * (BoxCenter+0.5*BoxLength+DiffLength-xit)) +
			   ( xit>x3 ) * 0. );
     #elif FITQPROF==7	
       box = norm*(  
                          0.5*(TMath::TanH(1000.*(xit-x1)) - TMath::TanH(1000.*(xit-x2)))
			+ 0.5*(TMath::TanH(1000.*(xit-x2)) - TMath::TanH(1000.*(xit-x3))) *
			        (1./DiffLength * (BoxCenter+0.5*BoxLength+DiffLength-xit))
			 ) ;
     #endif	   
     sum+=Norm * TMath::Gaus(xit,x[0],sigma,kTRUE) * box * bw ;
   }

   sum+=Bline  ;

   
   return sum ;

}
//------------------------------------------------------------------------

Double_t GaussExpBox( Double_t *x, Double_t *par ) {

//par[0]=box center, par[1]=Normalization, par[2]=box width  , par[3]=Gaussian width  , par[4]=Diffussion Region

   Double_t bw         = 0.001  ;  //Bin width
   Double_t BoxCenter  = par[0] ;
   Double_t Norm       = par[1] ;
   Double_t BoxLength  = par[2] ;
   Double_t sigma      = par[3] ;
   Double_t bexpo      = par[4] ;
   Double_t DiffLength = -1./bexpo;
   
   
   //Double_t Xmin  = BoxCenter-1.5*BoxLength ;
   //Double_t Xmax  = BoxCenter+0.5*BoxLength+DiffLength ;
   Double_t Xmin  = -0.2 ;
   Double_t Xmax  = 0.2 ;
   Int_t  Nbins   = 1+TMath::Nint((Xmax-Xmin)/bw) ; 
   
         
   Double_t x1    = BoxCenter-0.5*BoxLength , x2 = BoxCenter+0.5*BoxLength , sum = 0. ;
   Double_t aexpo = -bexpo*x2 ;
   Double_t Sexp  = TMath::Exp(aexpo)/bexpo*( TMath::Exp(bexpo*Xmax) - TMath::Exp(bexpo*x2) );
   Double_t iarea = 1./(1.0*BoxLength+Sexp) ;
   
   for ( Int_t j=1;j<=Nbins;j++  ) {

     Double_t xit = Xmin+(j-1)*bw ;

     Double_t box = iarea*( ( xit < x1) * 0. + 
     			    ( x1 <= xit && xit < x2 ) *  1.0 + 
                            ( xit>=x2 ) * ( TMath::Exp(aexpo+bexpo*xit) ) );
     
     sum+=Norm * TMath::Gaus(xit,x[0],sigma,kTRUE) * box * bw ;
   }

   //std::cout << "x1="<<x1<<" x2=" << x2 << " Sexp="<<Sexp  <<std::endl;
   return sum ;

}
//------------------------------------------------------------------------

Double_t GaussExpBoxBline( Double_t *x, Double_t *par ) {

//par[0]=box center, par[1]=Normalization, par[2]=box width  , par[3]=Gaussian width  , par[4]=Diffussion Region

   Double_t bw         = 0.001  ;  //Bin width
   Double_t BoxCenter  = par[0] ;
   Double_t Norm       = par[1] ;
   Double_t BoxLength  = par[2] ;
   Double_t sigma      = par[3] ;
   Double_t bexpo      = par[4] ;
   Double_t Bline      = par[5] ;
   //Double_t slope      = par[6] ;
   
   
   Double_t Xmin  = -0.2 ;
   Double_t Xmax  = 0.2 ;
   Int_t  Nbins   = 1+TMath::Nint((Xmax-Xmin)/bw) ; 
   
         
   Double_t x1    = BoxCenter-0.5*BoxLength , x2 = BoxCenter+0.5*BoxLength , sum = 0. ;
   Double_t aexpo = -bexpo*x2 ;
   Double_t Sexp  = TMath::Exp(aexpo)/bexpo*( TMath::Exp(bexpo*Xmax) - TMath::Exp(bexpo*x2) );
   Double_t iarea = 1./(1.0*BoxLength+Sexp) ;
   
   Double_t box ;
   for ( Int_t j=1;j<=Nbins;j++  ) {

     Double_t xit = Xmin+(j-1)*bw ;

     #if   FITQPROF==6
       box = iarea*( ( xit < x1) * 0. + 
     			    ( x1 <= xit && xit < x2 ) *  1.0 + 
                            ( xit>=x2 ) * ( TMath::Exp(aexpo+bexpo*xit) ) );
     #elif FITQPROF==8	
       box = iarea*(  
                          0.5*(TMath::TanH(1000.*(xit-x1)) - TMath::TanH(1000.*(xit-x2)))*1.0
			+ 0.5*(TMath::TanH(1000.*(xit-x2)) - TMath::TanH(1000.*(xit-Xmax))) *
			        ( TMath::Exp(aexpo+bexpo*xit))
			 ) ;
     #endif	   
     
     sum+= Norm * TMath::Gaus(xit,x[0],sigma,kTRUE) * box * bw ;
   }
   //sum+=Bline + slope*x[0] ;
   sum+=Bline  ;
   
   
   //std::cout << "x1="<<x1<<" x2=" << x2 << " Sexp="<<Sexp  <<std::endl;
   return sum ;

}
