#ifndef TSCANHEADER_H
#define TSCANHEADER_H

#include <TObject.h>

#define EXE 2             //1 for root .so, 2 for executable

#define SCALEVD 0         //0 do not scale neither efficiency nor vdrift
                          //1 Scale efficiency and vdrift by SumOfWeights
			  //2 Scale efficiency and vdrift to N0 (Ner of eh pairs)
			  //3 Scale vdrift by collection time
			  
#define FIXTHICKNESS 0   //0             Use calculated detector thickness in PlotEff()
                          //Non-zero [um] Use this value as thickness Int_{0}{val} Q(z)dz in PlotEff
			  //		  THIS VALUE HAS TO BE INTEGER (preprocessor constraint), so given in microns
			  //Examples: 90 n-irrad HVCMOSv3, 140 p-irradiated HVCMOS
			  
#define DECONV	0         //1 RC-deconvolute the current pulses	
			  //NOTE THAT SINCE APRIL 2013 WE CAN NOW CREATE THE ROOT FILE WITH THE DECONVOLUTED DATA INSIDE
			  //THEREFORE THIS OPTION HERE IS NOT NEEDED ANYMORE	  

#define FITQPROF 0        //0 Do not fit
			  //1 Fit to Gaussian
			  //2 Fit to Gaussian convoluted with a box
			  //3 Fit to Gaussian convoluted with assymetric box
			  //4 Fit to Gaussian + pol1(=background)
			  //5 Fit to Gaussian convoluted with box ending in a falling exponential
			  //6 Fit to Gaussian convoluted with box ending in a falling exponential + constant for background
			  //7 same as 3 but piecewise with tanh()+ constant for background
			  //8 same as 6 but piecewise with tanh()

#define I2CORR	0	  //0 no correction of waveforms
			  //1 waveforms are weighted by wvf=wvf*1/Ilaser^2 (TPA)
			  
#define WEIGHTINGFIELD 0  //0 no weighting field needed
			  //1 HVCMOS Ew calculated using TRACS
			  //2 HVCMOS Ew calculated with TCAD

#define BULKCOORD 1       //

#define OHMCM 10

#if SCALEVD==3
  #define TCollStrat 1      //1 Applies correction voltage by voltage
                            //2 Applies correction bin by bin			  
#endif

#define vZ0 0.0033        //Used by tcollection correction. Below this value, tright-tleft seems to be too high
 
//#define ZLIMIT "z<7.15"   //Useful to fit double-Erf of sensors with a difficult structure. COMMENT IT IF NOT NEEDED!!!

// GRAPHIC OUTPUT OPTIONS ................................................

#define NPX 1             //Number of plots in X (2D plots)
#define NPY 1             //Number of plots in Y (2D plots)
#define PORTRAIT  0       //=1 portrait, otherwise is landscape

#define HIGHQPDF 1        //1 Creates a big pdf out of the 2D distros, otherwise say 0


#define ATRL 15          //0 in case we want the program to calculate the time window for plotting
			  //non-zero in case we want a fixed time zoom in plots

#define TCOLL 20          //0 in case we want the collection time 2D plot to have no bonds
			  //non-zero in case we want a fixed time in CollTime(x,y)


#define ALLVTZ 0         //0 plot fewer V(t) for particular zi in same plot
                         //1 plot all V(t)_zi in same plot

#define AUTO 1                   //To be used in some cases where AutoPilot cannot help
#if AUTO==0                      //Then, we advance by hand the parameters that AutoPilot has to guess
  #define Z0         7.028
  #define T0         5.017
  #define THICKNESS  0.07845
  #define ZL	     7.04
  #define ZR	     7.105
  #define POLARITY   1
#endif

//----------------END OF USER EDITABLE ----------------------

#define Q0      1.602e-19   //C
#define AMP     100         //C 41 db in power~112 in voltage (http://www.muzique.com/schem/gain.htm)
//#define AMP     8.62         //C eTCT Miteq after mulfunctioning
//#define AMP     1.0        //C eTCT Miteq after mulfunctioning
#define ATRISE  0.4         //ns (this is the time over which we integrate the vdrift [ formula (4.4), Nicolaś thesis ]
#define ROSC    50.0  //Ohm
                     
#if EXE==1 || EXE==2
#include <cstdarg>
#include <string.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <libgen.h>

#include "TFile.h"
#include "TTree.h"
#include "THStack.h"
#include "TCanvas.h"
#include "TClassTable.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TH1.h"
#include "TH2D.h"
#include "TF1.h"
#include "TGraph.h"
#include "TAttLine.h"
#include "TTreeFormula.h"
#include "TEventList.h"
#include "TEntryListArray.h"
#include "TLegend.h"
#include "TString.h"
#include "TMath.h"
#include "TLatex.h"
#include "TLine.h"
#include "TVectorD.h"
#include "TProfile.h"
#include "TROOT.h"
#include "TPaletteAxis.h"
#include "TPaveText.h"

#include "TMeas.h"
#include "TMeasHeader.h"
#include "TWaveform.h"
#include "TScan.h"
#include "TMinuit.h"

#endif



class TScan : public TObject {   //TScan derives from TObject, so we can use some of its methods
  
  public:
  
     TScan(  )  ;
     TScan( char *fnm , Double_t Cend )  ;
     ~TScan() ;
     
     int      Polarity  ;     //Polarity of file
     Double_t BLMean    ;     //File average of BlineMean
     Double_t Cend      ;     //!Estimated end capacitance (20% above Vdep)
     Double_t Temp      ;     //Average temperature for current bias
     Double_t Tset      ;     //Intended temperature (Ex.:-19.76C=-20C)

     Double_t Vbmax     ;     //Maximum bias voltage
     Double_t Vbmin     ;     //Minimum bias voltage
     Double_t Vbhalf    ;     //!This is the middle voltage in the Vbias vector
     Double_t Vmax      ;     //Maximum signal in the scope, for all measurements
     Double_t Vmin      ;     //Minimum signal in the scope, for all measurements
     Double_t Vdep      ;     //!Estimated depletion voltage (from CCE measurements)
     Double_t Vbias     ;     //Bias voltages


     Double_t t0        ;     //Origin of times
     Double_t tlavg     ;     //!<tleft> at Vmax
     Double_t travg     ;     //!<tright> at Vmax
     Double_t Atrl      ;     //<tr-tl>

     Double_t qfwhm[3]	;    //!Rough: FWHM of Q50 at Vbmax

     Int_t    NV        ;     //Number of Bias voltages
     Int_t    Nx        ;     //Number of steps in X direction 
     Int_t    Ny        ;     //Number of steps in Y direction 
     Int_t    Nz        ;     //Number of steps in Z direction 

     Double_t x0        ;     //Fitted position where sensor starts 
     Double_t y0        ;     //Fitted position where sensor starts
     Double_t z0        ;     //Fitted position where sensor starts

     Double_t xl        ;     //!Position where Q50@Vmax starts
     Double_t yl        ;     //!Position where Q50@Vmax starts
     Double_t zl        ;     //!Position where Q50@Vmax starts

     Double_t xr        ;     //!Position where Q50@Vmax ends
     Double_t yr        ;     //!Position where Q50@Vmax ends
     Double_t zr        ;     //!Position where Q50@Vmax ends

     Double_t x0file    ;     //Origin of z array in the file
     Double_t y0file    ;     //Origin of z array in the file
     Double_t z0file    ;     //Origin of z array in the file

     Double_t thicknessz;     //Thickness
     Double_t thicknessx;     //Thickness
     Double_t thicknessy;     //Thickness

     Double_t dx        ;     //Step in x
     Double_t dy        ;     //Step in y
     Double_t dz        ;     //Step in z

     Int_t    ScanDir[4];     //(X,Y,Z,V)=(1,1,0,0) is a XY scan ; (0,0,1,1)=ZVscan and so on
     			      //!Note: Nicola's eTCT always has Z scan involved
			      //!TCT+: Always Z involved except for Normal TCT: Vscan or XY(V) scan

     Double_t dt        ;     //Step in times (scope)
     Double_t N0avg     ;     //!Average number of e-h pairs produced

     Double_t iann      ;     //annealing step [allowing steps like "3.5"]
     Double_t tann      ;     //Specific annealing time at TempAnn
     Double_t Etann     ;     //Accumulated annealing time at 60C for this sample
     Double_t TempAnn   ;     //Annealing temperature
     Double_t Fluence   ;     //Fluence
     

     string name        ;     //!Sensor name
     char   *fnm        ;     //!Filename with path
     char   *bnm        ;     //!Filename without path
     char   *dnm        ;     //!Path
     TString filename   ;     //Filename

     Double_t cce       ; //Charge collection Efficiency vs Vbias
     Double_t Itot      ; //Leakage current vs Vbias
     
     Double_t FwQx      ; //Fitted box width=depleted region for current bias w(V)
     Double_t FwQy      ; //Fitted box width=depleted region for current bias w(V)
     Double_t FwQz      ; //Fitted box width=depleted region for current bias w(V)

     Double_t FwQx0     ; //Starting value given for depl width fit
     Double_t FwQy0     ; //Starting value given for depl width fit
     Double_t FwQz0     ; //Starting value given for depl width fit

     Double_t FWHMQx	; //FWHM of charge collection region, for current bias (x direction) [was deplFWHM]
     Double_t FWHMQy	; //FWHM of charge collection region, for current bias (y direction)
     Double_t FWHMQz	; //FWHM of charge collection region, for current bias (z direction)

     Double_t RMSQx    ; //FWHM of charge collection region, for current bias (x direction) [was deplFWHM]
     Double_t RMSQy    ; //FWHM of charge collection region, for current bias (y direction)
     Double_t RMSQz    ; //FWHM of charge collection region, for current bias (z direction)

     Double_t sFWHMQx   ; //FWHM_eTCT-FWHM_theory parametrization y=TMath::Exp(-0.2484-169*x) [sigma_laser=10 um] 
     Double_t sFWHMQy   ; //FWHM_eTCT-FWHM_theory parametrization y=TMath::Exp(-0.2484-169*x) [sigma_laser=10 um] 
     Double_t sFWHMQz   ; //FWHM_eTCT-FWHM_theory parametrization y=TMath::Exp(-0.2484-169*x) [sigma_laser=10 um] 

     Double_t cQx        ; //Fitted center (mu) of depletion region for current bias
     Double_t cQy        ; //Fitted center (mu) of depletion region for current bias
     Double_t cQz        ; //Fitted center (mu) of depletion region for current bias

     Double_t cQx0       ; //Starting value for cdepl
     Double_t cQy0       ; //Starting value for cdepl
     Double_t cQz0       ; //Starting value for cdepl

     Double_t Ldiffx     ; //Diffusion length (exponential decay)
     Double_t Ldiffy     ; //Diffusion length (exponential decay)
     Double_t Ldiffz     ; //Diffusion length (exponential decay)

     Double_t sigmax     ; //Gaussian width of the beam
     Double_t sigmay     ; //Gaussian width of the beam
     Double_t sigmaz     ; //Gaussian width of the beam
     
     Int_t    FitResultx ; //0 if converged (See: "Access to the fit status" in TH1)
     Int_t    FitResulty ; //0 if converged (See: "Access to the fit status" in TH1)
     Int_t    FitResultz ; //0 if converged (See: "Access to the fit status" in TH1)

     Int_t   FitFunc     ;  //Fitting function used (=FITEFF)


     
     Double_t *x        ; //[Nx]
     Double_t *Qx       ; //[Nx]
     Double_t *tcollx   ; //[Nx]
     Double_t *rQx      ; //[Nx] Running charge along the x coordinate (not the time coordinate)
     
     Double_t *y        ; //[Ny]
     Double_t *Qy       ; //[Ny]
     Double_t *tcolly   ; //[Ny]
     Double_t *rQy      ; //[Ny] Running charge along the y coordinate (not the time coordinate)

     Double_t *z        ; //[Nz]
     Double_t *Qz       ; //[Nz] (was called "eff" before. I prefer now this Q(z))
     Double_t *tcollz   ; //[Nz]
     Double_t *rQz      ; //[Nz] Running charge along the z coordinate (not the time coordinate)
     Double_t *vdrz     ; //[Nz] (was called "vdr" before)
     Double_t *tcoll    ; //[Nz]

     Double_t  *Vb      ; //![NV] 
     Double_t  *vItot   ; //![NV] Average leakage current vs bias voltage
     Double_t  *vTemp   ; //![NV] Average temperature vs bias voltage

     Double_t  *vFwQx   ;    //![NV] Fitted depleted thickness (depends on voltage). thickness should be wV[NV-1] (renamed from dofV)
     Double_t  *vFwQy   ;    //![NV] 
     Double_t  *vFwQz   ;    //![NV] 

     Double_t  *vFwQx0  ;    //![NV] Starting value for dofV fit (dofV0)
     Double_t  *vFwQy0  ;    //![NV] 
     Double_t  *vFwQz0  ;    //![NV] 
     
     Double_t  *vFWHMQx ;    //![NV] FWHM of depleted thickness (depends on voltage)
     Double_t  *vFWHMQy ;    //![NV] 
     Double_t  *vFWHMQz ;    //![NV] 

     Double_t  *vRMSQx  ;    //![NV] Differs from FWHM for assymetric distributions
     Double_t  *vRMSQy  ;    //![NV] Note that FWHM=2.355*RMS (if gaussian)
     Double_t  *vRMSQz  ;    //![NV] 

     Double_t  *vcQx    ;    //![NV] Center of depleted thickness (depends on voltage). 
     Double_t  *vcQy    ;    //![NV] 
     Double_t  *vcQz    ;    //![NV] 

     Double_t  *vcQx0	;    //![NV] Starting value of the fit for cdepV. 
     Double_t  *vcQy0	;    //![NV] 
     Double_t  *vcQz0	;    //![NV] 
     
     Double_t  *vLdiffx  ;    //![NV] Center of depleted thickness (depends on voltage). 
     Double_t  *vLdiffy  ;    //![NV] 
     Double_t  *vLdiffz  ;    //![NV] 

     Double_t  *vGwidthx ;    //![NV] Fitted gaussian width of beam (depends on voltage). 
     Double_t  *vGwidthy ;    //![NV] 
     Double_t  *vGwidthz ;    //![NV] 
     
     Int_t     *vFitResx ;    //![NV] Convergence of the fit (for each voltage) [0=OK]
     Int_t     *vFitResy ;    //![NV] 
     Int_t     *vFitResz ;    //![NV] 

     TVectorD   N0 	 ;     //![NV] Vector with number of electron hole pairs produced by laser 
                              //(calculated from Q(z)) 
     
     vector<THStack*> vhsV_zi  ;    //!NV   thstacks of V(t) at fixed z	, different Vbias
     vector<THStack*> vhsV_Vbi  ;   //!NV   thstacks of V(t) at fixed Vbias, different z
     
     vector<TLegend*> vlhsV_zi ;    //!NV   Legends
     vector<TLegend*> vlhsV_Vbi ;   //!NV   Legends
     
     vector<TH2D*>    vhvtx  ;    //!NV   hVtx[2]->Draw("colz") (Get them from the stack)
     vector<TH2D*>    vhvty  ;    //!NV   hVty[2]->Draw("colz") (Get them from the stack)
     vector<TH2D*>    vhvtz  ;    //!NV   hVtz[2]->Draw("colz") (Get them from the stack)
     
     vector<TH2D*>    vhQ2D      ; //!NV   2D Q(x,y)
     vector<TH2D*>    vhtcoll2D  ; //!NV   2D collection time map   
     vector<TH2D*>    vhvd2D     ; //!NV   2D vd(x,y)
     vector<TH2D*>    vhvdEw2D   ; //!NV   2D vd(x,y)*WeightingField
     vector<Double_t> CCE        ; //!NV
     
     
     vector<THStack*>  vhvdx ;    //!Using a vector for the case we divide vdrift plots in 2 sets
     vector<THStack*>  vhvdy ;    //!Using a vector for the case we divide vdrift plots in 2 sets
     vector<THStack*>  vhvdz ;    //!Using a vector for the case we divide vdrift plots in 2 sets
     
     vector<THStack*>  vhsQx ; //! Stack of Q(x): 1 or 2 elements, in case we plot all Eff's in 1 or 2 plots
     vector<THStack*>  vhsQy ; //! Stack of Q(y): 1 or 2 elements, in case we plot all Eff's in 1 or 2 plots
     vector<THStack*>  vhsQz ; //! Stack of Q(z): 1 or 2 elements, in case we plot all Eff's in 1 or 2 plots
          
     vector<THStack*>  vhstcollx ; //! Stack of tcoll(x): 1 or 2 elements, in case we plot all Eff's in 1 or 2 plots
     vector<THStack*>  vhstcolly ; //! Stack of tcoll(y): 1 or 2 elements, in case we plot all Eff's in 1 or 2 plots
     vector<THStack*>  vhstcollz ; //! Stack of tcoll(z): 1 or 2 elements, in case we plot all Eff's in 1 or 2 plots
          
     vector<THStack*> hstcoll ; //!Collection time for each Z, at each voltage
     


     //----> Methods 
     
     void     RoughTRTL( TString vselection , TTree *tree , TWaveform *wv , Double_t *trtl ) ;   //!
     void     RoughTRTL( TString vselection , TTree *tree , Double_t *trtl , Double_t Cend ) ;   //!
     void      VoltVsTime_ci ( Int_t coord ) ;   //!
     void      VoltVsTime ( )   ;   //!
     void	 TimeCoordVolt2D ( Int_t coord ) ;  //!
     void     TimeVbiasVolt2D() ;   //!
     void     Maps2D( Int_t iwhat ) ; //!
     void     Q2D_IntegrationTime( ) ; //!
     void                Vdt( Double_t toffset , char *cond , Int_t coord ) ;            //!
     void            PlotVdt( Double_t toffset , vector<THStack *> hvd , Int_t coord ) ;              //!
     void            FitChargeProfile( Int_t c );                                        //! Fits Eff(z) for each bias
     void            PlotCCExyz_vs_Vbias( Int_t c ) ;			                 //!	  
     void            PlotCCEt_vs_Vbias( TTree *tree ) ;                                  //! Normal TCT, no Z coordinate
     void            PlotVar( Int_t iwhat , Int_t c ) ;                                  //! Normal TCT, no Z coordinate
     void      GetCollectionTimeStack( TString xyz ) ;		                         //!	  
     void      ScaleTHStack( THStack *hs , TString rfnm , Double_t toffset , Int_t c ) ; //! 
     void      DumpToTree( );                                              //!
     void      StackToArray( THStack *hs , Int_t iVb , Double_t *arr );    //!
     void      StackAbscissaToArray( THStack *hs , Int_t iVb , Double_t *arr );   //!
     void      StackRunningQcToArray( THStack *hs , Int_t iVb , Int_t coord );   //!
     string    GetName()              ;   //! Return sensor name
     void      VsBias( TTree *tree , TMeas *em );
     void      CalcWhereSignal( Int_t coord ) ; //!
     Int_t     FitGaussian( TH1D *h1 , Double_t *pars ) ;          //!
     Int_t     FitGaussianBox( TH1D *h1 , Double_t *pars ) ;       //!
     Int_t     FitGaussianAssyBox( TH1D *h1 , Double_t *pars ) ;   //!
     Int_t     FitGaussianPol1( TH1D* h1 , Double_t *pars ) ;         //!
     Int_t     XYZ;            //!1=x, 2=y, 12=xy, 3=z, 13=xz, 23=yz
     Int_t     XOrY();
     
            
private:

     Double_t  GetN0AtVbias( TH1D * h1 ) ;                                 //!
     Double_t  GetPlateauN0vsV( Int_t c ) ;                                       //!
     void      TCollNormalization( TList *hls,  TFile *fout ) ;            //!
     void      TCollVoltageWise( Int_t iVbias, TH1D *hiv , TH1D *h1 ) ;    //!

     void      TCollBinByBin( Int_t iVbias , TH1D *hiv, TH1D *hv  ) ;	   //!
     
     Double_t  FindErfStart( TTree *tree, TString what , TString selection ) ;

protected:

     void      Print_t0 ( TTree *tree , Double_t t0, Double_t Vbmax , Double_t Cend, Double_t dt, TString dnm , TString bnm ) ; //!
     void      WeightedMean( vector<Double_t> tl , Double_t &Mean , Double_t &RMS ) ;   //!
     void      WeightedMean( vector<Double_t> tl , Double_t min , Double_t max, Double_t &Mean , Double_t &RMS) ; //!
     void      RenameHistosStack( THStack *hs ) ;                                                 //!
     void      StackSetMaxMin( vector<THStack *> hsv  ) ;                                         //!
     void      PlotSplit( vector<THStack *> hsv , TString fnm ) ;                                 //!
     void      vd_times_Ew() ; //!   

ClassDef(TScan,1) ; //ROOT RTTI

};

#endif
