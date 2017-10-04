// Hay que grabar (tcollx,tcolly,tcollz) =(tright-tleft,...,..)
//

/*

   To avoid problems with rounding in VoltVsTime_ci, calculate what is the step width in the scan and use multiples of it!!!

   Analyzes eTCT measurements without external intervention.

   Use (t0,z0)   
   
   List of plots wished:
   ---------------------
   
   	1) Volt:time at different depths   (for each Vbias)
	2) Volt:time at different voltages (at selected depths)
	3) time:z:volt-BlineMean           (at all voltages)
	4) Drift velocity profiles         (at all voltages)
	5) Efficicenct                     (at all voltages)
   
   Input:   Any number of files
   ------	
   
   Output:  All the histograms generated (?) in root format and pdf format
   -------

   Note: I prefer not to add this analysis to the original root-file, cause
   that one can be O(100 MB) big, and we dont want to move this around everytime
   we update the analysis.

   TODO: 
   
   1) plotvd uses plot1D to plot "Sum$((volt-BlineMean)*(...):time". However we should better plot:
   "Sum$((volt-Average(BlineMean))*(...)):z" for each voltage. This is overall important for low voltages. 
   Note that we can not use Sum$(BlineMean) cause BlineMean is not an array.
   So we finally need to loop for each voltage, calculate vd and fill a THStack.
   As of today, 30th Nov 2012, I added the option in plotvd to use the file average of BLineMean to calculate
   the drift velocity. It seems to pack better efficiency curves Q=Q(z) at high Vbias (see BLINEMEAN==1 or 2 in plotvd.C)
   
   2) When I use more than 1 argument, then argv[2] (2nd filename) gets corrupted. Does it have to do with va_list 
   somewhere?
   
   2) WE NEED TO CHANGE Atrl for each voltage when calculating efficiency?
   3) In cases where there is no depletion... do we need to include diffusion??

Investigate:
   dlopen error: TMeas_cpp.so: cannot open shared object file: No such file or directory
   Load Error: Failed to load Dynamic link library /home/mfg/hpk/ETCT-Analyse/cpp/./TWaveform_cpp.so


*/

#include "TScan.h"
ClassImp(TScan);

void separator() ;

//------------------------------------------------------------------------
void separator() {

     cout<<endl;
     cout<<"---------------------------------------------------------------------------"<<endl;
     cout<<endl;
}

//------------------------------------------------------------------------

int main( int argc , char *argv[] ) {

   /* Fonts and sizes */
   gStyle->SetOptStat(0) ;
   gStyle->SetTextFont(22)        ; gStyle->SetTextSize(0.045)  ; 
   gStyle->SetTitleFont(22,"x")   ; gStyle->SetTitleFont(22,"y"); gStyle->SetTitleFont(22,"z");
   gStyle->SetTitleSize(0.045,"x"); gStyle->SetTitleSize(0.045,"y");
   gStyle->SetLabelFont(22,"x")   ; gStyle->SetLabelSize(0.045,"x");
   gStyle->SetLabelFont(22,"y")   ; gStyle->SetLabelSize(0.045,"y");
   gStyle->SetLabelFont(22,"z")   ; gStyle->SetLabelSize(0.045,"z");
   gStyle->SetHistLineWidth(2)    ; gStyle->SetFuncWidth(2); 
   
   TH1::AddDirectory(kFALSE);

   
   Int_t istep = 1 ;
   #if DECONV == 1
     istep=2 ;
     if (argc!=3) {
       cout << "You have requested to deconvolute the RC," <<endl;
       cout << "but you forgot to provide Cend in the argument list:"<<endl;
       cout << "AutoPilot_eTCT fnm Cend"<<endl;
       return 0;
     }
   #endif
   
   char empty[]="";
   Double_t Cend = 0. ;
   for (Int_t im = 1 ; im <argc ; im=im+istep ) {

    
     #if DECONV == 1   //Provisional, only checking if Ideconv makes sense. Later on, look up in file
       Cend = atof( argv[2] ) ;  //Input in pF, then changing to nF, so RC=ns
       Cend = Cend / 1.e3;
     #endif
     TScan *Scan = new TScan( argv[im] , Cend ) ;
     
     //Encode the movement into a variable
     Scan->XYZ=0;  //No movement

     /*
         				        x				   xy					xz
         				        y				   yz
         				        z
        If 0<XYZ<10 we have  1D scan
	If XYZ>10 we have  2D scan
	If XYZ>100 we have 3D scan
     */
     if      ( Scan->ScanDir[0] ) { Scan->XYZ = 1 ; if (Scan->ScanDir[1]) Scan->XYZ=12 ; if (Scan->ScanDir[2]) Scan->XYZ=13 ; if (Scan->ScanDir[1] && Scan->ScanDir[2]) Scan->XYZ=123 ;} 
     else if ( Scan->ScanDir[1] ) { Scan->XYZ = 2 ; if (Scan->ScanDir[2]) Scan->XYZ=23 ; }
     else if ( Scan->ScanDir[2] ) { Scan->XYZ = 3 ; } 
     
     
     /* Normal TCT without any coordinate motion */
     if ( Scan->XYZ == 0 ) { 
        
	Scan->VoltVsTime( )	  ; separator() ;
        Scan->TimeVbiasVolt2D( )  ; separator() ;
        Scan->DumpToTree( ) ;	  
      
        delete Scan ;

	continue;

     } else {

        //1D and 2D scans
	for ( Int_t c=0 ; c<3 ; c++ ) {  //HERE:Vdt fixed (see TString Coord). I think fails because of thstacks having 2 instances in a 2D scan
          if ( Scan->ScanDir[c] ) {
             
	     //1D plots: V(t) at different values of the moving coord (centered in the other coord)
	     Scan->VoltVsTime_ci( c )   ; separator() ;

	     //2D plots: V(t,c) (centered in the other coord)
	     Scan->TimeCoordVolt2D( c ) ; separator() ;
             
             //Q(coord) (centered in the other coord)
	     Scan->Vdt( Scan->Atrl  , "" , c )  ; 
	     if (c==0) Scan->PlotVdt( Scan->Atrl , Scan->vhsQx , c ) ;   
	     if (c==1) Scan->PlotVdt( Scan->Atrl , Scan->vhsQy , c ) ;   
	     if (c==2) Scan->PlotVdt( Scan->Atrl , Scan->vhsQz , c ) ;   
	     separator() ; 
	     
	     //Drift velocity=f(coord)
	     Scan->Vdt( 0.  , "" , c)  ; 
	     if (c==0) Scan->PlotVdt( 0. , Scan->vhvdx , c ) ;    
	     if (c==1) Scan->PlotVdt( 0. , Scan->vhvdy , c ) ;    
	     if (c==2) Scan->PlotVdt( 0. , Scan->vhvdz , c ) ;   
	     separator() ;   
	     
	     //CCE
	     Scan->PlotCCExyz_vs_Vbias( c ) ;  separator() ;    

	     //CCE
	     Scan->PlotVar( 1, c ) ;  separator() ;    //0 (empty) , 1=collection time

	     #if FITQPROF!=0
	       Scan->FitChargeProfile( c )      ;  separator() ;    
	     #endif
	     
          }
	  
	}

	//Only for 2D scans
	if ( Scan->XYZ > 10 ) {
	  Scan->Q2D_IntegrationTime();
	  //Scan->Maps2D( 0 ) ;  // 0=Charge (Deprecated: see Q2D_IntegrationTime)
	  Scan->Maps2D( 1 ) ;    // 1=Collection Time
	  Scan->Maps2D( 2 ) ;    // 2=Drift velocity
	}
	
     }
     
     Scan->DumpToTree( ) ;    
     
     delete Scan ;

   }
   
   return 0;     

}
