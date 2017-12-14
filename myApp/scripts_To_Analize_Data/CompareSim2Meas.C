// Example:
// root -l 
//  .L CompareSim2Meas.C
//  CompareSim2Meas(file1,file2)
//

void CompareSim2Meas( TString fmeas , TString fsim ){
  
  TFile *file0 = TFile::Open( fmeas );
  TFile *file1 = TFile::Open( fsim );
  
  TTree *trm=(TTree*) file0->Get("edge") ;
  TTree *trs=(TTree*) file1->Get("edge") ;
  
  trm->SetLineColor(1)    ; trs->SetLineColor(2)    ;  
  trm->SetLineWidth(2)    ; trs->SetLineWidth(2)    ;  
  trs->SetMarkerStyle(20) ; trs->SetMarkerColor(2)  ; 
  
  //Compare all in one go
  TString pdfnm0=fsim+".pdf[" , pdfnm=fsim+".pdf" , pdfnmf=fsim+".pdf]";
  trm->Draw("volt:time","","l") ;
  trs->Draw("volt:time","","psame") ;
  gPad->Print( pdfnm0 ); gPad->SetGrid(1);
  //htemp->GetXaxis()->SetTitle("Time [ns]") ;   htemp->GetYaxis()->SetTitle("Signal [A]") ; 
  gPad->Print( pdfnm );
  
  
  for (Int_t i=0 ; i<4;i++) {
    TString how=TString( Form("event==%d" , i)  ) ;
    trm->Draw("volt:time", how ,"l") ;
    trs->Draw("volt:time", how ,"lpsame") ;
    gPad->SetGrid(1);
    //htemp->GetXaxis()->SetTitle("Time [ns]") ;   htemp->GetYaxis()->SetTitle("Signal [A]") ; 
    gPad->Print( pdfnm );
  }   
  gPad->Print( pdfnmf );
  
  
}
