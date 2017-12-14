// root -l 
// .L CompareTRACS2Meas.C
// CompareTRACS2Meas("2017-07-24_16-14-41_dmil_PIN_8622_W5_E3_4.txt.root","NOirrad_dt100ps_4pF_tNOtrappingns_dz5um_dy5dV5V_0nns_bottom_0_rc.hetct.root")

void CompareTRACS2Meas( char *fmeas , char *fsim , Double_t Norm=1.){

  TFile *fm = new TFile( fmeas ) ;
  TFile *fs = new TFile( fsim )  ;
  TTree *trm = (TTree *) fm->Get("edge") ;
  TTree *trs = (TTree *) fs->Get("edge") ;
  trm->SetLineColor(1) ; trm->SetLineWidth(2) ;
  trs->SetLineColor(2) ; trs->SetLineWidth(2) ;

  trs->Draw("event", "" ,"goff");
  Int_t Nsim = trs->GetSelectedRows() ;
  
  TCanvas *c1 = new TCanvas();
  c1->Print("CompareTRACS2Meas.pdf[");
  c1->SetGrid();
  gStyle->SetOptDate();
  for ( Int_t i=0 ; i<Nsim ; i++ ) {
   
    TString tev=Form("event==%d",i);
    trs->Draw("Vbias", tev ,"goff");
    Double_t Volt= TMath::Mean(trs->GetSelectedRows(),trs->GetV1());
    tev=Form( "TMath::Abs(Vbias)==%f",TMath::Abs(Volt) ) ;
    
    if (Norm==1.) trm->Draw( "(volt-BlineMean)/Q50:time" , tev , "l"  ) ; 
    else         trm->Draw( "(volt-BlineMean):time" , tev , "l"      ) ; 
    

    TGraph *gm = new TGraph(trm->GetSelectedRows(),trm->GetV2(), trm->GetV1());
    tev=Form( "%d" , TMath::Nint(Volt) ) ;
    gm->SetNameTitle( tev.Data() , tev.Data() );
    TAxis *axis = gm->GetXaxis();
    axis->SetLimits(-5.,35.);    
    gm->Draw("awl");
    gm->GetXaxis()->SetTitle("Time [ns]"); gm->GetYaxis()->SetTitle("Signal [a.u.]"); 
                
    tev=Form( "TMath::Abs(Vbias)==%f",TMath::Abs(Volt) ) ;

    if (Norm==1) trs->Draw( "(volt-BlineMean)/Q50:time" , tev , "lsame" ) ; 
    else {
      TString what = Form( "%f*(volt-BlineMean):time",Norm );
      trs->Draw( what , tev , "lsame"     ) ; 
    }
    
    
    TGraph *gs = new TGraph(trs->GetSelectedRows(),trs->GetV2(), trs->GetV1());
    tev=Form( "%d" , TMath::Nint(Volt) ) ;
    gs->SetNameTitle( tev.Data() , tev.Data() );
    
    //c1->BuildLegend() ;
    
    c1->Print("CompareTRACS2Meas.pdf");
    
  }

  c1->Print("CompareTRACS2Meas.pdf]");
  
}
