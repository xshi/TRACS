//Creates a new TF evaluated each 50 ps (as the measurements are)

root -l /home/jcalvopi/TRACS_Concurrency/myApp/Centered_100ps_TransferFunction_Cividec_06052014.root 

Int_t Nbins = shtf->GetNbinsX() ;
Double_t x1 = shtf->GetXaxis()->GetXmin() , x2 = shtf->GetXaxis()->GetXmax() ;
TH1D *htf = new TH1D( "htf" , "C2 TF 50 ps" , 2*Nbins , x1 , x2 ) ;

for ( Int_t i=1 ; i<=2*Nbins ; i++ ) {
  Double_t xval = htf->GetBinCenter( i ) ;
  htf->SetBinContent( i , shtf->Interpolate( xval ) ) ;
}

TFile *fout = new TFile( "/home/jcalvopi/TRACS_Concurrency/myApp/Centered_50ps_TransferFunction_Cividec_14112017.root" , "RECREATE" );
htf->Write();
fout->Close();
