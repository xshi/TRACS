#include "vector"
#include <algorithm>
#include <iostream>
#include <fstream>
#include <sstream> //ostringstream
#include <string>  //strcpy

void PlotCarrierFile( char *fnm ){

  string what , line ;
  stringstream myStream ;
  
  vector<Double_t> q,y,z ;   //Using TRACS reference frame
  ifstream in( fnm );
  if (!in || in.bad()) {
    cout << "Error opening " << fnm << endl ;
  }

  Double_t qval , yval , zval, tval ;
  while (!in.eof() && !in.bad()) {
    getline(in,line) ; 
    myStream.str(""); myStream.clear() ; myStream << line ;
    myStream >> what >> qval >> yval >> zval >> tval;  
    q.push_back(qval/1.602e-19);
    y.push_back(yval);
    z.push_back(zval);
    //std::cout<<yval<<std::endl;
  }
  Double_t ymin=*std::min_element(y.begin(),y.end()) ;
  Double_t ymax=*std::max_element(y.begin(),y.end()) ;
  Double_t zmin=*std::min_element(z.begin(),z.end()) ;
  Double_t zmax=*std::max_element(z.begin(),z.end()) ;
//   Double_t Ax=0.04;
//   Double_t Ay=0.04;
//   xmin=xmin-0.5*Ax ; xmax=xmax+0.5*xmax;
//   ymin=ymin-0.5*Ay ; ymax=ymax+0.5*ymax;
  //Int_t Ny=TMath::Nint((ymax-ymin)/Ay);
  //Int_t Nz=TMath::Nint((zmax-zmin)/Az);

  TH2D *hq = new TH2D( "hq" , "carriers" , 40 , ymin,ymax, 40,zmin,zmax) ;
  for (Int_t ii=0 ; ii<y.size() ;ii++ ) hq->SetBinContent( hq->FindBin(y[ii],z[ii]),q[ii] );
  
  
  hq->Draw("colz");
  gStyle->SetOptStat(0);
  gPad->SetGrid(1);
  
  TFile *f = new TFile("SimHisto.root","RECREATE");
  hq->Write();
  f->Close();

}
