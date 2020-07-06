#include <iostream>
#include <fstream>
#include <sstream>
#include <string>


int plot_data(std::string filename)
{
    gStyle->SetOptStat(000);

    TFile *f = new TFile(filename.c_str());
    
    TTree* tree = dynamic_cast<TTree*>(gDirectory->Get("edge"));

    TCanvas *MyC = new TCanvas("MyC", "",800,600);    
    tree->Draw("volt:time", "event==120");

    MyC->Print("HPK_w13_SE2_event120.gif");

    
    return 0;
}
