#include <iostream>


int plot_sim_data()
{
    gStyle->SetOptStat(000);

    std::stringstream hist_title_meas;
    std::stringstream hist_title_sim_spline;

    hist_title_meas << "i_total_meas" ;
    hist_title_sim_spline << "i_total_sim_spline" ;

    TFile *f1 = new TFile("./out/fit_current_latest.root");

    TH1D* htmp_meas = dynamic_cast<TH1D*>(gDirectory->Get(hist_title_meas.str().c_str()));
    TH1D* htmp_sim_spline = dynamic_cast<TH1D*>(gDirectory->Get(hist_title_sim_spline.str().c_str()));


    int nbin = htmp_meas->GetXaxis()->GetNbins();
    //std::cout << "nbin = " << nbin << std::endl;

    int nbin1 = htmp_sim_spline->GetXaxis()->GetNbins();
    //std::cout << "nbin1 = " << nbin1 << std::endl;

    TH1D * htmp_meas_conv = (TH1D*)htmp_meas->Clone();
    TH1D * htmp_sim_spline_conv = (TH1D*)htmp_sim_spline->Clone();

    TCanvas *MyC_tmp = new TCanvas("MyC_tmp","",1000,600);    
    gPad->SetBottomMargin(0.17);
    gPad->SetLeftMargin(0.16);

    MyC_tmp->SetTicks(1,1);
  
    htmp_meas_conv->Scale(-1.0);
    htmp_sim_spline_conv->Scale(-1.0);

    htmp_meas_conv->SetMarkerColor(2);
    htmp_sim_spline_conv->SetMarkerColor(4);

    htmp_meas_conv->SetMarkerStyle(2);
    htmp_sim_spline_conv->SetMarkerStyle(2);    
   
    htmp_meas_conv->SetMarkerSize(1);
    htmp_sim_spline_conv->SetMarkerSize(1);    
    
    htmp_meas_conv->Draw("hist,p");    
    htmp_sim_spline_conv->Draw("hist, p, same");

    htmp_meas_conv->SetTitle("");

    htmp_meas_conv->GetXaxis()->SetTitleFont(22);
    htmp_meas_conv->GetYaxis()->SetTitleFont(22);
    htmp_meas_conv->GetXaxis()->SetTitleOffset(1.2);
    //htmp_meas_conv->GetYaxis()->SetTitleOffset(1.6);
    htmp_meas_conv->GetYaxis()->SetTitleOffset(0.9);
    htmp_meas_conv->GetXaxis()->SetTitleSize(0.05);
    htmp_meas_conv->GetYaxis()->SetTitleSize(0.05);
    
    htmp_meas_conv->GetXaxis()->SetLabelSize(0.05);
    htmp_meas_conv->GetYaxis()->SetLabelSize(0.04);
    
    
    htmp_meas_conv->GetXaxis()->SetTitle("Time [ns]");
    htmp_meas_conv->GetYaxis()->SetTitle("Volt");
    
    htmp_meas_conv->GetXaxis()->CenterTitle();
    htmp_meas_conv->GetYaxis()->CenterTitle();
    
    htmp_meas_conv->GetXaxis()->SetNdivisions(505);
    htmp_meas_conv->GetYaxis()->SetNdivisions(505);

    MyC_tmp->Print("com_sim_meas.pdf");
    


    

    return 0;
}
