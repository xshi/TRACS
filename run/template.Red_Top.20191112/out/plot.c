#include <iostream>


int plot()
{
    gStyle->SetOptStat(000);
    
    std::stringstream hist_title_meas;
    std::stringstream hist_title_sim_spline;
    std::stringstream hist_title_sim;
    std::stringstream hist_title_sim_woshift;    

    hist_title_meas << "i_total_meas" ;
    hist_title_sim_spline << "i_total_sim_spline" ;
    hist_title_sim << "i_total_sim" ;
    hist_title_sim_woshift << "i_total_sim_woshift" ;        

    
    //TFile *f1 = new TFile("fit_current_29.root");
    
    TH1D* htmp_meas = dynamic_cast<TH1D*>(gDirectory->Get(hist_title_meas.str().c_str()));
    TH1D* htmp_sim_spline = dynamic_cast<TH1D*>(gDirectory->Get(hist_title_sim_spline.str().c_str()));
    TH1D* htmp_sim = dynamic_cast<TH1D*>(gDirectory->Get(hist_title_sim.str().c_str()));        
    TH1D* htmp_sim_woshift = dynamic_cast<TH1D*>(gDirectory->Get(hist_title_sim_woshift.str().c_str()));    


    htmp_meas->SetMarkerColor(1);
    htmp_sim_spline->SetMarkerColor(2);
    htmp_sim_woshift->SetMarkerColor(4);

   
    htmp_meas->SetMarkerStyle(20);
    htmp_meas->SetMarkerSize(0.2);

    htmp_sim_spline->SetMarkerStyle(20);
    htmp_sim_spline->SetMarkerSize(0.2);

    htmp_sim_woshift->SetMarkerStyle(20);
    htmp_sim_woshift->SetMarkerSize(0.2);        

    
    htmp_meas->Draw("hist p");
    htmp_sim_woshift->Draw("hist p same");
    htmp_sim_spline->Draw("hist p same");


    return 0;
}
