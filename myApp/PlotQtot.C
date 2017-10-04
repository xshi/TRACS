{
    // root -l NOirrad_dt200ps_5pF_tNOtrappingns_dz10um_dy5dV5V_0nns_edge_0_rc.hetct.root
    TTree *tree=(TTree *) _file0->Get("edge");

    Int_t icol=1;
    TGraph *gr; 
    for (Int_t iv=50 ; iv<=5 ; iv=iv-5) {
      TString sel=Form( "Vbias==%d" , iv);
      if ( iv==50 )tree->Draw("Qtot:z",sel,"l"); 
      if ( iv!=50 ) tree->Draw("Qtot:z",sel,"lsame"); 
      TString tsel=Form( "%d V" , iv);
      gr=(TGraph *)gPad->GetListOfPrimitives()->Last() ; gr->SetNameTitle(tsel.Data(),tsel.Data()); gr->SetLineColor(icol);icol++; gr->SetLineWidth(2);
    }  

}
