{
   TNtuple t1("t1","Diffusion","x:z");
   t1.SetMarkerStyle(20) ; t1.SetMarkerColor(1);  t1.SetLineColor(1);
   t1.ReadFile("XZs.dat");
   
   TH2D *fr=new TH2D("fr","",100,0,400,100,0.,300); fr->Draw();gStyle->SetOptStat(0);
   t1.Draw("z:x","","lpsame");
   htemp->GetXaxis()->SetTitle("X [#mum]");
   htemp->GetYaxis()->SetTitle("Z [#mum]");
   
   gPad->Print("NerOfCarriers.pdf");
   
   
}
