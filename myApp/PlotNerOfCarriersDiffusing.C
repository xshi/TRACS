{
   TNtuple t1("t1","Diffusion","z:N");
   TNtuple t2("t2","Diffusion","z:N");
   TNtuple t3("t3","Diffusion","z:N");
   t1.SetMarkerStyle(20) ; t1.SetMarkerColor(1);  t1.SetLineColor(1);
   t2.SetMarkerStyle(24) ; t2.SetMarkerColor(2);  t2.SetLineColor(2);
   t3.SetMarkerStyle(21) ; t3.SetMarkerColor(3);  t3.SetLineColor(3);
   t1.ReadFile("Diffusion.dat");
   t1.Draw("N:z","","lp");
   
   t2.ReadFile("Diffusion2.dat");
   t2.Draw("N:z","","lpsame");
   
   t3.ReadFile("Diffusion3.dat");
   t3.Draw("N:z","","lpsame");
   
   htemp->GetXaxis()->SetTitle("Z [#mum]");
   htemp->GetYaxis()->SetTitle("Fraction of carriers entering depleted region [%]");
   
   gPad->Print("NerOfCarriers.pdf");
   
   
}
