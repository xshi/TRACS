2017-07-24_16-14-41_dmil_PIN_8622_W5_E3_4.txt.root
NOirrad_dt50ps_4pF_tNOtrappingns_dz5um_dy5dV5V_0nns_bottom_0_rc_Vdep58_280um_DIFF1.hetct.root
NOirrad_dt50ps_4pF_tNOtrappingns_dz5um_dy5dV5V_0nns_bottom_0_rc_Vdep58_285um_DIFF2.hetct.root
NOirrad_dt50ps_4pF_tNOtrappingns_dz5um_dy5dV5V_0nns_bottom_0_rc_Vdep58_280um_DIFF1.hetct.root
NOirrad_dt50ps_4pF_tNOtrappingns_dz5um_dy5dV5V_0nns_bottom_0_rc_Vdep38_285um_DIFF1.hetct.root

root -l 2017-07-24_16-14-41_dmil_PIN_8622_W5_E3_4.txt.root NOirrad_dt50ps_4pF_tNOtrappingns_dz5um_dy5dV5V_0nns_bottom_0_rc_Vdep58_285um_DIFF1.hetct.root

root -l 2017-07-24_16-14-41_dmil_PIN_8622_W5_E3_4.txt.root NOirrad_dt50ps_4pF_tNOtrappingns_dz5um_dy5dV5V_0nns_bottom_0_rc_Vdep58_280um_DIFF1.hetct.root NOirrad_dt50ps_4pF_tNOtrappingns_dz5um_dy5dV5V_0nns_bottom_0_rc_Vdep38_285um_DIFF1.hetct.root

TH2D *fr=new TH2D("fr","",40,0,700,40,0,1.2) ; 
fr->Draw()
vtree[0]->Draw("Sum$((volt-BlineMean)*(time>10.53 && time<100))/2.88144:-Vbias","Vbias>-40","lsame")
vtree[0]->Draw("Sum$((volt-BlineMean)*(time>10.53 && time<50))/2.88144:-Vbias","Vbias<-40","lsame")
vtree[1]->Draw("Qtot/0.00145856:Vbias","","lsame")
vtree[2]->Draw("Qtot/0.00145856:Vbias","","lsame")

vtree[0]->Draw("volt-BlineMean:time-10.53","Vbias==-100","l")
vtree[2]->Draw("85*(volt-BlineMean):time-10.","Vbias==100","same")

//This one is GOOD
root -l 2017-07-24_16-14-41_dmil_PIN_8622_W5_E3_4.txt.root NOirrad_dt50ps_4pF_tNOtrappingns_dz5um_dy5dV5V_0nns_bottom_0_rc_Vdep58_280um_DIFF1.hetct.root 
vtree[0]->Draw("(volt-BlineMean)/Q50:time-10.53","Vbias==-100","l")
vtree[1]->Draw("(volt-BlineMean)/Qtot:time-10.","Vbias==100","lsame")

vtree[0]->Draw("(volt-BlineMean)/Q50:time-10.53","Vbias==-200","l") ; vtree[1]->Draw("(volt-BlineMean)/Qtot:time-10.","Vbias==200","lsame")
vtree[0]->Draw("(volt-BlineMean)/Q50:time-10.53","Vbias==-300","l") ; vtree[1]->Draw("(volt-BlineMean)/Qtot:time-10.","Vbias==300","lsame")

