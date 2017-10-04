# root -l  irrad_dt200ps_5pF_t5ns_dz3um_dy5dV100V_2nns_edge_0_rc_nodiff.hetct.root irrad_dt200ps_5pF_t5ns_dz3um_dy5dV100V_2nns_edge_0_rc_diff.hetct.root

{
   
   TTree *tnd=(TTree*)_file0->Get("edge"); 
   TTree *td=(TTree*)_file1->Get("edge"); 
   
   tnd->SetLineColor(1) ; tnd->SetMarkerColor(1) ;  tnd->SetMarkerStyle(20) ; 
   td->SetLineColor(2)  ; td->SetMarkerColor(2)  ;  td->SetMarkerStyle(24) ;
   
   tnd->Draw("Sum$(volt*((time>0.5 && time<50))):z","","p")    ;
   td->Draw("Sum$(volt*((time>0.5 && time<50))):z","","psame") ;
   
   
   
   td->Draw("volt:time-tleft" , "", "l") ; 
  
   


}
