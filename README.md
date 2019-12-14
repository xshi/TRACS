# TRACS

### Note
Add the avalanche feature in the TRACS found at https://github.com/CERN-SSD/TRACS .   
The purpose of this development is to simulate so-called LGAD, which has "gain layer" beneath the charge collection electrodes. 
One important policy with this work (so far) is that minimizing the change as mush as possible
though re-arrangement of the code structure could fit more in future development.

### 1. Install
login to cepcvtx.    
    
user@cepcvtx:~$ cd {some directory under which you want to setup the code = Somewhere}  
user@cepcvtx:Somewhere$ git clone https://github.com/rkiuchi/TRACS  
user@cepcvtx:Somewhere$ cd TRACS   

- ###### Necessary libraries, such as FEniCS, ROOT, Eigen, QT4, are assummed to be already installed if we want to run TRACS at the other PC circumstances.   

### 2. Setup
user@cepcvtx:Somewhere/TRACS$ ./setup.sh

### 3. Build Code
user@cepcvtx:Somewhere/TRACS$ make

### 4. Run Template 
user@cepcvtx:Somewhere/TRACS$ cd ./run/template.IR_Bottom.20190901/   
user@cepcvtx:Somewhere/TRACS$ ../../myApp/DoTracsOnly 1 MyConfigTRACS  

- ###### Note that the above "DoTracsOnly" executable requests the form as "DoTracsOnly [number of threads] [config file]". However, number of threads == 1, is assumed for this version, ( the original code as well ) and no gurantee if we try to run with more than one thread.

### 4a. (Optional) Output root files. -- after successfully run "DoTracsOnly"  . 

You can check several distribution in the ROOT files. 
Notice that those are bi-products during development/debugging procedures, 
thus could be changed/modified in future. 

- ###### Weighting/Electric fields ("120" is a setting value in the config file and would be different)
user@cepcvtx:Somewhere/TRACS$ root wf120V   
root>  h_w_u->Draw();          // Weighting potential  
root>  h_w_f_grad->Draw();     // Weighting field  
root>  h_d_f_grad->Draw();     // E-field  
root>  h_d_f_grad_Y->Draw();   // E-field along with Z-axis  
  
- ###### Transient Current  ("120" is a setting value in the config file and would be different)   
user@cepcvtx:Somewhere/TRACS$ root current120V_scan0    
root>  i_total->Draw();          // Total current    
root>  i_init_elec->Draw();      // Contribution from initial electrons     
root>  i_init_hole->Draw();      // Contribution from initial holes    
root>  i_gen_elec->Draw();       // Contribution from secondary electrons    
root>  i_gen_hole->Draw();       // Contribution from secondary holes    
  
    
- ###### Timing information. Generation time of carriers.      
user@cepcvtx:Somewhere/TRACS$ root ncarrier.root    
root>  e_gentime->Draw();          // Generation time of all of electrons   
root>  h_gentime->Draw();          // Generation time of all of holes   
   
- ###### Effective Doping Profile along with detector vertical direction.   
user@cepcvtx:Somewhere/TRACS$ root Neff_dist.root    
root>  neff_dist->Draw();          // Effective doping profile  
  
   
### 4b. (Optional) "Edge_tree" . -- after successfully run "DoTracsOnly" .   
##### user@cepcvtx:Somewhere/TRACS$ ../../myApp/Edge_tree NOirrad_dt0ps_4pF_tNOtrappingns_dz5um_dy5dV20V_0nns_bottom_0_rc.hetct  
##### user@cepcvtx:Somewhere/TRACS$ root NOirrad_dt0ps_4pF_tNOtrappingns_dz5um_dy5dV20V_0nns_bottom_0_rc.hetct.root  
  
root> .x loadlib.c  
root> edge->StartViewer();  
or   
root> TBrowser a;  


- ###### The file name "NOirrad...." would be different, depending on the setting parameters in "MyConfigTRACS"   
- ###### The functionality/usage is not fully tested and that's why it is marked as "Optional" now. (09/02/2019)

### 5. Config File
Here, newly introduced parameters in the config file ("MyConfigTRACS" in the template) are mainly described.  

- ###### SetAvalanche = yes   # yes: set this doping profile.  no: turn-off, even though parameters are given bellow.  

The shape of "Gain Layer" is now assumed as a simple Gaussian for which we need to specify some parameters:  
- ###### DopingProfile_PeakHeight   = 3.0e16       # unit in  cm^-3  
- ###### DopingProfile_PeakPosition = 1.5          # Z (vertical) position from the top. unit in micron-meter   
- ###### DopingProfile_GaussSigma  = 0.3           # unit in micron-meter  

The (carrier loop)calculation time scales with the number of carriers, including the secondary generated ones due to the Avalanche effect.
Thus, stopping the loop calculation iteration, if the number of generated carriers reaches or exceeds the ratio bellow. The ratio
is defined as the number between the total number of carriers generated and the number of initial carriers read from a carrier file.  
- ###### MaxMultiplicationRatio = 10.0         # no unit.   

In addion to above, the effective doping level in the bulk except the "gain layer" is also important.
This is determined from the detector thickness and the full depletion voltage, which is already used in the original TRACS.
Therefore, you need to set the full depletion voltage correctly. ( In fact, the impact ionization effect happens only 
at the high doping area, and less effect if one sets the full depletion voltage == nominal doping level wrongly )   
One of point is that this full depletion voltage is that one without the "gain layer", to estimate normal doping level correctly.
- ###### DepletionVoltage = 10.0 # in volts   

The other factor related with the impact ionization effect is the temperature of the detector. Do not forget about its value.  

### ChangeLog (from the original repo.)

V1.1)  
- ###### Addition of the impact ionization effect (mainly in Carrier.cpp)      
- ###### Carrier container is changed so that secondary carriers can be stored as well.      
- ###### Note that, above treatment might be dangerous because the memory(stack) overflow can be happen.   
- ###### Introducing the effective doping profile option for the avalanche region.   
- ######    

### Reference 
- #### Impact ionization effect  

The model and parameters are taken from :   
M. Valdinoci, D. Ventura, M. C. Vecchi, M. Rudan, G. Baccarani, F. Illien, A. Stricker, L. Zullino, "Impact-ionization in silicon at large operating temperature", SISPAD '99, Sept. 6-8, 1999, Kyoto, Japan.


