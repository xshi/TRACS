#include "Source.h"

// Static member
const double Source::_elementary_charge = 1.60217662e-19;    // [ C ]
const double Source::_vacuum_permittivity = 8.85418782e-12;  // [ F/m ]
const double Source::_relative_permittivity_silicon = 11.9;  // no unit.


void Source::save_Neff_dist(double zmin, double zmax)
{
    constexpr int nstep = 500;
    
    TFile *fout = new TFile("Neff_dist.root","RECREATE");
    TH1D hist = TH1D("neff_dist", "Effective Doping Profile", nstep, zmin, zmax);
    
    for( int i=0; i<nstep; i++)
    {            
        double tmp_z = zmin + (zmax-zmin)/static_cast<double>(nstep) * i;
        
        auto gauss_term1 = _peak_height[0] * std::exp(- std::pow((_peak_position[0] - tmp_z), 2.0)/(2*std::pow(_gauss_sigma[0], 2.0))) ;  // [ /cm^3 ]
        auto gauss_term2 = _peak_height[1] * std::exp(- std::pow((_peak_position[1] - tmp_z), 2.0)/(2*std::pow(_gauss_sigma[1], 2.0))) ;  // [ /cm^3 ]        

        auto permittivity = _vacuum_permittivity * _relative_permittivity_silicon * 1e-2;  // [ F/cm ]        
        auto base_term = ( _f_poisson * 1e4*1e4 ) * permittivity / _elementary_charge ;  // [ /cm^3 ]

        hist.SetBinContent(i+1, std::abs(base_term)+std::abs(gauss_term1)+std::abs(gauss_term2) );
    }

    hist.Write();    
    fout->Close();
}


void Source::eval(Array<double>& values, const Array<double>& x) const
{
    if (NeffApproach == "Triconstant") 
    {
        /*
         * 3 ZONE constant space distribution
         *
         * We define here a Neff distribution consisting in 3 different zones
         * each zone is defined as a constant within the given region.
         * It uses all but the last parameter (y3 = Neff(z3)). It takes zX 
         * values as boundaries of the zones and the three first yX as the 
         * value of Neff in each region
         *
         * Even though a function like this is generally not continuous, we 
         * add the hyperbolic tangent bridges to ensure not only continuity 
         * but also derivability.
         *
         */
        double neff_1 = y0;
        double neff_2 = y1;
        double neff_3 = y2;
        
        // For continuity and smoothness purposes
        double bridge_1 = tanh(1000*(x[1]-z0)) - tanh(1000*(x[1]-z1));
        double bridge_2 = tanh(1000*(x[1]-z1)) - tanh(1000*(x[1]-z2));
        double bridge_3 = tanh(1000*(x[1]-z2)) - tanh(1000*(x[1]-z3));
        
        double neff = 0.5*((neff_1*bridge_1)+(neff_2*bridge_2)+(neff_3*bridge_3));
        values[0] = neff*0.00152132;
    }
    else if (NeffApproach == "Linear") 
    {
        /*
         * 1 ZONE approximatin
         *
         * First aproximation to the after-irradiation space charge distribution
         * Consists on a simple straight line defined by the points (z0, y0) and 
         * (z3, y3) and neglects the rest of the values.
         *
         */
        
        double neff = ((y0-y3)/(z0-z3))*(x[1]-z0) + y0;
        values[0] = neff*0.00152132;
    }
    else if ( NeffApproach == "AvalancheMode" )
    {
        /*
         * Gaussian approximation. 
         * Pedestal part is "_f_poisson"  parameter
         */
        
        auto permittivity = _vacuum_permittivity * _relative_permittivity_silicon;  // [ F/m ]
        
        auto poisson_term1 =  ((std::signbit(_f_poisson)== false) ? +1.0 : -1.0) * ( _elementary_charge * _peak_height[0] * 1e6 / permittivity );  // [ V/m/m ]
        auto poisson_term_unit_in_microm1 = poisson_term1 * 1e-12;  // [ V/um/um ]

        auto poisson_term2 =  ((std::signbit(_f_poisson)== false) ? +1.0 : -1.0) * ( _elementary_charge * _peak_height[1] * 1e6 / permittivity );  // [ V/m/m ]
        auto poisson_term_unit_in_microm2 = poisson_term2 * 1e-12;  // [ V/um/um ]
                
        auto gauss_term1 = poisson_term_unit_in_microm1 * std::exp(- std::pow((_peak_position[0] - x[1]), 2.0)/(2*std::pow(_gauss_sigma[0], 2.0))) ;
        auto gauss_term2 = poisson_term_unit_in_microm2 * std::exp(- std::pow((_peak_position[1] - x[1]), 2.0)/(2*std::pow(_gauss_sigma[1], 2.0))) ;
        
        
        auto base_term = _f_poisson ;  // [ V/um/um ]
        
        values[0] = base_term + gauss_term1 + gauss_term2;
    }
    else 
    {
        /*
         * 3 ZONE space distribution
         *
         * It consists in 3 different straight lines corresponding to 3 different
         * charge distributions. It uses all 8 parameters to compute the Neff.
         *
         * Continuity is assumed as straight lines have common points, continuity 
         * is ensured by the hyperbolic tangent bridges
         */
        double neff_1 = ((y0-y1)/(z0-z1))*(x[1]-z0) + y0;
        double neff_2 = ((y1-y2)/(z1-z2))*(x[1]-z1) + y1;
        double neff_3 = ((y2-y3)/(z2-z3))*(x[1]-z2) + y2;
        
        // For continuity and smoothness purposes
        double bridge_1 = tanh(1000*(x[1]-z0)) - tanh(1000*(x[1]-z1));
        double bridge_2 = tanh(1000*(x[1]-z1)) - tanh(1000*(x[1]-z2));
        double bridge_3 = tanh(1000*(x[1]-z2)) - tanh(1000*(x[1]-z3));
        
        double neff = 0.5*((neff_1*bridge_1)+(neff_2*bridge_2)+(neff_3*bridge_3));
        values[0] = neff*0.00152132;
        
    }
    // Fix units from the PdC version
    
    
}
