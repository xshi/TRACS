/*
 * @ Copyright 2014-2017 CERN and Instituto de Fisica de Cantabria - Universidad de Cantabria. All rigths not expressly granted are reserved [tracs.ssd@cern.ch]
 * This file is part of TRACS.
 *
 * TRACS is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public License as published by the Free Software Foundation,
 * either version 3 of the Licence.
 *
 * TRACS is distributed in the hope that it will be useful , but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * See the GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License along with TRACS. If not, see <http://www.gnu.org/licenses/>
 */

#include <cmath>
#include <array>

#include <TH1D.h>
#include <TFile.h>

#include <dolfin.h>

using namespace dolfin;


/*
 * SOURCE TERM
 *
 * This class contains the source term for the poisson equation
 * Space charge distribution is parametrized here.
 *
 *
 */
 class Source : public Expression
 {
 public:


	 //Declaration and default values. They will taken from steering file.
	// Concentration in transition and extremal points
	double y0 = -25.; // Neff(z0)
	double y1 = 0.02; // Neff(z1)
	double y2 = y1+y1*10; // Neff(z2)
	double y3 = 33; // Neff(z3)

	// Define diferent zones in depth
	double z0 = 0.;
	double z1 = 120.;
	double z2 = 220.;
	double z3 = 300.;
	std::string NeffApproach = "Triconstant";

    
	void eval(Array<double>& values, const Array<double>& x) const;

    std::string get_NeffApproach() const
    {
        return NeffApproach;
    }        
	void set_NeffApproach(std::string Neff_type)
	{
		NeffApproach = Neff_type;
	}

    void set_avalanche_doping_param( const std::array<std::array<double, 3>, 2>& param )
    {
        _peak_height[0]    = param[0][0];
        _peak_position[0] = param[0][1];
        _gauss_sigma[0]  = param[0][2];

        _peak_height[1]    = param[1][0];
        _peak_position[1] = param[1][1];
        _gauss_sigma[1]  = param[1][2];
    }
    void set_bulk_doping_param( const double& param )
    {
        _f_poisson = param;
    }
    
	void set_y0(double newValue)
	{
		y0 = newValue;
	}

	void set_y1(double newValue)
	{
		y1 = newValue;
	}

	void set_y2(double newValue)
	{
		y2 = newValue;
	}

	void set_y3(double newValue)
	{
		y3 = newValue;
	}

	void set_z0(double newValue)
	{
		z0 = newValue;
	}

	void set_z1(double newValue)
	{
		z1 = newValue;
	}

	void set_z2(double newValue)
	{
		z2 = newValue;
	}

	void set_z3(double newValue)
	{
		z3 = newValue;
	}

    void save_Neff_dist(double ymin, double ymax);
        
private:

    // For the effective doping parameters
    std::array<double, 2> _peak_height;
    std::array<double, 2> _peak_position;
    std::array<double, 2> _gauss_sigma;
    double _f_poisson;

    static const double _elementary_charge;
    static const double _vacuum_permittivity;
    static const double _relative_permittivity_silicon; 
    
 };

 
