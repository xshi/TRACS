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


#include "CarrierMobility.h"

/*
 ******************************CARRIER MOBILITY*************************
 *
 * This function calculates the value of the mobility using data from 
 * desy experiment.
 *
 * Inputs: Carrier type ('e' for electrons / 'h' for holes)
 * 	   Temperature of the diode
 *
 * 	   - Field at desired point (only to obtain mobility values)
 *
 * "Outputs": Desired mobility.
 *
 * Data Source:  http://www.desy.de/~beckerj/cc/files/Thesis.pdf
 *
 */
/**
 *
 * @param carrier_type
 * @param T
 */
JacoboniMobility::JacoboniMobility( char carrier_type, double T)
{
  _T = T;
  if (carrier_type == 'e') // Electrons
  {
    //_mu0 = 1440.e8 * std::pow(_T/300., -2.260); //um**2/ Vs
    //_vsat = 1.054e11  * std::pow(_T/300., -0.602); //um**2/ Vs
    //_beta = 0.992 * std::pow(_T/300., 0.572); // <100> orientation

	//Extended Canali model. Following Sentaurus manual. Pag 217-224.
    _mu0  = 1417.e8 * std::pow(_T/300., -2.5); //um**2/ Vs
    _vsat = 1.07e11  * std::pow(_T/300., -0.87); //um**2/ Vs
    _beta = 1.109 * std::pow(_T/300., 0.66); // <100> orientation
  }

  //Extended Canali model. Following Sentaurus manual. Pag 217-224.
  else if (carrier_type == 'h') // Holes
  {
    //_mu0 = 474.e8 * std::pow(_T/300., -2.619); //um**2/ Vs
    //_vsat = 0.940e11  * std::pow(_T/300., -0.226); //um**2/ Vs
    //_beta = 1.181 * std::pow(_T/300., 0.633 ); // <100> orientation

	  //Extended Canali model. Following Sentaurus manual. Pag 217-224.
    _mu0  = 470.e8 * std::pow(_T/300., -2.2); //um**2/ Vs
    _vsat = 0.837e11  * std::pow(_T/300., -0.52); //um**2/ Vs
    _beta = 1.213 * std::pow(_T/300., 0.17 ); // <100> orientation
  }
}

/*
 * OBTAIN MOBILITY - Getter method
 *
 * This method provides the value of the mobility using all the 
 * varibles initialized in the invokation.
 *
 */
/**
 *
 * @param e_field_mod
 * @return
 */
double JacoboniMobility::obtain_mobility(double e_field_mod)
{
  return _mu0/std::pow(1.0+std::pow(_mu0*e_field_mod/_vsat,_beta), 1.0/_beta); // mum**2/ Vs
}
/**
 *
 * @return
 */
double JacoboniMobility::obtain_mu0(){

  return _mu0;

}

JacoboniMobility::~JacoboniMobility()
{

}
JacoboniMobility::JacoboniMobility(){}
