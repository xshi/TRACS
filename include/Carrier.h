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

#ifndef CARRIER_H
#define CARRIER_H

#include  <valarray>
#include  <mutex>
#include <iostream>
#include <fstream>

#ifndef Q_MOC_RUN  // See: https://bugreports.qt-project.org/browse/QTBUG-22829
#include <boost/numeric/odeint/stepper/runge_kutta4.hpp>
#endif
#include <TRandom3.h>

#include <CarrierTransport.h>
#include <SMSDetector.h>
#include <Constants.h>
#include <Global.h>

using namespace boost::numeric::odeint;
//using namespace dolfin;

/************************************CARRIER***********************************
 *
 * Used for instansciating every carrier read from the file.
 * Each carrier is treated independently, calculating its position, diffusion and drifting, then, its contribution to the global induced current is taken into account.
 *
 * Carriers have different properties like carrier_type (e or h), trapping time, etc...
 *
 */

class Carrier
{
	//It is not very clear why, but this class does not accept more variables.
	//The size is limited somewhere by someone, maybe the carrier vector list is touching the maximum space of the program's stack...
private:

	char _carrier_type;
	double _q; // charge
	double _dy;
	double _dx; //to calculate diffusion step
	double _gen_time; // instant of generation of the carrier
	std::array< double,2> _x; // carrier position array
	std::array< double,2> _e_field; // electric field at the carrier position
	std::array< double,2> _w_field; // weighting field at the carrier positions
	double _e_field_mod;
	int _sign; // sign to describe if carrier moves in e field direction or opposite
	mutable std::mutex safeRead;
	SMSDetector * _detector;
	double _myTemp; // Temperature of the detector
	DriftTransport _drift;
	JacoboniMobility _mu;
	double _trapping_time;
	double diffDistance;
	bool _crossed;
	TRandom3 gRandom;

public:



	Carrier( char carrier_type, double q, double x_init, double y_init, SMSDetector * detector, double gen_time);
	Carrier(Carrier&& other); // Move declaration
	Carrier& operator = (Carrier&& other); // Move assignment
	Carrier(const Carrier& other); // Copy declaration
	Carrier& operator = (const Carrier& other); // Copy Assignment
	~Carrier();

	char get_carrier_type();
	//		double get_gen_time();
	//    std::array< double,2> get_e_field;
	//    std::array< double,2> get_w_field;
	//		double get_e_field_mod;
	//    int get_sign;
	//    SMSDetector *get_detector;
	//    double get_myTemp;
	//    DriftTransport get_drift;
	//    JacoboniMobility get_mu;
	//    double get_trapping_time;
	//
	std::array< double,2> get_x();
	double get_q();
	double get_diffDistance();
	bool crossed();

	//double get_diffH(){return diffH;};

	void calculateDiffusionStep(double dt);

	//void calculateDiffusionH(double dt);
	std::valarray<double> simulate_drift( double dt, double max_time);
	std::valarray<double> simulate_drift(double dt, double max_time, double x_init, double y_init );
};

#endif // CARRIER_H
