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

/************************************CARRIER***********************************
 *
 * One of the most important classes fo the program, used for instansciate every carrier read from the file.
 * Each carrier is treated independently, calculating its position, diffusion and drifting, then, its contribution to the global induced current is taken into account.
 *
 * Carriers have different properties like carrier_type (e or h), trapping time, etc...
 *
 */

#include "../include/Carrier.h"

std::mutex mtn;

/*
 * Constructor for Carrier.cpp that sets and stores the values given in their respective places.
 *
 *
 *
 */
/**
 * @param carrier_type
 * @param q
 * @param x_init
 * @param y_init
 * @param detector
 * @param gen_time
 */
Carrier::Carrier( char carrier_type, double q,  double x_init, double y_init , SMSDetector * detector, double gen_time = 1.e-9):

								_carrier_type(carrier_type), // Charge carrier(CC)  type. Typically  electron/positron
								_q(q), //Charge in electron units. Always positive.
								_dy(0.),
								_dx(0.),
								_gen_time(gen_time), // Instant of CC generation
								_e_field_mod(0.),
								_detector(detector), // Detector type and characteristics
								//	_electricField(_detector->get_d_f_grad(),
								//	_weightingField(_detector->get_w_f_grad(),
								_myTemp(_detector->get_temperature()), // Temperature of the diode
								_drift(_carrier_type, detector->get_d_f_grad(), _myTemp, _detector->diffusionON(), _detector->get_dt()), // Carrier Transport object
								_mu(_carrier_type, _myTemp),// Mobility of the CC
								_trapping_time(_detector->get_trapping_time()),
								diffDistance(0.),
								_crossed(false)

{

	_x[0] = x_init; // Starting horizontal position
	_x[1] = y_init; // Starting vertical position

	if (_carrier_type == 'e')
	{ // If electron-like
		_sign = -1; // Negative charge
	}
	else
	{ // it's hole-like
		_sign = 1; // Positive charge
	}
}


/*Diffusion method: When calling obtain_mobility method, the _mu object has been instanciated knowing whether it is a H or an E.
 * Thus, both obtain_mobility are the same method but with different variables.
 * The result of the diffusion is added to both coordenates Z(Y) and X modifying the position of the carrier.
 *
 */
/**
 *
 * @param dt
 */
void Carrier::calculateDiffusionStep(double dt){

	TRandom3 Rand(0);

	diffDistance = pow(2*_mu.obtain_mobility(_e_field_mod)*kB*_myTemp/(ECH)*dt,0.5);
	_dx = diffDistance * Rand.Gaus(0,1);
	_dy = diffDistance * Rand.Gaus(0,1);

	_x[0] += _dx;
	_x[1] += _dy;

	/*Becker Thesis, pag. 50. Lutz book, pag 18. Ejercicio de Lutz, pag. 36.
	diffDistance = pow(2*(1440*cm*cm/(V*s))*kB*_myTemp/(ECH)*dt,0.5);

	Other approach could be following MCTSi, Tim Janssen
	diff = sqrt(2*coefficientElectron*dt) * sin(2*pi*gRandom.Uniform()) * sqrt(-2*log(gRandom.Uniform()));*/


}


/*
 ******************** CARRIER DRIFT SIMULATION METHOD**************************
 *
 * Simulates how the CC drifts inside the detector in the 
 * desired number of steps
 *
 */
/**
 *
 * @param dt
 * @param max_time
 * @param x_init
 * @param y_init
 * @return
 */
std::valarray<double> Carrier::simulate_drift(double dt, double max_time, double x_init/*y_pos*/, double y_init /*z_pos*/ )
{
	_x[0] = x_init;
	_x[1] = y_init;
	//std::ofstream fileDiffDrift;

	bool regularCarrier = true;
	double tDiff=0.;
	double tDep = 0;

	// get number of steps from time
	int max_steps = (int) std::floor(max_time / dt);
	std::valarray<double>  i_n(max_steps); // valarray to save intensity
	runge_kutta4<std::array< double,2>> stepper;
	// wrapper for the arrays using dolphin array class
	Array<double> wrap_x(2, _x.data());
	Array<double> wrap_e_field(2, _e_field.data());
	Array<double> wrap_w_field(2, _w_field.data());


	/*Carrier is in NO depleted area*/
	if ( (_x[1] > _detector->get_depletionWidth()) && (_x[1] < _detector->get_y_max()) && (_x[0] > _detector->get_x_min()) && (_x[0] < _detector->get_x_max()) && (_detector->diffusionON()) ){
		regularCarrier = false;
		//Four times the trapping time represent almost 100% of the signal.
		while( (tDiff < (4 * _trapping_time)) &&  (tDiff<(max_time))  ){
			_e_field_mod = 0;
			calculateDiffusionStep(dt); //Carrier movement due to diffusion
			tDiff += dt;
			if ((_x[1] < _detector->get_depletionWidth()) && (_x[1] < _detector->get_y_max()) && (_x[0] > _detector->get_x_min()) && (_x[0] < _detector->get_x_max())){
				regularCarrier = true;
				//To count carriers passing to the depleted region
				_crossed = true;

				break;
			}

		}

	}
	/*End NO depleted area*/
	if ((regularCarrier) && (_x[1] <= _detector->get_depletionWidth())){

		int it0 = ( _detector->diffusionON() ) ? TMath::Nint( (_gen_time + tDiff)/dt ) : TMath::Nint( _gen_time/dt ) ;
		//bool saleOut = false;

		//fileDiffDrift << "New carrier at: "<< "x= " <<_x[0] << "; z= " << _x[1] << std::endl;
		for ( int i = it0 ; i < max_steps; i++)
		{

			if (_detector->is_out(_x)) // If CC outside detector
			{
				i_n[i] = 0;

				//Take into account if it is not depleted area.
				//Check if is inside depletion in less than n*trapping time.
				//If yes, that particles starts to feel diffusion and electric field.
				break; // Finish (CC gone out)
			}
			else
			{


				//safeRead.lock();
				_detector->get_d_f_grad()->eval(wrap_e_field, wrap_x);
				_detector->get_w_f_grad()->eval(wrap_w_field, wrap_x);
				//safeRead.unlock();

				_e_field_mod = sqrt(_e_field[0]*_e_field[0] + _e_field[1]*_e_field[1]);

				i_n[i] = _q *_sign* _mu.obtain_mobility(_e_field_mod) * (_e_field[0]*_w_field[0] + _e_field[1]*_w_field[1]);

				stepper.do_step(_drift, _x, tDep, dt); //Carrier movement due to drift


				if  (_detector->diffusionON()){
					calculateDiffusionStep(dt); //Carrier movement due to diffusion

				}


				// Trapping effects due to radiation-induced defects (traps) implemented in CarrierColleciton.cpp
			}
			tDep+=dt;
		}


		return i_n;
	}
	return i_n=0.;

}

/************************************************************************
 *************************************************************************
 ***                                                                   ***
 ***                                     ***
 ***                                                                   ***
 *************************************************************************
 *************************************************************************/

/*
 * Getter for the type of the CC (electro / hole)
 */
/**
 *
 * @return
 */
char Carrier::Carrier::get_carrier_type()
{
	return _carrier_type; // electron or hole
}

/*
 * Getter for the position of the CC
 */

std::array< double,2> Carrier::get_x()
{
	return _x;
}

/*
 * Getter for the charge of the CC
 */

double Carrier::get_q()
{
	return _q;
}

double Carrier::get_diffDistance(){

	return diffDistance;
}

bool Carrier::crossed(){

	return _crossed;
}


/*
 ********************** DESTRUCTOR OF THE CLASS CARRIER	**************************
 */
Carrier::~Carrier()
{

}

/*
 **************************** PARALLEL **********************************
 */

/*
 * Copy constructor
 */
Carrier::Carrier(const Carrier& other)
{
	_carrier_type = other._carrier_type;
	_q = other._q;
	_gen_time = other._gen_time;
	_x = other._x; 
	_e_field = other._e_field; 
	_w_field = other._w_field;
	_e_field_mod = other._e_field_mod;
	_sign = other._sign; 
	_detector = other._detector;
	_myTemp = other._myTemp;
	_drift = other._drift;
	_mu = other._mu;
	_trapping_time = other._trapping_time;
	_dx = other._dx;
	_dy = other._dy;
	diffDistance = other.diffDistance;
	_crossed = other._crossed;
	std::lock_guard<std::mutex> lock(other.safeRead);
}

/*
 * Copy assignment
 */
Carrier& Carrier::operator = (const Carrier& other) 
{
	std::lock(safeRead, other.safeRead);
	std::lock_guard<std::mutex> self_lock(safeRead, std::adopt_lock);
	std::lock_guard<std::mutex> other_lock(other.safeRead, std::adopt_lock);
	_carrier_type = other._carrier_type;
	_q = other._q;
	_gen_time = other._gen_time;
	_x = other._x; 
	_e_field = other._e_field; 
	_w_field = other._w_field;
	_e_field_mod = other._e_field_mod;
	_sign = other._sign; 
	_detector = other._detector;
	_myTemp = other._myTemp;
	_drift = other._drift;
	_trapping_time = other._trapping_time;
	_dx = other._dx;
	_dy = other._dy;
	diffDistance = other.diffDistance;
	_crossed = other._crossed;
	return *this;
}


/*
 * Move constructor
 */
Carrier::Carrier(Carrier&& other)
{
	_carrier_type = std::move(other._carrier_type);
	_q = std::move(other._q);
	_gen_time = std::move(other._gen_time);
	_x = std::move(other._x); 
	_e_field = std::move(other._e_field); 
	_w_field = std::move(other._w_field);
	_e_field_mod = std::move(other._e_field_mod);
	_sign = std::move(other._sign); 
	_detector = std::move(other._detector);
	_myTemp = std::move(other._myTemp);
	_drift = std::move(other._drift);
	_mu = std::move(other._mu);
	_trapping_time = std::move(other._trapping_time);
	_dx = std::move(other._dx);
	_dy = std::move(other._dy);
	diffDistance = std::move(other.diffDistance);
	_crossed = std::move(other._crossed);
	std::lock_guard<std::mutex> lock(other.safeRead);
}

/*
 * Move assignment
 */
Carrier& Carrier::operator = ( Carrier&& other) 
{
	std::lock(safeRead, other.safeRead);
	std::lock_guard<std::mutex> self_lock(safeRead, std::adopt_lock);
	std::lock_guard<std::mutex> other_lock(other.safeRead, std::adopt_lock);
	_carrier_type = std::move(other._carrier_type);
	other._carrier_type = '\0';
	_q = std::move(other._q);
	other._q = 0;
	_gen_time = std::move(other._gen_time);
	other._gen_time = 0;
	_x = std::move(other._x); 
	other._x = {0,0};
	_e_field = std::move(other._e_field); 
	other._e_field = {0,0};
	_w_field = std::move(other._w_field);
	other._w_field = {0,0};
	_e_field_mod = std::move(other._e_field_mod);
	other._e_field_mod = 0;
	_sign = std::move(other._sign); 
	other._sign = 0;
	_detector = std::move(other._detector);
	other._detector = NULL;
	_myTemp = std::move(other._myTemp);
	other._myTemp = 0;
	_drift = std::move(other._drift);
	_mu = std::move(other._mu);
	_trapping_time = std::move(other._trapping_time);
	other._trapping_time = 1e-300;
	_dx = std::move(other._dx);
	other._dx = 0.;
	_dy = std::move(other._dy);
	other._dy = 0.;
	diffDistance = std::move(other.diffDistance);
	other.diffDistance = 0.;
	_crossed = std::move(other._crossed);
	other._crossed = false;
	return *this;
}
