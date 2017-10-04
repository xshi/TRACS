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

#ifndef CARRIER_COLLECTION_H
#define CARRIER_COLLECTION_H

#include "Carrier.h"
#include <CarrierMobility.h>

#include <string>
#include <sstream>
#include <fstream>
#include <vector>

#include <QString>
#include <TH2D.h>
#include <TString.h>
#include <TMath.h>
#include <TRandom3.h>

/*
 ***********************************CARRIER COLLECTION***********************************
 *
 *Carrier collection get carriers information from a file and store it in a vector of carriers.
 *This vector of carriers will be used by each of the threads in each of the steps of the TCT scan.
 *
 */

class CarrierCollection
{
private:

	std::vector<Carrier> _carrier_list_sngl;
	SMSDetector * _detector;

	TRandom3 gRandom;



public:
	CarrierCollection(SMSDetector * detector);
	~CarrierCollection();

	double beamy = 0. , beamz = 0.; //Mean position of the injected carriers in detector plane (y,z)

	void add_carriers_from_file(const std::string& filename,const std::string& scanType, double depth);
	void simulate_drift( double dt, double max_time, double shift_x, double shift_y,  std::valarray<double> &curr_elec, std::valarray<double> &curr_hole, int &totalCrosses);

	TH2D get_e_dist_histogram(int n_bins_x, int n_bins_y, TString hist_name = "e_dist", TString hist_title ="e_dist");
	TH2D get_e_dist_histogram(int n_bins_x, int n_bins_y, double shift_x, double shift_y, TString hist_name = "e_dist", TString hist_title ="e_dist");

};


#endif // CARRIER_COLLECTION_H


