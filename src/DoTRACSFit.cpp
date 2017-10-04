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

/************************************DoTRACSFit***********************************
 *
 * TRACS can run multiple modes of operation, one of them is fitting some parameters related to the characterization of the detector.
 * This file contains the main function which manages the files passed by the user and the Minuit callings.
 * It implements operator as well, the steering method of the loop used in the minimization.
 *
 * There is a common part in every main function used in TRACS which correspond to the threads initialization.
 *
    DoTRACSFit MeasurementFile TRACS.conf "Vbias==200 && Tset == 20"

 */

//#include <TApplication.h>
#include <Minuit2/FunctionMinimum.h>
#include <Minuit2/MnUserParameterState.h>
#include <Minuit2/MnUserParameters.h>
#include <Minuit2/MnMachinePrecision.h>
#include <Minuit2/MnPrint.h>
#include <Minuit2/MnMigrad.h>
#include <Minuit2/MnMinos.h>
#include <Minuit2/MnSimplex.h>
#include <Minuit2/MnPlot.h>
#include <Minuit2/MinosError.h>
#include <Minuit2/FCNBase.h>

#include <boost/asio.hpp>

#include <TRACSFit.h>
#include <TRACSInterface.h>
#include <TString.h>
#include <stdio.h>

std::vector<TRACSInterface*> TRACSsim;
std::vector<std::thread> t;
boost::posix_time::time_duration total_timeTaken ;
TRACSFit *fit ;
std::vector<double> vector_zValues, vector_voltValues, vector_yValues;
std::vector<std::string> carrierThread_fileNames;
std::string scanType;
double dTime;
double max_time;

void spread_into_threads();

using namespace ROOT::Minuit2;

int main( int argc, char *argv[]) {

	std::string carrierFile;
	std::valarray<std::valarray<double>> i_total;
	TH1D *i_rc;
	int counted_numLines = 0;
	int number_of_lines = 0;
	int nns, count2;
	std::string line;
	std::string temp;


	//TApplication theApp("DoTRACSFit", 0, 0);

	//Number of threads
	num_threads = atoi(argv[1]);

	//Measurement file
	TString FileMeas = TString( argv[2] ) ;
	//fnm = FileMeas.Data();

	//Configuration file
	TString FileConf = TString( argv[3] ) ;
	fnm = argv[3];
	//std::string lfnm(argv[3]) ;
	//fnm = lfnm;

	//Restrictions for fits
	TString how="";
	if (argc>2) how = TString( argv[4] ) ;



	//file creation
	for (int i = 0; i < num_threads; ++i) {
		temp = std::to_string(i);
		carrierThread_fileNames.push_back("file" + temp);
	}

	utilities::parse_config_file(fnm, carrierFile);

	std::ifstream infile(carrierFile);
	while (std::getline(infile, line))
	{
		++number_of_lines;
	}
	int resto_carriers = number_of_lines % num_threads;
	int carriers_per_thr = number_of_lines / num_threads;
	infile.close();

	counted_numLines = counted_numLines - resto_carriers;
	std::ifstream in(carrierFile);
	for (int i = 0; i < num_threads; ++i) {
		std::ofstream out(carrierThread_fileNames[i]);
		while (in){
			counted_numLines++;
			std::getline(in, line);
			out << line << std::endl;
			if (counted_numLines == carriers_per_thr){
				out.close();
				counted_numLines = 0;
				break;
			}

		}
	}
	in.close();



	spread_into_threads();
	double timeSteps = (int) std::floor(max_time / dTime);

	if (scanType == "edge"){
		i_rc_array.resize(vector_voltValues.size());
		vItotals.resize(vector_voltValues.size()* vector_zValues.size());
		for (int i = 0; i < vector_voltValues.size() ; i++)
			vItotals[i].resize(timeSteps);
		//i_rc_array[i].resize(vector_zValues.size());
	}


	if (scanType == "top" || scanType == "bottom"){
		i_rc_array.resize(vector_voltValues.size());
		vItotals.resize(vector_voltValues.size()* vector_yValues.size());
		for (int i = 0; i < vector_voltValues.size() ; i++)
			//i_rc_array[i].resize(vector_yValues.size());
			vItotals[i].resize(timeSteps);

	}

	for (int i = 0 ; i < vItotals.size() ; i++){
		for (int j = 0 ; j < vItotals[i].size() ; j++)
			vItotals[i][j] = 0;
	}


	TRACSsim.resize(num_threads);
	t.resize(num_threads);
	for (int i = 0; i < num_threads; ++i) {
		t[i] = std::thread(call_from_thread, i, std::ref(carrierThread_fileNames[i]), vector_zValues, vector_yValues, vector_voltValues);
	}
	for (int i = 0; i < num_threads; ++i) {
		t[i].join();
	}


	fit = new TRACSFit( FileMeas, FileConf , how ) ;

	//Define parameters and their errors to Minuit

	vector<Double_t> parIni = TRACSsim[0]->get_NeffParam();
	parIni.push_back(TRACSsim[0]->get_fitNorm());
	Int_t parIniSize = parIni.size() ;

	vector<Double_t> parErr(parIniSize, 1.) ;

	//Pass parameters to Minuit
	MnUserParameters upar(parIni,parErr) ;
	for ( int i=0 ; i < parIniSize ; i++ ) {
		char pname[parIniSize]; sprintf( pname , "p%d" , i);
		upar.SetName( i , pname );
	}

	//Fix parameters
	upar.Fix(0) ;
	upar.Fix(1) ; upar.Fix(2) ;
	//upar.Fix(3) ;
	upar.Fix(4) ; upar.Fix(5); upar.Fix(6) ; upar.Fix(7);
	//upar.Fix(8);

	std::cout << "=============================================" << std::endl;
	std::cout<<"Initial parameters: "<<upar<<std::endl;
	std::cout << "=============================================" << std::endl;

	//Do the minimization

	MnMigrad mn( *fit , upar ) ;
	FunctionMinimum min = mn() ;

	//Status report
	if (min.IsValid()) std::cout << "Fit success"         << std::endl ;
	else               std::cout << "Fit failed"   << std::endl ;
	std::cout << "Total time: " << total_timeTaken.total_seconds() << std::endl ;
	std::cout << "MINIMIZATION OUTCOME: " <<  min  << std::endl ;

	//Release parameter 3 and 0, minimize again
	//upar.SetValue(0, min.UserState().Value(0) ) ;
	//upar.SetValue(3, min.UserState().Value(3) ) ;
	//upar.Release(0) ; upar.Release(3) ;
	//MnMigrad mnr( *fit , upar ) ;
	//min = mnr() ;

	//Status report
	//if (min.IsValid()) std::cout << "Fit success"         << std::endl ;
	//else               std::cout << "Fit failed"   << std::endl ;
	//std::cout << "Total time: " << total_timeTaken.total_seconds() << std::endl ;
	//std::cout << "MINIMIZATION OUTCOME: " <<  min  << std::endl ;

	//Get the fitting parameters
	for (uint i=0; i<parIniSize;i++) {
		parIni[i]=min.UserState().Value(i);
		parErr[i]=min.UserState().Error(i);
	}

	//Calculate TCT pulses with the fit output parameters
	for (int i = 0; i < num_threads; ++i) {
		t[i] = std::thread(call_from_thread_FitPar, i, std::ref(carrierThread_fileNames[i]), vector_zValues, vector_yValues, vector_voltValues ,parIni);
	}

	for (int i = 0; i < num_threads; ++i) {
		t[i].join();
	}

	//Dump tree to disk
	TFile fout("output.root","RECREATE") ;
	TTree *tout = new TTree("edge","Fitting results");

	TMeas *emo = new TMeas( );
	emo->Nt   = TRACSsim[0]->GetnSteps() ;
	emo->volt = new Double_t [emo->Nt] ;
	emo->time = new Double_t [emo->Nt] ;
	emo->Qt   = new Double_t [emo->Nt] ;

	// Create branches
	tout->Branch("raw", &emo,32000,0);

	//Read RAW file
	TRACSsim[0]->DumpToTree( emo , tout ) ;

	//TRACSsim[0]->GetTree( tsim );
	fout.Write();
	delete tout ;
	fout.Close();
	delete emo ;

	//Clean
	for (uint i = 0; i < TRACSsim.size(); i++)	{
		delete TRACSsim[i];
	}

	delete fit;
	std::quick_exit(1);
}

//_____________________________________________________________________
/**
 *
 * @param par
 * @return
 */
Double_t TRACSFit::operator() ( const std::vector<Double_t>& par  ) const {

	static int icalls ;
	boost::posix_time::ptime start = boost::posix_time::second_clock::local_time();

	for (int i = 0; i < num_threads; ++i) {
		t[i] = std::thread(call_from_thread_FitPar, i, std::ref(carrierThread_fileNames[i]), vector_zValues, vector_yValues, vector_voltValues ,par);
	}

	for (int i = 0; i < num_threads; ++i) {
		t[i].join();
	}

	Double_t chi2 = fit->LeastSquares( ) ;
	boost::posix_time::ptime end = boost::posix_time::microsec_clock::local_time();
	boost::posix_time::time_duration timeTaken = end - start;
	total_timeTaken += timeTaken;

	std::cout << "-------------------------------------------------------------------------------------> " << std::endl;
	std::cout << "----------------------------> Time taken for chi2 calculation (milliseconds): " << timeTaken.total_milliseconds() << std::endl;
	std::cout << "-------------------------------------------------------------------------------------> " << std::endl;
	std::cout << "----------------------------> icalls="<<icalls<<" chi2=" << chi2 << "\t" ;
	for (uint ipar=0 ; ipar<par.size() ; ipar++) std::cout << "p["<<ipar<<"]="<<par[ipar]<<"\t" ; std::cout << std::endl;
	std::cout << "-------------------------------------------------------------------------------------> " << std::endl;
	icalls++;

	for (uint i = 0; i < TRACSsim.size(); i++)	{
		delete TRACSsim[i];
	}

	return chi2 ;

}

void spread_into_threads(){

	double vInit, deltaV, vMax, v_depletion, zInit, zMax, deltaZ, yInit, yMax, deltaY;
	int zSteps, ySteps, vSteps, seeker;
	bool lastElement = false;

	uint index = 1;

	//Reading file to get values that will determine scan vectors
	utilities::parse_config_file(fnm, scanType, vInit, deltaV, vMax, v_depletion, zInit, zMax, deltaZ, yInit, yMax, deltaY, dTime, max_time);

	//Reserve to 300 since is a common measurement for a detector Z dimension.
	//Reserve does not implies anything, better for performance.
	//For sure there is at least one value of Z coordinates in an edge-TCT scan
	vector_zValues.reserve(300);
	vector_zValues.push_back(zInit);
	seeker = zInit + deltaZ;

	//Filling Z scan vector
	while (seeker < zMax){
		vector_zValues.push_back(seeker);
		++index;
		seeker = seeker + deltaZ;
		if (seeker > zMax){
			vector_zValues.push_back(zMax);
			lastElement = true;
		}

	}
	if (seeker >= zMax && !lastElement && zMax !=zInit) vector_zValues.push_back(zMax);
	lastElement = false;

	//Reserve to 100 for differente voltages.
	//There will be at least one bias voltage
	vector_voltValues.reserve(300);
	vector_voltValues.push_back(vInit);
	seeker = vInit + deltaV;

	//Filling voltage scan vector
	while (seeker < vMax){
		vector_voltValues.push_back(seeker);
		++index;
		seeker = seeker + deltaV;
		if (seeker > vMax) {
			vector_voltValues.push_back(vMax);
			lastElement = true;
		}
	}
	if (seeker >= vMax && !lastElement && vMax != vInit) vector_voltValues.push_back(vMax);
	lastElement = false;

	//Reserve to 400 since is a common measurement for a detector 4 dimension.
	//For sure there is at least one value of Z coordinates in an edge-TCT scan
	vector_yValues.reserve(400);
	vector_yValues.push_back(yInit);
	seeker = yInit + deltaY;

	//Filling Y scan vector
	while (seeker < yMax){
		vector_yValues.push_back(seeker);
		++index;
		seeker = seeker + deltaY;
		if (seeker > yMax) {
			vector_yValues.push_back(yMax);
			lastElement = true;
		}
	}
	if (seeker >= yMax && !lastElement && yMax != yInit) vector_yValues.push_back(yMax);
	lastElement = false;

	//Checking for a right user configuration
	if ((vector_yValues.size() > 1) && (vector_zValues.size() > 1)){
		std::cout << "TRACS does not allow this configuration, please, choose zScan OR yScan" << std::endl;
		std::quick_exit(1);
	}



}


