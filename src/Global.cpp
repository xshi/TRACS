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

/*Global.cpp is used to define some global variables that are used during the whole execution by the threads in different functions and methods.
 *The steering file and the number of threads are declared here
 */

#include <mutex>
#include <TRACSInterface.h>

//Main variables of TRACS to store the induced current during the whole execution
std::valarray<std::valarray <double> > vItotals;
//std::valarray<std::valarray <double> >i_temp;

vector<vector <TH1D*> >  i_ramo_array, i_conv_array;//, i_rc_array;
vector<TH1D*> i_rc_array;

//Define here the steering file you want to use. Store it in myApp folder.
std::string fnm="MyConfigTRACS";
//For mutex areas
std::mutex mtx;
int num_threads;
std::ofstream fileDiffDrift;



