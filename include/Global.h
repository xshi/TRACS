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

#ifndef GLOBAL_H
#define GLOBAL_H

#include <vector>
#include <valarray>
#include <mutex>
#include <atomic>

#include <TH1D.h> // 1 Dimensional ROOT histogram


//extern std::vector<std::vector <TH1D*> >  i_ramo_array, i_conv_array;//, i_rc_array;
extern std::vector<TH1D*> i_rc_array, i_ramo_vector, i_conv_vector;
extern int num_threads;
extern std::mutex mtx;
extern std::string fnm;
extern std::string global_TF;
extern std::valarray<std::valarray <double> > vItotals;
//extern std::valarray<std::valarray <double> > i_temp;
extern std::ofstream fileDiffDrift;
//extern std::atomic<double> numberDs;
//extern std::atomic<int> tempNumberDs;
//extern bool printa;

#endif // GLOBAL_H
