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


#ifndef SRC_THREADING_H_
#define SRC_THREADING_H_

#include <vector>
#include <TRACSFit.h>

void call_from_thread(int, std::string&, const std::vector<double>& z, const std::vector<double>& y, const std::vector<double>& volt);
void call_from_thread_FitPar(int, std::string&, const std::vector<double>& z, const std::vector<double>& y, const std::vector<double>& volt, const std::vector<Double_t>& par);
void call_from_thread_FitNorm(int, std::string&, const std::vector<double>& z, const std::vector<double>& y, const std::vector<double>& volt, const std::vector<Double_t>& par);

#endif /* SRC_THREADING_H_ */
