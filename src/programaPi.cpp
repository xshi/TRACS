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

/*Test program to test the modularization and independency calling of the main function of TRACS*/

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <vector>
#include <thread>


#include <TRACSInterface.h>
#include "../include/Threading.h"


std::vector<TRACSInterface*> TRACSsim;
std::vector<std::thread> t;


int main(int argc, char *argv[])
{
	num_threads = 2;
	TRACSsim.resize(num_threads);
	t.resize(num_threads);

	for (int i = 0; i < num_threads; ++i) {
		t[i] = std::thread(call_from_thread, i);
	}

	for (int i = 0; i < num_threads; ++i) {
		t[i].join();
	}

	std::vector<double> neff_test = TRACSsim[0]->get_NeffParam();
	std::cout << "Neff param.: " << std::endl;
	for (int i = 0; i < 8; i++)
	{
		std::cout << neff_test[i] << std::endl;
	}

	TRACSsim[0]->write_to_file(0);

	for (int i = 0; i < TRACSsim.size(); i++)
	{
		delete TRACSsim[i];
	}


	printf("Hola\n");
	return 0;
}

