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

/************************************CARRIER COLLECTION***********************************
 *
 * Carrier collection get carriers information from a file and store it in a vector of carriers.
 * This vector of carriers will be used by each of the threads in each of the steps of the TCT scan.
 *
 *
 */

#include <memory>

#include "../include/CarrierCollection.h"

#include <TFile.h>
#include <TH1D.h>

double extra_y;

CarrierCollection::CarrierCollection(SMSDetector * detector) :
						_detector(detector)
{

}
/*
 * Adds carriers from the file. The name of the file must be set in the steering file, in the field CarrierFile, more instructions can be found there.
 * The calculation of average positions is neccesary for fitting purposes, DumpToTree function.
 */
/**
 *
 * @param filename
 */
void CarrierCollection::add_carriers_from_file(const std::string& filename, const std::string& scanType, double depth)
{
    std::vector<std::unique_ptr<Carrier>>_carrier_list_sngl;
        
	// get char representation and make ifstream
	//char * char_fn = filename.toLocal8Bit().data();
	std::ifstream infile(filename);
	bool once = true;
	char carrier_type;
	double q, x_init, y_init, gen_time;

	// process line by line
	std::string line;

    
	//Preprocessing for fitting bottom scan_type: Fixing mismatch between detector depth and y_init.
	//Outside the loop for performance purpose
	if (scanType == "bottom"){
		std::getline(infile, line);
		std::istringstream iss(line);
		if (!(iss >> carrier_type >> q >> x_init >> y_init >> gen_time)) {
			//std::cout << "Error while reading file" << std::endl;
		}
		//Extra in micrometers to shift y_init position
		extra_y = depth - y_init;

		//Calculate average beam position
		beamy += x_init;
		beamz += y_init + extra_y;

        //std::cout << "Charge is registered at x=" << x_init << ", y = " << y_init + extra_y << ", gen_time = " << gen_time << std::endl;
        std::unique_ptr<Carrier> ptr( new Carrier(carrier_type, q, x_init, y_init + extra_y, _detector, gen_time) );
        _carrier_list_sngl.push_back( std::move(ptr) );
        
		if ( _carrier_list_sngl.size()!=0 ) {
			beamy = beamy / _carrier_list_sngl.size();
			beamz = beamz / _carrier_list_sngl.size();
		}
	}
	while (std::getline(infile, line))
	{

		std::istringstream iss(line);
		if (!(iss >> carrier_type >> q >> x_init >> y_init >> gen_time)) { 
			//std::cout << "Error while reading file" << std::endl;
			break;
		} 

		//Calculate average beam position
		beamy += x_init;
		beamz += y_init + extra_y;

        std::unique_ptr<Carrier> ptr( new Carrier(carrier_type, q, x_init, y_init + extra_y, _detector, gen_time) );
        _carrier_list_sngl.push_back( std::move(ptr) );        
	}
	if ( _carrier_list_sngl.size()!=0 ) {
		beamy = beamy / _carrier_list_sngl.size();
		beamz = beamz / _carrier_list_sngl.size();
	}

    // Push the default carrier list to the queue
    _queue_carrier_list.push_back( std::move(_carrier_list_sngl) );
    
}

/*
 * From the carrier list taken from the file of carriers, this method runs through the vector of carriers and calls the simulate.drift for each one of them, depending if the carrier
 * is a hole or an electron. Trapping effects are included at the end directly in the valarry of currents.
 * Info related to diffusion is displayed using this method as well. Whenever a carrier crosses to depleted region, the acumulative variable totalCross grows.
 */
/**
 *
 * @param dt
 * @param max_time
 * @param shift_x
 * @param shift_y
 * @param curr_elec
 * @param curr_hole
 * @param totalCrosses
 *
 * Modified version !
 */
void CarrierCollection::simulate_drift( double dt, double max_time, double shift_x /*yPos*/, double shift_y /*zPos*/,
                                        std::valarray<double>&curr_elec, std::valarray<double> &curr_hole,
                                        std::valarray<double>&curr_gen_elec, std::valarray<double> &curr_gen_hole, double max_mul_factor,                               
                                        int &totalCrosses, const std::string &scanType, std::string skip_event_loop)
{
    double x_init, y_init;
    int totalCross = 0;
    bool control = true;

    //auto itr_queue = _queue_carrier_list.begin();
    //std::vector<std::unique_ptr<Carrier>> gen_carrier_list;    // container for carriers which will be generated by the impact ionization effect
    bool continue_loop = true;

    // Skip the following loop, if the corresponding flag is set
    if( skip_event_loop == "yes")
    {
        std::cout << std::endl;
        std::cout << "Event loop is going to be skipped ! " << std::endl;
        std::cout << "If you want to run normally, please change the setting (\"SkipEventLoop\") in config file " << std::endl;        
        std::cout << std::endl;
        
        continue_loop = false;
    }
    int loop = 0;
    
    //auto initial_n_carrier = (*itr_queue).size();
    auto generated_n_carrier = 0;
    
    while( continue_loop )
    {        
        std::vector<std::unique_ptr<Carrier>> gen_carrier_list;    // container for carriers which will be generated by the impact ionization effect

        // Declare the iterator and increment it here to avoid the invalidation of the iterator.
        // Probably, it is fine, as long as std::deque with "puch_back" insertion method, thus it is kinds of safety reason.
        auto itr_queue = _queue_carrier_list.begin();
        auto initial_n_carrier = (*itr_queue).size();
        for( int i=0; i<loop; i++)
        {
            ++itr_queue;
            //std::cout << "Loop=" << loop << " : Size of std::vector<std::unique_ptr<Carrier> =  " << (*itr_queue).size() << std::endl;            
            //std::cout << "Loop=" << loop << " : Size of std::deque<std::vector<std::unique_ptr<Carrier>>> =  " << _queue_carrier_list.size() << std::endl;
        }
        
        for (auto& carrier : *itr_queue )   // modified            
        {               
            char carrier_type = carrier->get_carrier_type();  
            // simulate drift and add to proper valarray
            if (carrier_type == 'e')
            {                
                // get and shift carrier position
                std::array<double,2> x = carrier->get_x();     
                x_init = x[0] + shift_x;
                y_init = x[1] + shift_y;

                //std::cout << "x[0]=" << x[0] << std::endl;
                
                //x[0] represents the X position read from the carriers file
                //shift_x represents the shift applied to X read from the steering file, namely, where the laser points to.
                //x[1] represents the Y position read from the carriers file. Y is seen as Z.
                //shift_y represents the Z (y) position read from the steering file. Since the program (usually) does edge-TCT, Z can be defined in one or more steps.                

                if( loop == 0 ) // Current made by initial carriers
                {
                    curr_elec += carrier->simulate_drift( dt , max_time, x_init, y_init, scanType,  gen_carrier_list); 
                }
                else                //  Current made by secondary carriers
                {
                    curr_gen_elec += carrier->simulate_drift( dt , max_time, x_init, y_init, scanType,  gen_carrier_list); 
                }
            }            
            else if (carrier_type =='h')
            {
                // get and shift carrier position
                std::array<double,2> x = carrier->get_x();      
                double x_init = x[0] + shift_x;
                double y_init = x[1] + shift_y ;               

                if( loop == 0 ) // Current made by initial carriers
                {
                    curr_hole += carrier->simulate_drift( dt , max_time, x_init, y_init, scanType,  gen_carrier_list);   
                }
                else
                {
                    curr_gen_hole += carrier->simulate_drift( dt , max_time, x_init, y_init, scanType,  gen_carrier_list);    
                }
            }
            //Let's see how many carriers from the carrier list for this step in Z have crossed to the depleted region.
            //A flag is switched to true on carrier.simulate_drift when a carrier filfull the requirements. See carrier.simulate_drift
            if (carrier->crossed()){    
                totalCross += 1;
            }
            
        } // end of for( auto carrier : *itr_queue )

        
        // push back the newly generated carrier list to the Queue
        if( !gen_carrier_list.empty() )
        {
            generated_n_carrier += gen_carrier_list.size();
            std::cout << "Initial number of carriers = " << initial_n_carrier << std::endl;
            std::cout << "Number of carriers generated in this loop = " << gen_carrier_list.size() << std::endl;            
            std::cout << "Total number of generated carriers = " << generated_n_carrier << std::endl;
            std::cout << "loop number = " << loop << std::endl;
            std::cout << std::endl;
            
            if( generated_n_carrier > initial_n_carrier * max_mul_factor )  // If number of carriers exceeds a thredhols, stop the loop. (to avoid endless loop)
            {
                std::cout << "Number of generated charges exceeds the maxmimum limit set in the configuration file !" << std::endl;
                std::cout << "Therefore, finishing the iteration of carrier drift loop. The results obtained by now is successfully recorded. " << std::endl;
                continue_loop = false;
            }
            else
            {
                continue_loop = true;
                _queue_carrier_list.push_back( std::move(gen_carrier_list) );                
                loop++;
            }
        }
        else
            continue_loop = false;
            
    }   // end of while loop 

    
	//fileDiffDrift.close();
	//std::cout << "Number of carriers crossed to DR in last Z step with Height " << shift_y << ": " << totalCross << std::endl;
	totalCrosses += totalCross;

	double trapping_time = _detector->get_trapping_time();

	for (double i = 0.; i < curr_hole.size(); i ++)
	{
		double elapsedT = i*dt;
		curr_elec[i] *= exp(-elapsedT/trapping_time);
		curr_hole[i] *= exp(-elapsedT/trapping_time);

        curr_gen_elec[i] *= exp(-elapsedT/trapping_time);
		curr_gen_hole[i] *= exp(-elapsedT/trapping_time);
	}

    // Record number of carriers vs gen_time
    record_carrier_gen_time( max_time, curr_hole.size() );
}

/**
 * Record number of carriers
 **/
void CarrierCollection::record_carrier_gen_time(double max_time, int n_time_slice)
{
    TFile *f = new TFile("ncarrier.root", "RECREATE");

    TH1D *h_e_gen_time = new TH1D("e_gentime", "electron gen_time", n_time_slice, 0.0, max_time);
    TH1D *h_h_gen_time = new TH1D("h_gentime", "hole gen_time", n_time_slice, 0.0, max_time);    
   
    for( auto itr = _queue_carrier_list.begin();  itr != _queue_carrier_list.end(); itr++)
    {
        for (auto& carrier : *itr ) 
        {
            char carrier_type = carrier->get_carrier_type();
            double gen_time = carrier->get_gen_time();

            if( carrier_type == 'e') 
                h_e_gen_time->Fill( gen_time );
            else
                h_h_gen_time->Fill( gen_time );
        }
    }

    h_e_gen_time->Write();
    h_h_gen_time->Write();

    f->Close();  
}



/**
 *
 * @param n_bins_x
 * @param n_bins_y
 * @param hist_name
 * @param hist_title
 * @return
 */
TH2D CarrierCollection::get_e_dist_histogram(int n_bins_x, int n_bins_y,  TString hist_name, TString hist_title)
{
	// get detector limits
	double x_min = _detector->get_x_min();
	double x_max = _detector->get_x_max();
	double y_min = _detector->get_y_min();
	double y_max = _detector->get_y_max();

	// create histogram object
	TH2D e_dist = TH2D(hist_name, hist_title, n_bins_x , x_min, x_max, n_bins_y, y_min, y_max);

    auto itr_queue = _queue_carrier_list.begin();
	for (auto& carrier : *itr_queue)   
	{
		char carrier_type = carrier->get_carrier_type();
		if (carrier_type == 'e')
		{
			std::array< double,2> x = carrier->get_x();
			double q = carrier->get_q();
			e_dist.Fill(x[0], x[1], q);
		}
	}
	return e_dist;
}
/**
 *
 * @param n_bins_x
 * @param n_bins_y
 * @param shift_x
 * @param shift_y
 * @param hist_name
 * @param hist_title
 * @return
 */
TH2D CarrierCollection::get_e_dist_histogram(int n_bins_x, int n_bins_y, double shift_x, double shift_y, TString hist_name, TString hist_title)
{
	// get detector limits
	double x_min = _detector->get_x_min();
	double x_max = _detector->get_x_max();
	double y_min = _detector->get_y_min();
	double y_max = _detector->get_y_max();

	// create histogram object
	TH2D e_dist = TH2D(hist_name, hist_title, n_bins_x , x_min, x_max, n_bins_y, y_min, y_max);

    auto itr_queue = _queue_carrier_list.begin();
	for (auto& carrier : *itr_queue)   
	{
		char carrier_type = carrier->get_carrier_type();
		if (carrier_type == 'e')
		{
			std::array< double,2> x = carrier->get_x();
			double q = carrier->get_q();
			e_dist.Fill(x[0]+shift_x, x[1]+shift_y, q);
		}
	}
	return e_dist;
}

/*
 ********************** DESTRUCTOR OF THE CLASS CARRIER	COLLECTION **************************
 */
CarrierCollection::~CarrierCollection()
{

}
