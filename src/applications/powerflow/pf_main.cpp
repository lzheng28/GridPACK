
/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
/**
 * @file   pf_main.cpp
 * @author Bruce Palmer
 * @date   2016-07-14 14:23:07 d3g096
 *
 * @brief
 */
// -------------------------------------------------------------

#include "mpi.h"
#include <ga.h>
#include <macdecls.h>
#include "gridpack/include/gridpack.hpp"
#include "gridpack/applications/modules/powerflow/pf_app_module.hpp"

#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <sstream>

#ifdef USE_HELICS
//#include "helics/ValueFederates.hpp"
//#include <helics/shared_api_library/ValueFederate.h>
#include <helics/helics.hpp>
#endif

const char* help = "GridPACK power flow application";

int main(int argc, char **argv)
{
  // Initialize libraries (parallel and math)
  gridpack::Environment env(argc,argv,help);

  if (1) {
    gridpack::utility::CoarseTimer *timer =
      gridpack::utility::CoarseTimer::instance();
    gridpack::parallel::Communicator world;

    // read configuration file
    gridpack::utility::Configuration *config =
      gridpack::utility::Configuration::configuration();
    if (argc >= 2 && argv[1] != NULL) {
      char inputfile[256];
      sprintf(inputfile,"%s",argv[1]);
      config->open(inputfile,world);
    } else {
      config->open("input.xml",world);
    }

    gridpack::utility::Configuration::CursorPtr cursor;
    cursor = config->getCursor("Configuration.Powerflow");
    bool useNonLinear = false;
    useNonLinear = cursor->get("UseNonLinear", useNonLinear);
    bool exportPSSE23 = false;
    std::string filename23;
    exportPSSE23 = cursor->get("exportPSSE_v23",&filename23);
    bool exportPSSE33 = false;
    std::string filename33;
    exportPSSE33 = cursor->get("exportPSSE_v33",&filename33);
    bool noPrint = false;
    cursor->get("suppressOutput",&noPrint);

    // setup and run powerflow calculation
    boost::shared_ptr<gridpack::powerflow::PFNetwork>
      pf_network(new gridpack::powerflow::PFNetwork(world));

    gridpack::powerflow::PFAppModule pf_app;
    if (noPrint) {
      pf_app.suppressOutput(noPrint);
    }
    pf_app.readNetwork(pf_network,config);
    pf_app.initialize();


// #ifdef USE_HELICS
#if 0
    //std::cout << "-------------!!!helics test: HELICS Version: " << helics::versionString << std::endl;
    std::shared_ptr<helics::ValueFederate> fed;
    helics::Publication pub;
    helics::Input sub;
    double helics_requestTime;
    double helics_grantime;
    std::complex<double> subvalue = {1.0, 1.0};
    
    //to get publication definitions
    int pubCount;
    int subCount;

    int nBus = pf_network->numBuses();
    int outnodeIndex = 118 - 1; //use node 118 to transfer data

    std::cout << "nBus:              " << nBus << std::endl;

    for(int i = 0; i < nBus; i++){
      if(pf_network->getGlobalBusIndex(i) == outnodeIndex && pf_network->getActiveBus(i)){

        std::cout << "-------------!!!helics test: HELICS Version: " << helics::versionString << std::endl;
        std::string configFile = "/home/lei/Desktop/test-helics/federates/helics_config_1.json";
        // string configFile = "/home/huan495/gridpack-dev/src/build/applications/dynamic_simulation_full_y/testcase/helics_39bus_3.json";
        //   helics::ValueFederate fed(configFile);
        // helics::Publication pub;
        // helics::Input sub;
        // double helics_requestTime = 0.0;
        
        //to get publication definitions
        fed = std::make_shared<helics::ValueFederate>(configFile);
        pubCount = (*fed).getPublicationCount();
        
        printf("-------------helics test: num of pub: %d \n", pubCount);
        for(int j = 0; j < pubCount; j++) {
          pub = (*fed).getPublication(j);
          std::string pubInfo = pub.getInfo();
          // do stuff to tie pub to GridPACK object property
        }
        
        //to get subscription definitions
        subCount = (*fed).getInputCount();
        printf("-------------helics test: num of sub: %d \n", subCount);
        
        for(int j = 0; j < subCount; j++) {
          sub = (*fed).getInput(j);
          std::string subInfo = sub.getInfo();
          // do stuff to tie pub to GridPACK object property
        }

        //let helics broker know you are ready to start simulation 
        (*fed).enterExecutingMode();	
      }
    }
#endif  //end if of HELICS
    // double temp = 10.0;
    int simu_k = 288; //288
    double h_sol1 = 300;

    std::complex<double> pub_info = {2.0, 2.0};

#if 1
    // read file
    std::ifstream inFile("power_load.csv", std::ios::in);
    std::string lineStr;
    std::vector<std::vector<double>> strArray;
    while (getline(inFile, lineStr))
    {
      std::stringstream ss(lineStr);
      std::string str;
      std::vector<double> lineArray;
      while (getline(ss, str, ','))
          lineArray.push_back(atof(str.c_str()));
      strArray.push_back(lineArray);
    }
#endif
    std::vector<std::complex<double>> v_voltage_cosim;

    for(int I_Steps = 0; I_Steps < simu_k; I_Steps++){

      // Set actual demand
#if 1
      for(int k = 0; k < pf_network->numBuses(); k++){
        gridpack::powerflow::PFBus *bus = dynamic_cast<gridpack::powerflow::PFBus*>(pf_network->getBus(k).get());

        if((bus->p_pl).size() > 0){
          // bus->p_pl[0] += 0.001;
          // bus->p_ql[0] += 0.001;
          bus->p_pl[0] = strArray[I_Steps][pf_network->getGlobalBusIndex(k)];
          bus->p_ql[0] = strArray[I_Steps][pf_network->getGlobalBusIndex(k)] * tan(acos(0.85));
          // std::cout << "&&&&& bus->p_pl[0] = " << bus->p_pl[0] << "  &&&&&bus->p_ql[0] = " << bus->p_ql[0] << std::endl;
        }
      }
#endif
      //collect voltage
      std::complex<double> voltage_cosim = {0, 0};
#if 0
      for(int k = 0; k < pf_network->numBuses(); k++){
        if(pf_network->getGlobalBusIndex(k) == 117 && pf_network->getActiveBus(k)){
          gridpack::powerflow::PFBus *bus = dynamic_cast<gridpack::powerflow::PFBus*>(pf_network->getBus(k).get());
          voltage_cosim = {bus->getVoltage() * 138.0 * 1000, 0};
          // std::complex<double> voltage_cosim{142300, 0};
          // std::cout << "$$$$$voltage_cosim = " << voltage_cosim << std::endl;
          // std::complex<double> pub_info = {};
          if(I_Steps == 0){
            voltage_cosim *= 1.043;
          }
          v_voltage_cosim.push_back(voltage_cosim);
        }
      }
#endif
// #ifdef USE_HELICS
#if 0
      for(int i = 0; i < nBus; i++){
        if(pf_network->getGlobalBusIndex(i) == outnodeIndex && pf_network->getActiveBus(i)){
          //pub.publish(widearea_deltafreq);
          std::cout << "Receive&send helics info" << std::endl;
          for(int j = 0; j < pubCount; j++) {
            pub = (*fed).getPublication(j);
            std::string pubInfo = pub.getInfo();
            //std::cout << "-------------!!!helics test: HELICS pub info: " << pubInfo << std::endl;
            // pub.publish(widearea_deltafreq);
            gridpack::powerflow::PFBus *bus = dynamic_cast<gridpack::powerflow::PFBus*>(pf_network->getBus(i).get());
            std::complex<double> voltage_cosim {bus->getVoltage() * 138.0 * 1000, 0};
            // std::complex<double> voltage_cosim{142300, 0};
            std::cout << "$$$$$voltage_cosim = " << voltage_cosim << std::endl;
            // std::complex<double> pub_info = {};
            if(I_Steps == 0){
              voltage_cosim *= 1.043;
            }
            pub.publish(voltage_cosim);
            // do stuff to tie pub to GridPACK object property
          }

          helics_requestTime = double (I_Steps*h_sol1);
          printf("-------------!!!Helics request time: %12.6f \n", helics_requestTime);
          std::cout << "I_Steps: " << I_Steps << "h_sol1: " << h_sol1 << std::endl; 
          //  double helics_grantime;
          helics_grantime = (*fed).requestTime(helics_requestTime);
          //  printf("-------------!!!Helics grant time: %12.6f \n", helics_grantime); 
          
          // subvalue = 1.0;
          
          for(int j = 0; j < subCount; j++) {
            sub = (*fed).getInput(j);
            //printf("-------------!!!helics debug entering  sub loop\n"); 
            //if(sub.isUpdated()) {
            sub.getValue(subvalue);
            // subvalue = (*fed).getDouble(sub);
            std::cout << "subvalue:" << subvalue << std::endl; std::cout << "line:340" << " world rank: " << std::endl;
              // printf("-------------!!!Helics sub value: %12.6f \n", subvalue);
                                    //update GridPACK object property with value
                //}
          }
          //printf("-------------!!!Outside Helics def sub value: %12.6f \n", subvalue);

          int load_amplification_factor = 15;

          gridpack::powerflow::PFBus *bus = dynamic_cast<gridpack::powerflow::PFBus*>(pf_network->getBus(i).get());
          bus->p_pl[0] = subvalue.real() * load_amplification_factor / 1000000;
          bus->p_ql[0] = subvalue.imag() * load_amplification_factor / 1000000;
        }
      }
#endif  //end if of HELICS

#if 0
      for(int k = 0; k < pf_network->numBuses(); k++){
        gridpack::powerflow::PFBus *bus = dynamic_cast<gridpack::powerflow::PFBus*>(pf_network->getBus(k).get());

        if((bus->p_pl).size() > 0){
          bus->p_pl[0] = strArray[I_Steps][k];
          bus->p_ql[0] = strArray[I_Steps][k] * tan(acos(0.85));
          std::cout << "&&&&& bus->p_pl[0] = " << bus->p_pl[0] << "  &&&&&bus->p_ql[0] = " << bus->p_ql[0] << std::endl;
        }
      }
#endif

      if (useNonLinear) {
        pf_app.nl_solve();
      } else {
        pf_app.solve();
      }
      //pf_app.write();
    }

#if 0
      for(int k = 0; k < pf_network->numBuses(); k++){
        if(pf_network->getGlobalBusIndex(k) == 117 && pf_network->getActiveBus(k)){
          std::ofstream file("voltage_gp.csv");
          if (file)
          {
            file << "voltage" << "\n";
            for(int i = 0; i < v_voltage_cosim.size(); i++){
              if(i != 0){
                for(int j = 0; j < 5; j++){
                  file << v_voltage_cosim[i].real() << "+" << v_voltage_cosim[i].imag() << "j" << "\n";
                }
              } else {
                file << v_voltage_cosim[i].real() << "+" << v_voltage_cosim[i].imag() << "j" << "\n";
              }
            }
          }
          for(int i = 0; i < 4; i++){
            file << v_voltage_cosim[v_voltage_cosim.size() - 1].real() << "+" << v_voltage_cosim[v_voltage_cosim.size() - 1].imag() << "j" << "\n";
          }
          file.close();
        }
      }
#endif
    pf_app.write();
    pf_app.saveData();
    if (exportPSSE23) {
      pf_app.exportPSSE23(filename23);
    }
    if (exportPSSE33) {
      pf_app.exportPSSE33(filename33);
    }
    if (!noPrint) {
      timer ->dump();
    }
  }

  return 0;
}

