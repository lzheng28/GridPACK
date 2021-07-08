/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
/**
 * @file   hw_app.cpp
 * @author Bruce Palmer
 * @date   2014-01-28 10:31:00 d3g096
 * 
 * @brief  
 * 
 * 
 */
// -------------------------------------------------------------

#include <iostream>
#include "boost/smart_ptr/shared_ptr.hpp"
#include "gridpack/include/gridpack.hpp"
#include "hw_app.hpp"
#include "hw_factory.hpp"

// #include "json.hpp"

#ifdef USE_HELICS
// #include "helics/ValueFederates.hpp"
// #include <helics/shared_api_library/ValueFederate.h>
#include <helics/helics.hpp>
#endif


// Calling program for hello world application

/**
 * Basic constructor
 */
gridpack::hello_world::HWApp::HWApp(void)
{
}

/**
 * Basic destructor
 */
gridpack::hello_world::HWApp::~HWApp(void)
{
}

/**
 * Execute application
 * @param argc number of arguments
 * @param argv list of character strings
 */
void gridpack::hello_world::HWApp::execute(int argc, char** argv)
{
  // load input file
  gridpack::parallel::Communicator world;
  boost::shared_ptr<HWNetwork> network(new HWNetwork(world));

  // read configuration file
  std::string filename = "10x10.raw";

  // Read in external PTI file with network configuration
  gridpack::parser::PTI23_parser<HWNetwork> parser(network);
  parser.parse(filename.c_str());

  int nTotalBus = network->totalBuses();

  std::cout << "number of total buses: " << nTotalBus << std::endl;

  // partition network
  network->partition();

  // create factory
  gridpack::hello_world::HWFactory factory(network);
  factory.load();

  // leizheng debug:
  // std::cout << "network" << network->totalBuses() << std::endl;
  int nBus = network->numBuses();
  std::cout << "nBus: " << nBus << std::endl;

  int outnodeIndex = 6; //use node 0 to transfer data

  int myid;
  MPI_Comm_rank(MPI_COMM_WORLD, &myid);

  std::cout << "myid:" << myid << std::endl; 

#ifdef USE_HELICS
  std::string configFile = "/home/lei/Desktop/test-helics/federates/helics_config_1.json";
  // helics::ValueFederate fed(configFile);
  std::shared_ptr<helics::ValueFederate> fed;
  helics::Publication pub;
  helics::Input sub;
  double helics_requestTime;
  double helics_grantime;
  
  //to get publication definitions
  int pubCount;
  int subCount;

#endif //end if of HELICS

#ifdef USE_HELICS

  for(int i = 0; i < nBus; i++){
    // std::cout << "Global Index:"<<network->getGlobalBusIndex(i) << std::endl;
    if(network->getGlobalBusIndex(i) == outnodeIndex && network->getActiveBus(i)){
      //do helics here.
      fed = std::make_shared<helics::ValueFederate>(configFile);
      std::cout << "*****Do helics****" <<  outnodeIndex << std::endl;
      
      //std::cout << "-------------!!!helics test: HELICS Version: " << helics::versionString << std::endl;
      std::cout << "-------------!!!helics test: HELICS Version: " << helics::versionString << std::endl;
      std::string configFile = "/home/lei/Desktop/test-helics/federates/helics_config_1.json";
      // helics::ValueFederate fed(configFile); // int a = new int();
      // fed = helics::ValueFederate(configFile);
      helics_requestTime = 0.0;

      // helics::Publication pub;
      // helics::Input sub;
      
      //to get publication definitions
      pubCount = (*fed).getPublicationCount();
      
      printf("-------------helics test: num of pub: %d \n", pubCount);
      for(int i = 0; i < pubCount; i++) {
          pub = (*fed).getPublication(i);
          std::string pubInfo = pub.getInfo();
          // do stuff to tie pub to GridPACK object property
      }
        
      //to get subscription definitions
      subCount = (*fed).getInputCount();
      printf("-------------helics test: num of sub: %d \n", subCount);
      
      for(int i = 0; i < subCount; i++) {
          sub = (*fed).getInput(i);
          // std::string subInfo = sub.getInfo(); printf("----------------leizheng debug, subInfo: %s", subInfo);
          // do stuff to tie pub to GridPACK object property
      }
      //let helics broker know you are ready to start simulation 
      (*fed).enterExecutingMode();
#endif  //end if of HELICS
    }
  }


  //Process

  for(int k = 4; k < 6; k++){
    for(int i = 0; i < nBus; i++){ // for each processor, do
      // std::cout << "Global Index:"<<network->getGlobalBusIndex(i) << std::endl;
      if(network->getGlobalBusIndex(i) == outnodeIndex && network->getActiveBus(i)){
        //do helics here.
        std::cout << "*****Do helics receive&send ****" <<  outnodeIndex << std::endl;

  #ifdef USE_HELICS
    for(int i = 0; i < pubCount; i++) {
          pub = (*fed).getPublication(i);
          std::string pubInfo = pub.getInfo();
          // double pub_info = 1.0;
          // auto pub_info = helics::helicsGetComplex(helics::helicsComplexString(1, 1));
          std::complex<double> pub_info {k, k};
          // pub_info.real = 1, pub_info.imag = 1;
          pub.publish(pub_info);
      }

      // helics_requestTime = 0.05;
      // double helics_grantime;
      helics_requestTime += 0.005;
      helics_grantime = (*fed).requestTime(helics_requestTime);

      // double subvalue = 0.0;
      std::complex<double> subvalue {0, 0};
      for(int i = 0; i < subCount; i++) {
          sub = (*fed).getInput(i);
          printf("-------------!!!helics debug entering  sub loop\n"); 
          // if(sub.isUpdated()) {
          // auto subvalue = fed.getDouble(sub);
          sub.getValue(subvalue);
          // subvalue = fed.getDouble(sub);
          // printf("-------------!!!Helics sub value: %s \n", subvalue);
          std::cout << "fed 1 subvalue:  " << subvalue << std::endl;
                              //update GridPACK object property with value
          // }
      }
  #endif  //end if of HELICS

      }
    }
  }
// #ifdef USE_HELICS
// 	fed.finalize();
// #endif

  // BusPtr zero_bus = network->getBus(0);

  // Create serial IO object to export data from buses
  gridpack::serial_io::SerialBusIO<HWNetwork> busIO(128,network);
  // busIO.header("\nMessage from buses\n");
  // busIO.write();


  // Create serial IO object to export data from branches
  gridpack::serial_io::SerialBranchIO<HWNetwork> branchIO(128,network);
  // branchIO.header("\nMessage from branches\n");
  // branchIO.write();
}
