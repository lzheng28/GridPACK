
/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
/**
 * @file   hadrec_main.cpp
 * @author Bruce Palmer
 * @date   2020-04-23 13:26:55 d3g096
 *
 * @brief
 */
// -------------------------------------------------------------

#include "mpi.h"
#include <ga.h>
#include <macdecls.h>
#include "gridpack/include/gridpack.hpp"
#include "gridpack/applications/modules/hadrec/hadrec_app_module.hpp"
#include "gridpack/applications/modules/dynamic_simulation_full_y/dsf_app_module.hpp"
#include <vector>
#include <string>

#define ARRAY_LEN 2

#ifdef USE_HELICS
//#include "helics/ValueFederates.hpp"
//#include <helics/shared_api_library/ValueFederate.h>
#include <helics/helics.hpp>
#endif

const char* help = "HADREC Test application";

int main(int argc, char **argv)
{
  // Initialize libraries (parallel and math)
  gridpack::NoPrint *noprint_ins = gridpack::NoPrint::instance();
  //  noprint_ins->setStatus(true);

  //bool bnoprint = gridpack::NoPrint::instance()->status();
  //printf ("------------- hadrec_main function test 1  bnoprint: %d \n", bnoprint);

  gridpack::Environment env(argc,argv,help);
  
  gridpack::parallel::Communicator world;
  int me = world.rank();

  //bnoprint = gridpack::NoPrint::instance()->status();
  //printf ("------------- hadrec_main function test 2  bnoprint: %d \n", bnoprint);

  gridpack::utility::CoarseTimer *timer =
    gridpack::utility::CoarseTimer::instance();
  int t_total = timer->createCategory("Dynamic Simulation: Total Application");

  //allocate the memory for the HADRECAppModule
  boost::shared_ptr<gridpack::hadrec::HADRECAppModule>
    hadrec_app_sptr (new gridpack::hadrec::HADRECAppModule() );

  // solve power flow
  std::string file;
  if (argc > 1) {
    file = argv[1];
  } else {
    file = "input.xml";
  }

  //bnoprint = gridpack::NoPrint::instance()->status();
  //printf ("------------- hadrec_main function test 3  bnoprint: %d \n", bnoprint);

  hadrec_app_sptr->solvePowerFlowBeforeDynSimu(const_cast<char *>(file.c_str()));

  //hadrec_app_sptr->solvePowerFlowBeforeDynSimu(const_cast<char *>(file.c_str()));

  // transfer power flow results to dynamic simulation
  hadrec_app_sptr->transferPFtoDS();

  std::vector<gridpack::dynamic_simulation::Event> BusFaults;
  BusFaults.clear();
  
  //printf ("------------- hadrec_main function before initializeDynSimu \n");

  // initialize dynamic simulation
  hadrec_app_sptr->initializeDynSimu(BusFaults);

  int outnodeIndex = 1;

  boost::shared_ptr<gridpack::dynamic_simulation::DSFullNetwork> network = hadrec_app_sptr->ds_network;

  double time_step = hadrec_app_sptr->getTimeStep();

  bool use_helics = false;

  use_helics = hadrec_app_sptr->useHelicsStatus();

// #ifdef USE_HELICS
// #if 0
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

  int nBus = network->numBuses();
  // int outnodeIndex = 2 - 1; //use node 118 to transfer data

  std::cout << "nBus:              " << nBus << std::endl;
  // std::vector<int> localIndices;
  // localIndices = network->getLocalBusIndices(504);

  int connectedBusID = hadrec_app_sptr->getHelicsConnectNode();
  std::vector<int> localIndices = network->getLocalBusIndices(connectedBusID);
  std::cout << "connectedBusID: " << connectedBusID << ", me" << std::endl;

  for(int i = 0; i < localIndices.size(); i++){
    if(network->getActiveBus(localIndices[i])){
      std::cout << "Active Local indics, me: "<< localIndices[i] << " " << me << std::endl;
    } else {
      std::cout << "Inactivate Local indics, me: " << localIndices[i] << " " << me << std::endl;
    }
  }

  MPI_Win mpi_win;
  MPI_Aint mpi_size;
  double *baseptr;
  if (me == 0)
   {
      mpi_size = ARRAY_LEN * sizeof(double);
      MPI_Win_allocate_shared(mpi_size, sizeof(double), MPI_INFO_NULL,
                              MPI_COMM_WORLD, &baseptr, &mpi_win);
   }
   else
   {
      int disp_unit;
      MPI_Win_allocate_shared(0, sizeof(double), MPI_INFO_NULL,
                              MPI_COMM_WORLD, &baseptr, &mpi_win);
      MPI_Win_shared_query(mpi_win, 0, &mpi_size, &disp_unit, &baseptr);
   }
   double *arr = baseptr;

  if(use_helics){
    if(localIndices.size() > 0 && network->getActiveBus(localIndices[0])){
      outnodeIndex = localIndices[0];
      std::cout << "-------------!!!helics test: HELICS Version: " << helics::versionString << std::endl;
      std::string configFile = hadrec_app_sptr->getHelicsConfigFile();
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
// #endif  //end if of HELICS

  //printf ("------------- hadrec_main function after initializeDynSimu \n");

  bool debugoutput = false; // whether print out debug staffs
  bool boutputob = false;
  double lp, lq, pg, qg;
  int busno = 5;
  bool btmp;

  //-----test get load and get generator function-----------
if (debugoutput){

/*
  busno = 5;
  btmp = hadrec_app_sptr->getBusTotalLoadPower(busno, lp, lq);
  if (me==0) printf("------------test hadrec_main, load at bus %d, has P: %f, Q: %f\n", busno, lp, lq);
  busno = 7;
  btmp = hadrec_app_sptr->getBusTotalLoadPower(busno, lp, lq);
  if (me==0) printf("------------test hadrec_main, load at bus %d, has P: %f, Q: %f\n", busno, lp, lq);
  busno = 9;
  btmp = hadrec_app_sptr->getBusTotalLoadPower(busno, lp, lq);
  if (me==0) printf("------------test hadrec_main, load at bus %d, has P: %f, Q: %f\n", busno, lp, lq);
*/

  busno = 36;
  std::string genid = "1 ";
  btmp = hadrec_app_sptr->getGeneratorPower(busno, genid, pg, qg);
  if (me==0) printf("------------test hadrec_main, %d find generator at bus %d, has P: %f, Q: %f\n", btmp, busno, pg, qg);

  busno = 36;
  genid = "1";
  btmp = hadrec_app_sptr->getGeneratorPower(busno, genid, pg, qg);
  if (me==0) printf("------------test hadrec_main, %d find generator at bus %d, has P: %f, Q: %f\n", btmp, busno, pg, qg);

}


  gridpack::hadrec::HADRECAction loadshedact;
  loadshedact.actiontype = 0;
  loadshedact.bus_number = 5;
  loadshedact.componentID = "1";
  loadshedact.percentage = -0.2;

  gridpack::hadrec::HADRECAction loadshedact1;
  loadshedact1.actiontype = 0;
  loadshedact1.bus_number = 7;
  loadshedact1.componentID = "1";
  loadshedact1.percentage = -0.2;

  gridpack::hadrec::HADRECAction linetrip;
  linetrip.actiontype = 1;
  //linetrip.bus_number = 6;
  linetrip.brch_from_bus_number = 6;
  linetrip.brch_to_bus_number = 7;
  linetrip.branch_ckt = "1 ";
  
  gridpack::hadrec::HADRECAction loadpchange;
  loadpchange.actiontype = 3;
  loadpchange.bus_number = 501;
  loadpchange.percentage = 200.0;


  int isteps = 0;
  bool bApplyAct_LoadShedding = false;  // whether apply the load shedding action in the simulation steps
  bool bApplyAct_LineTripping = false;  // whether apply the line tripping action in the simulation steps
  bool bApplyAct_LoadPchange = false;  // apply a sudden load P change, could mimic fault
  std::vector<double> ob_vals;
  int idxtmp;

  // test getOblist
  std::vector<int> obs_genBus;
  std::vector<std::string> obs_genIDs;
  std::vector<int> obs_loadBus;
  std::vector<std::string> obs_loadIDs;
  std::vector<int> obs_vBus;
  hadrec_app_sptr->getObservationLists(obs_genBus, obs_genIDs,
      obs_loadBus, obs_loadIDs, obs_vBus);

  if (boutputob && me == 0){
    printf("-----------renke debug, getObservationLists------------\n");
    printf("-----------ob gen bus list, ");
    for (idxtmp=0; idxtmp<obs_genBus.size(); idxtmp++){
      printf(" %d, ", obs_genBus[idxtmp]);
    }
    printf(" \n");

    printf("-----------ob gen ID list, ");
    for (idxtmp=0; idxtmp<obs_genIDs.size(); idxtmp++){
      printf(" %s, ", obs_genIDs[idxtmp].c_str());
    }
    printf(" \n");
	
	printf("-----------ob bus list, ");
    for (idxtmp=0; idxtmp<obs_vBus.size(); idxtmp++){
      printf(" %d, ", obs_vBus[idxtmp]);
    }
	printf(" \n");

    printf("-----------ob load bus list, ");
    for (idxtmp=0; idxtmp<obs_loadBus.size(); idxtmp++){
      printf(" %d, ", obs_loadBus[idxtmp]);
    }
    printf(" \n");

    printf("-----------ob load ID list, ");
    for (idxtmp=0; idxtmp<obs_loadIDs.size(); idxtmp++){
      printf(" %s, ", obs_loadIDs[idxtmp].c_str());
    }

    printf(" \n");
  }

  std::vector<int> zone_id;
  std::vector<double> tmp_p;
  std::vector<double> tmp_q;

  if (debugoutput && me == 0){
	hadrec_app_sptr->getZoneLoads(tmp_p, tmp_q, zone_id);
  
    printf("\n-------------------get zone load information, total zones: %d \n\n", zone_id.size());
    for (idxtmp=0; idxtmp<zone_id.size(); idxtmp++){

      printf(" zone number: %d, total load p: %f, total load q: %f,\n", zone_id[idxtmp], tmp_p[idxtmp], tmp_q[idxtmp]);

    }
  }
  
  if (debugoutput && me == 0){
	hadrec_app_sptr->getZoneGeneratorPower(tmp_p, tmp_q, zone_id);
  
    printf("\n-------------------get zone generation information, total zones: %d \n\n", zone_id.size());
    for (idxtmp=0; idxtmp<zone_id.size(); idxtmp++){

      printf(" zone number: %d, total generation p: %f, total generation q: %f,\n", zone_id[idxtmp], tmp_p[idxtmp], tmp_q[idxtmp]);

    }
  }

  if (boutputob && me == 0) {
    printf("observations,  ");
	for (idxtmp=0; idxtmp<obs_genBus.size(); idxtmp++){
      printf("gen-at-bus-%d-ID-%s-speed,  ", obs_genBus[idxtmp], obs_genIDs[idxtmp].c_str());
	 
    }
	
	for (idxtmp=0; idxtmp<obs_genBus.size(); idxtmp++){
      //printf("gen-at-bus-%d-ID-%s-speed,  ", obs_genBus[idxtmp], obs_genIDs[idxtmp].c_str());
	  printf("gen-at-bus-%d-ID-%s-angle,  ", obs_genBus[idxtmp], obs_genIDs[idxtmp].c_str());

    }
	
	for (idxtmp=0; idxtmp<obs_genBus.size(); idxtmp++){

	  printf("gen-at-bus-%d-ID-%s-P,  ", obs_genBus[idxtmp], obs_genIDs[idxtmp].c_str());

    }
	
	for (idxtmp=0; idxtmp<obs_genBus.size(); idxtmp++){

	  printf("gen-at-bus-%d-ID-%s-Q,  ", obs_genBus[idxtmp], obs_genIDs[idxtmp].c_str());
    }
	
    for (idxtmp=0; idxtmp<obs_vBus.size(); idxtmp++){
      printf("bus-%d-magnitude, ", obs_vBus[idxtmp]);
	  //printf("bus-%d-angle, ", obs_vBus[idxtmp]);
    }
	
	for (idxtmp=0; idxtmp<obs_vBus.size(); idxtmp++){
      //printf("bus-%d-magnitude, ", obs_vBus[idxtmp]);
	  printf("bus-%d-angle, ", obs_vBus[idxtmp]);
    }

    for (idxtmp=0; idxtmp<obs_loadBus.size(); idxtmp++){
      printf("remaining-load-at-bus-%d-ID-%s,  ", obs_loadBus[idxtmp], obs_loadIDs[idxtmp].c_str());
    }
    printf(" \n");

  }

  while(!hadrec_app_sptr->isDynSimuDone()){
  //   // if the dynamic simulation is not done (hit the end time)
  //   if ( bApplyAct_LoadShedding && (isteps == 2500 || isteps == 3000 ||
  //         isteps == 3500 || isteps == 4000 ||
  //         isteps == 4500 || isteps == 5000 || isteps == 5500 ) ){
  //     //apply action
  //     hadrec_app_sptr->applyAction(loadshedact);
  //     hadrec_app_sptr->applyAction(loadshedact1);
  //     //printf("----renke debug load shed, isteps: %d \n", isteps);
  //   }
	
	// if ( bApplyAct_LoadPchange && isteps == 400){  // apply a load change with 1000 MW, mimic fault
	// 	hadrec_app_sptr->applyAction(loadpchange);
	// }

  //   if ( bApplyAct_LineTripping && isteps == 400){
  //     if (me == 0) printf("----renke debug line trip, isteps: %d \n", isteps);
  //     hadrec_app_sptr->applyAction(linetrip);
  //   }
    //execute one dynamic simulation step
    hadrec_app_sptr->executeDynSimuOneStep();

// #ifdef USE_HELICS
// #if 0

    std::vector<int> buslist;
    std::vector<double> plist;
    std::vector<double> qlist;
    double ptmp = 700;
    double qtmp = 250;

    if(use_helics){
      if(localIndices.size() > 0 && network->getActiveBus(localIndices[0])){
        //pub.publish(widearea_deltafreq);

        for(int j = 0; j < pubCount; j++) {
          pub = (*fed).getPublication(j);
          std::string pubInfo = pub.getInfo();
          //std::cout << "-------------!!!helics test: HELICS pub info: " << pubInfo << std::endl;
          // pub.publish(widearea_deltafreq);
          gridpack::dynamic_simulation::DSFullBus *bus = dynamic_cast<gridpack::dynamic_simulation::DSFullBus*>(network->getBus(outnodeIndex).get());

          gridpack::ComplexType voltage = network->getBus(outnodeIndex)->getComplexVoltage();
          double rV = real(voltage);
          double iV = imag(voltage);
          double V = sqrt(rV*rV+iV*iV);
          double Ang = acos(rV/V);
          if (iV < 0) {
            Ang = -Ang;
          }

          double basevoltage = hadrec_app_sptr->getBaseVoltage(outnodeIndex);

          std::complex<double> voltage_cosim {V * basevoltage * 1000, 0};
          // std::complex<double> voltage_cosim{142300, 0};
          std::cout << "$$$$$voltage_cosim = " << voltage_cosim << "basevoltage" <<  basevoltage << std::endl;
          // std::complex<double> pub_info = {};
          pub.publish(voltage_cosim);
          // do stuff to tie pub to GridPACK object property
        }

        helics_requestTime = (double)(time_step * isteps);
        printf("-------------!!!Helics request time: %f \n", helics_requestTime);
        std::cout << "isteps: " << isteps << "time_step: " << time_step << std::endl; 
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
          std::cout << "subvalue:" << subvalue << std::endl;
            // printf("-------------!!!Helics sub value: %12.6f \n", subvalue);
                                  //update GridPACK object property with value
              //}
        }
        //printf("-------------!!!Outside Helics def sub value: %12.6f \n", subvalue);

        double load_amplification_factor = hadrec_app_sptr->getLoadAmplifier();

        // gridpack::powerflow::PFBus *bus = dynamic_cast<gridpack::powerflow::PFBus*>(p_network->getBus(i).get());
        ptmp = subvalue.real() * load_amplification_factor / 1000000;
        qtmp = subvalue.imag() * load_amplification_factor / 1000000;

        arr[0] = ptmp;
        arr[1] = qtmp;
      }
      MPI_Barrier(MPI_COMM_WORLD);
      ptmp = arr[0];
      qtmp = arr[1];

      double ppu_tmp = ptmp/100.0;
      double qpu_tmp = qtmp/100.0;

      std::cout << "line: 470, ppu_tmp, qpu_tmp: " << ppu_tmp << " " << qpu_tmp << ", me: " << me << std::endl;

      buslist.push_back(connectedBusID);
      plist.push_back(ppu_tmp);
      qlist.push_back(qpu_tmp);
      //  the scatterInjectionLoad will change the loads at the buslist to be the values in the plist and qlist
      //  the value of load P and Q should be p.u., based on 100 MVA system base
      hadrec_app_sptr->scatterInjectionLoadNew(buslist, plist, qlist);
    }
    
// #endif  //end if of HELICS
    ob_vals.clear();
    ob_vals = hadrec_app_sptr->getObservations();
    if (boutputob && me == 0) {
      printf("observations, ");
      for (idxtmp=0; idxtmp<ob_vals.size(); idxtmp++) {
        printf(" %16.12f, ", ob_vals[idxtmp]);
      }
      printf(" \n");
    }

    isteps++;
  }

// #ifdef USE_HELICS
// #if 0
  if(use_helics){
    if(localIndices.size() > 0 && network->getActiveBus(localIndices[0])){
      (*fed).finalize();
    }
  }
// #endif

  if (me == 0) printf("\n----------------finished first round of dynamic simulation----\n ");
  //timer->stop(t_total);
  //timer->dump();


  //start the reload and second time dynamic simulation here
  // transfer power flow results to dynamic simulation
  bool btest_2dynasimu = false;
  if (btest_2dynasimu) {

    gridpack::dynamic_simulation::Event busfault;
    busfault.start = 1.0;
    busfault.end = 1.2;
    busfault.step = 0.005;
    busfault.isBus = true;
    busfault.bus_idx = 7;

    BusFaults.clear();
    BusFaults.push_back(busfault);

    //hadrec_app_sptr->solvePowerFlowBeforeDynSimu(argc, argv);
    //hadrec_app_sptr->solvePowerFlowBeforeDynSimu(const_cast<char *>(file.c_str()), 1);
    hadrec_app_sptr->solvePowerFlowBeforeDynSimu(const_cast<char *>(file.c_str()), 0);

    if (me == 0) printf("\n---------------renke debug, hadrec main, second dyn starts----------------\n");
    hadrec_app_sptr->transferPFtoDS();

    // initialize dynamic simulation
    hadrec_app_sptr->initializeDynSimu(BusFaults);

    busno = 5;
    hadrec_app_sptr->getBusTotalLoadPower(busno, lp, lq);
    if (me == 0) printf("------------test hadrec_main, load at bus %d, has P: %f, Q: %f\n", busno, lp, lq);
    busno = 7;
    hadrec_app_sptr->getBusTotalLoadPower(busno, lp, lq);
    if (me == 0) printf("------------test hadrec_main, load at bus %d, has P: %f, Q: %f\n", busno, lp, lq);
    busno = 9;
    hadrec_app_sptr->getBusTotalLoadPower(busno, lp, lq);
    if (me == 0) printf("------------test hadrec_main, load at bus %d, has P: %f, Q: %f\n", busno, lp, lq);

    busno = 1;
    std::string genid = "1";
    btmp = hadrec_app_sptr->getGeneratorPower(busno, genid, pg, qg);
    if (me == 0) printf("------------test hadrec_main, %d find generator at bus %d, has P: %f, Q: %f\n", btmp, busno, pg, qg);

    busno = 2;
    genid = "1";
    btmp = hadrec_app_sptr->getGeneratorPower(busno, genid, pg, qg);
    if (me == 0) printf("------------test hadrec_main, %d find generator at bus %d, has P: %f, Q: %f\n", btmp, busno, pg, qg);

    busno = 3;
    genid = "1";
    btmp = hadrec_app_sptr->getGeneratorPower(busno, genid, pg, qg);
    if (me == 0) printf("------------test hadrec_main, %d find generator at bus %d, has P: %f, Q: %f\n", btmp, busno, pg, qg);

    isteps = 0;
    //bApplyAct_LoadShedding = true;  // whether apply the action in the simulation steps

    while(!hadrec_app_sptr->isDynSimuDone()){
      // if the dynamic simulation is not done (hit the end time)
      if ( bApplyAct_LoadShedding && (isteps == 2500 || isteps == 3000 ||
            isteps == 3500 || isteps == 4000 ) ){
        //apply action
        hadrec_app_sptr->applyAction(loadshedact);
        hadrec_app_sptr->applyAction(loadshedact1);
        //printf("----renke debug load shed, isteps: %d \n", isteps);
      }
      //execute one dynamic simulation step
      hadrec_app_sptr->executeDynSimuOneStep();

      ob_vals.clear();
      ob_vals = hadrec_app_sptr->getObservations();

      if (debugoutput && me == 0){
        printf("observations, ");
        for (idxtmp=0; idxtmp<ob_vals.size(); idxtmp++) {
          printf(" %16.12f, ", ob_vals[idxtmp]);
        }
        printf(" \n");
      }

      isteps++;
    }
  }

  if(me == 0){
  if(hadrec_app_sptr->isSecure() == -1){
    std::cout << "Power system is secure!" << std::endl;
  } else {
    std::cout << "Power system is not secure!" << std::endl;
  }

  printf(" ----------------- hadrec finished \n ");
  }
  // Make sure we could do it again with a new instance if we wanted to
  hadrec_app_sptr.reset(new gridpack::hadrec::HADRECAppModule());
  hadrec_app_sptr.reset();

  timer->stop(t_total);
  timer->dump();

  return 0;
}

