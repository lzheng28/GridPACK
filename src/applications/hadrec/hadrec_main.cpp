
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
#include <fstream>

#define ARRAY_LEN 512

#ifdef USE_HELICS
//#include "helics/ValueFederates.hpp"
//#include <helics/shared_api_library/ValueFederate.h>
#include <helics/helics.hpp>
#include <helics/helics98.hpp>
#endif

const char* help = "HADREC Test application";

void printBus(boost::shared_ptr<gridpack::dynamic_simulation::DSFullNetwork> network, int me){
  for(int i = 0; i < network->numBuses(); i++){
    std::cout << "Rank: " << me << "  " << network->getOriginalBusIndex(i);
    if(!network->getActiveBus(i)){
      std::cout << " GhostBus";
    }
    std::cout << std::endl;
  }
}

bool HELICS_DEBUG = true;

// helics value publication class
class helics_value_pub{
public:
  std::string key;
  std::string pubInfo;
  helicscpp::Publication HelicsPublication;
};

// helics value subscription class
class helics_value_sub{
public:
  std::string key;
  std::string subInfo;
  helicscpp::Input HelicsSubscription;
};

// clock class
class helics_clock{
public:
  double initial_time;
  double last_cosim_time;
  double next_cosim_time;
  double cosim_timestep;
  void step_forward();
};
// forward one step
void helics_clock::step_forward(){
  last_cosim_time += cosim_timestep;
  next_cosim_time += cosim_timestep;
}

// helics message class
class helics_msg{
public:
  std::vector<helics_value_pub*> helics_value_pubs;
  std::vector<helics_value_sub*> helics_value_subs;
  std::shared_ptr<helicscpp::CombinationFederate> helics_fed;
  helics_clock time;
  std::string configFile;
  int configure(std::string fileName);
  int init();
  void finalize();
  void publishVariables(std::vector<std::complex<double>>& pubVector);
  void subscribeVariables(std::vector<std::complex<double>>& subVector);
  void clockupdate(double current_time);
  int getPubCount();
  int getSubCount();
  void clockInit(double current_time, double helics_time_step);
};

int helics_msg::configure(std::string fileName){
  int returnValue = 1;
	if (fileName != "") {
		configFile = fileName;
	} else {
		std::cout << "Helics configuration file doesn't exist, please check it!" << std::endl;
		returnValue = 0;
	}
	return returnValue;
}

void helics_msg::finalize(){
  if(HELICS_DEBUG)
    std::cout << "helics_msg: Calling finalize\n" << std::endl;
  const helics_federate_state fed_state = helics_fed->getCurrentMode();
  if(fed_state != helics_state_finalize) {
    helics_fed->finalize();
  }
  helicscpp::cleanupHelicsLibrary();
}

int helics_msg::init(){
  if(helics_fed == NULL){
    try {
      helics_fed = std::make_shared<helicscpp::CombinationFederate>(configFile);
      int pub_count = helics_fed->getPublicationCount();
			int sub_count = helics_fed->getInputCount();
      int idx = 0;
      helics_value_pub *gpk_pub;
			helics_value_sub *gpk_sub;
      // pub of helics fed
      for( idx = 0; idx < pub_count; idx++ ) {
        helicscpp::Publication pub = helics_fed->getPublication(idx);
        if( pub.isValid() ) {
            gpk_pub = new helics_value_pub();
            gpk_pub->key = std::string(pub.getKey());
            gpk_pub->HelicsPublication = pub;
            gpk_pub->pubInfo = std::string(pub.getInfo());
            helics_value_pubs.push_back(gpk_pub);
        }
      }
      // sub of helics fed
      for( idx = 0; idx < sub_count; idx++ ) {
        helicscpp::Input sub = helics_fed->getSubscription(idx);
        if( sub.isValid() ) {
            gpk_sub = new helics_value_sub();
            gpk_sub->key = std::string(sub.getKey());
            gpk_sub->HelicsSubscription = sub;
            gpk_sub->subInfo = std::string(sub.getInfo());
            helics_value_subs.push_back(gpk_sub);
        }
      }

    } catch(const std::exception &e) {
      std::cout << "helics_msg::init(): configuration file error!" << std::endl;
    }
  }
  // TODO: pub type & sub type compare with complex data type

  // enterInitializationState
  helics_fed->enterInitializingMode();
  // enterExecutionState
  helics_fed->enterExecutingMode();
  // atexit(finalize);
  return 1;
}

void helics_msg::publishVariables(std::vector<std::complex<double>>& pubVector){
  std::complex<double> complex_temp = {0.0, 0.0};
  int idx = 0;
  int pubCount  = helics_value_pubs.size();
  if(pubVector.size() != pubCount){
    std::cout << "Publication size doesn't match conneceted node size, please check it!" << std::endl;
    return ;
  }
  for(idx = 0; idx < pubCount; idx++) {
    double real_part = pubVector[idx].real();
		double imag_part = pubVector[idx].imag();
    complex_temp = {real_part, imag_part};
    helics_value_pubs[idx]->HelicsPublication.publish(complex_temp);
  }
}

void helics_msg::subscribeVariables(std::vector<std::complex<double>>& subVector){
  std::complex<double> complex_temp = {0.0, 0.0};
  int idx = 0;
  int subCount  = helics_value_subs.size();
  if(subVector.size() != subCount){
    std::cout << "Subscription size doesn't match conneceted node size, please check it!" << std::endl;
    return ;
  }
  for(idx = 0; idx < subCount; idx++) {
    helicscpp::Input sub = helics_value_subs[idx]->HelicsSubscription;
    if(sub.isUpdated()) {
      complex_temp = sub.getComplex();
      subVector[idx].real(complex_temp.real());
      subVector[idx].imag(complex_temp.imag());
    }
  }
}

void helics_msg::clockupdate(double current_time){
  double helics_t = 0;
  if(current_time < time.next_cosim_time && current_time != 0){
    return ;
  }
  double helics_requestTime = current_time;
  double helics_grantime;
  helics_grantime = helics_fed->requestTime(helics_requestTime);
  if(HELICS_DEBUG)
    std::cout << "helics_msg: Granted " << (double)helics_grantime << std::endl;
  time.step_forward();
}

int helics_msg::getPubCount(){
  return helics_value_pubs.size();
}

int helics_msg::getSubCount(){
  return helics_value_subs.size();
}

void helics_msg::clockInit(double current_time, double helics_time_step){
  time.initial_time = current_time;
  time.last_cosim_time = -helics_time_step;
  time.next_cosim_time = current_time;
}


/*****************************************************************************************/
/************************************* hadrec.x main() ***********************************/
/*****************************************************************************************/
int main(int argc, char **argv)
{
  // Initialize libraries (parallel and math)
  gridpack::NoPrint *noprint_ins = gridpack::NoPrint::instance();

  gridpack::Environment env(argc,argv,help);
  // MPI communicator
  gridpack::parallel::Communicator world;
  int me = world.rank();
  int size = world.size();
  // Timer for GridPACK
  gridpack::utility::CoarseTimer *timer =
    gridpack::utility::CoarseTimer::instance();
  int t_total = timer->createCategory("Dynamic Simulation: Total Application");
  timer->start(t_total);

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
  // solve powerflow before DynamicSimualtion
  hadrec_app_sptr->solvePowerFlowBeforeDynSimu(const_cast<char *>(file.c_str()));

  // transfer power flow results to dynamic simulation
  hadrec_app_sptr->transferPFtoDS();
  // Declare BusFaults
  std::vector<gridpack::dynamic_simulation::Event> BusFaults;
  BusFaults.clear();
  
  // initialize dynamic simulation
  hadrec_app_sptr->initializeDynSimu(BusFaults);

  // <local index, original index> vector on each process 
  std::vector<std::pair<int, int>> outnodeIndexs;
  // network on each process
  boost::shared_ptr<gridpack::dynamic_simulation::DSFullNetwork> network = hadrec_app_sptr->ds_network;
  // dynamic simulation timestep
  double time_step = hadrec_app_sptr->getTimeStep();
  // use helics status
  bool use_helics = false;
  use_helics = hadrec_app_sptr->useHelicsStatus();

  // // helics federates shared pointer
  // std::shared_ptr<helicscpp::CombinationFederate> fed;
  // // helics publication
  // helicscpp::Publication pub;
  // // helics subscription
  // helicscpp::Input sub;
  // // send to helics broker to sync time
  // double helics_requestTime;
  // // helics broker return time
  // double helics_grantime;
  // complex number of subvalue
  std::complex<double> subvalue = {1.0, 1.0};
  // std::vector<std::complex<double>> subvalues;
  
  //to get publication definitions
  int pubCount;
  int subCount;

  int nBus = network->numBuses();

  std::vector<int> connectedBusIDs = hadrec_app_sptr->getHelicsConnectNodes();
  std::vector<std::vector<int>> localIndices_;
  for(int i = 0; i < connectedBusIDs.size(); i++){
    std::vector<int> tmp = network->getLocalBusIndices(connectedBusIDs[i]);
    if(tmp.size() > 0 && network->getActiveBus(tmp[0])){
      localIndices_.push_back(tmp); //Get local indices of connected nodes
      
      outnodeIndexs.push_back(std::make_pair(localIndices_.back()[0], connectedBusIDs[i])); 
    }
  }

  std::cout << "localIndices_.size() = " << localIndices_.size() << ", me = " << me << std::endl;

  int* sendcounts = (int*)malloc(sizeof(int) * size);
  sendcounts[me] = localIndices_.size(); //std::cout << "sendcounts.size() = " << sendcounts.size() << "  size = " << size << std::endl;
 
  // The displs
  int* displs = (int*)malloc(sizeof(int) * size);
  for(int ii = 1; ii < size; ii++){
    displs[ii] = displs[ii - 1] + sendcounts[ii - 1];
  }

  std::complex<double>* recvbuf = (std::complex<double>*)malloc(sizeof(std::complex<double>) * sendcounts[me]);

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

  helics_msg helicsMsg;
  if(use_helics){
    if(me == 0){
      std::cout << "-------------!!!helics test: HELICS Version: " << helics::versionString << std::endl;
      std::string configFile = hadrec_app_sptr->getHelicsConfigFile();

      helicsMsg.configure(configFile);
      helicsMsg.init();
      
      //to get publication definitions
      // fed = std::make_shared<helics::ValueFederate>(configFile);
      // pubCount = (*fed).getPublicationCount();

      pubCount = helicsMsg.getPubCount();
      printf("-------------helics test: num of pub: %d \n", pubCount);
      // for(int j = 0; j < pubCount; j++) {
      //   pub = (*fed).getPublication(j);
      //   std::string pubInfo = pub.getInfo();
      //   // do stuff to tie pub to GridPACK object property
      // }
      
      //to get subscription definitions
      // subCount = (*fed).getInputCount();
      subCount = helicsMsg.getSubCount();
      printf("-------------helics test: num of sub: %d \n", subCount);
      // std::vector<std::complex<double>> subvalues(subCount);
      // for(int j = 0; j < subCount; j++) {
      //   sub = (*fed).getInput(j);
      //   std::string subInfo = sub.getInfo();
      //   // do stuff to tie pub to GridPACK object property
      // }

      // //let helics broker know you are ready to start simulation 
      // (*fed).enterExecutingMode();
    }
  }

  bool debugoutput = false; // whether print out debug staffs
  bool boutputob = true;
  double lp, lq, pg, qg;
  int busno = 5;
  bool btmp;

  std::string observationFileName = hadrec_app_sptr->getObservationFileName();
  if(observationFileName == ""){
    boutputob = false;
  } else {
    boutputob = true;
  }

  //-----test get load and get generator function-----------
  int isteps = 0;

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

  std::vector<double> tmp_p;
  std::vector<double> tmp_q;

  std::ofstream outFile;
  
  char buf[1024];
  std::string outbuf;
  // Output observation
  if (boutputob && me == 0) {  
    outFile.open(observationFileName, std::ios::out);
    sprintf(buf, "time,  ");
    outbuf += buf;
    for (idxtmp=0; idxtmp<obs_genBus.size(); idxtmp++){
        sprintf(buf, "gen-at-bus-%d-ID-%s-speed,  ", obs_genBus[idxtmp], obs_genIDs[idxtmp].c_str());
        outbuf += buf;
    }
	
    for (idxtmp=0; idxtmp<obs_genBus.size(); idxtmp++){
        sprintf(buf, "gen-at-bus-%d-ID-%s-angle,  ", obs_genBus[idxtmp], obs_genIDs[idxtmp].c_str());
        outbuf += buf;
    }
	
	  for (idxtmp=0; idxtmp<obs_genBus.size(); idxtmp++){
        sprintf(buf, "gen-at-bus-%d-ID-%s-P,  ", obs_genBus[idxtmp], obs_genIDs[idxtmp].c_str());
        outbuf += buf;
    }
	
	  for (idxtmp=0; idxtmp<obs_genBus.size(); idxtmp++){
      sprintf(buf, "gen-at-bus-%d-ID-%s-Q,  ", obs_genBus[idxtmp], obs_genIDs[idxtmp].c_str());
      outbuf += buf;
    }
	
    for (idxtmp=0; idxtmp<obs_vBus.size(); idxtmp++){
      sprintf(buf, "bus-%d-magnitude, ", obs_vBus[idxtmp]);
      outbuf += buf;
    }
	
	  for (idxtmp=0; idxtmp<obs_vBus.size(); idxtmp++){
      sprintf(buf, "bus-%d-angle, ", obs_vBus[idxtmp]);
      outbuf += buf;
    }

    for (idxtmp=0; idxtmp<obs_loadBus.size(); idxtmp++){
      sprintf(buf, "remaining-load-at-bus-%d-ID-%s,  ", obs_loadBus[idxtmp], obs_loadIDs[idxtmp].c_str());
      outbuf += buf;
    }
    sprintf(buf, " \n");
    outbuf += buf;
    outFile << outbuf;
  }
  double co_sim_time_interval = hadrec_app_sptr->getCosimTimeInterval(); // transfer data time step
  double co_sim_next_step_time = 0.0;

  // std::cout << "co_sim_time_interval: " << co_sim_time_interval << std::endl;

  while(!hadrec_app_sptr->isDynSimuDone()){
    //execute one dynamic simulation step
    hadrec_app_sptr->executeDynSimuOneStep();

    std::vector<int> buslist;
    std::vector<double> plist;
    std::vector<double> qlist;
    double ptmp = 700;
    double qtmp = 250;
    
    if(use_helics && isteps * time_step == co_sim_next_step_time){
      co_sim_next_step_time += co_sim_time_interval;
      int* grecvcounts = (int*)malloc(sizeof(int) * size);
      grecvcounts[me] = outnodeIndexs.size() * 2;

      arr[me] = outnodeIndexs.size() * 2;
      MPI_Barrier(MPI_COMM_WORLD);

      for(int i = 0; i < size; i++){
        grecvcounts[i] = arr[i];
      }

      int* gdispls = (int*)malloc(sizeof(int) * size);
      gdispls[0] = 0;
      for(int i = 1; i < size; i++){
        gdispls[i] = gdispls[i - 1] + grecvcounts[i - 1];
      }
      double* gsendbuf = (double*)malloc(sizeof(double) * outnodeIndexs.size() * 2);
      double* grecvbuf;

      for(int i = 0; i < outnodeIndexs.size() * 2; i = i + 2){
        gridpack::dynamic_simulation::DSFullBus *bus = dynamic_cast<gridpack::dynamic_simulation::DSFullBus*>(network->getBus(outnodeIndexs[i / 2].first).get());

        gridpack::ComplexType voltage = network->getBus(outnodeIndexs[i / 2].first)->getComplexVoltage();
        double rV = real(voltage);
        double iV = imag(voltage);
        double V = sqrt(rV*rV+iV*iV);
        double Ang = acos(rV/V);
        if (iV < 0) {
          Ang = -Ang;
        }

        double basevoltage = hadrec_app_sptr->getBaseVoltage(outnodeIndexs[i / 2].first);
        gsendbuf[i] = outnodeIndexs[i / 2].second;
        gsendbuf[i + 1] = V * basevoltage * 1000;
      }
      if(me == 0){
        grecvbuf = (double*)malloc(sizeof(double) * connectedBusIDs.size() * 10);
        for(int i = 0; i < 10; i++){
          grecvbuf[i] = 135000;
        }
      }

      MPI_Gatherv(gsendbuf, grecvcounts[me], MPI_DOUBLE, grecvbuf, grecvcounts, gdispls, MPI_DOUBLE, 0, MPI_COMM_WORLD);    
      if(me == 0){
        for(int k = 0; k < outnodeIndexs.size() * 2; k++){
          // std::cout << "grecvbuf: " << grecvbuf[k] << std::endl;
        }
        std::vector<double> voltage_v(pubCount, 0);
        for(int k = 0; k < pubCount; k++){
          for(int l = 0; l < pubCount * 2; l = l + 2){
            if(round(grecvbuf[l]) == connectedBusIDs[k]){
              voltage_v[k] = grecvbuf[l + 1];
              break;
            }
          }
        }

        std::vector<std::complex<double>> pubValueVector(pubCount);
        for(int j = 0; j < pubCount; j++) {
          // pub = (*fed).getPublication(j);
          std::complex<double> voltage_cosim {voltage_v[j], 0};
          pubValueVector[j] = voltage_cosim;
          // pub.publish(voltage_cosim);
        }
        helicsMsg.publishVariables(pubValueVector);
        helicsMsg.clockupdate((double)(time_step * isteps));

        // helics_requestTime = (double)(time_step * isteps);
        // helics_grantime = (*fed).requestTime(helics_requestTime);
        //  printf("-------------!!!Helics grant time: %12.6f \n", helics_grantime); 
        
        // for(int j = 0; j < subCount; j++) {
        //   sub = (*fed).getInput(j);
        //   sub.getValue(subvalue);
        //   subvalues.push_back(subvalue);
        // }
        std::vector<std::complex<double>> subValueVector(subCount);
        helicsMsg.subscribeVariables(subValueVector);
        for(int j = 0; j < subCount && HELICS_DEBUG; j++){
          printf("-------------!!!Outside Helics def pub value: %12.6f + %12.6f i \n", pubValueVector[j].real(), pubValueVector[j].imag());
          printf("-------------!!!Outside Helics def sub value: %12.6f + %12.6f i \n", subValueVector[j].real(), subValueVector[j].imag());
        }
        arr[0] = subCount;
        for(int j = 1; j <= subCount; j++){
        double load_amplification_factor = hadrec_app_sptr->getLoadAmplifier();

        // gridpack::powerflow::PFBus *bus = dynamic_cast<gridpack::powerflow::PFBus*>(p_network->getBus(i).get());
        ptmp = subValueVector[j - 1].real() * load_amplification_factor / 1000000;
        qtmp = subValueVector[j - 1].imag() * load_amplification_factor / 1000000;

        arr[j] = ptmp;
        arr[j + subCount] = qtmp; std::cout << "line: 477, ptmp, qtmp: " << ptmp << " " << qtmp << ", me: " << me << std::endl;
        }
        
      }

      // MPI_Scatterv(sendbuf, sendcounts, displs, MPI_COMPLEX, recvbuf, sendcounts[me], MPI_COMPLEX, 0, MPI_COMM_WORLD);

      MPI_Barrier(MPI_COMM_WORLD);subCount = arr[0];
      for(int j = 1; j <= subCount; j++){
        ptmp = arr[j];
        qtmp = arr[j + subCount];

        double ppu_tmp = ptmp/100.0;
        double qpu_tmp = qtmp/100.0;

        // std::cout << "line: 488, ppu_tmp, qpu_tmp: " << ppu_tmp << " " << qpu_tmp << ", me: " << me << std::endl;

        buslist.push_back(connectedBusIDs[j - 1]);
        plist.push_back(ppu_tmp);
        qlist.push_back(qpu_tmp);
      }
      //  the scatterInjectionLoad will change the loads at the buslist to be the values in the plist and qlist
      //  the value of load P and Q should be p.u., based on 100 MVA system base
      for(int j = 0; j < buslist.size() && HELICS_DEBUG; j++){
        std::cout << "plist[" << j << "] = " << plist[j] << std::endl;
        std::cout << "qlist[" << j << "] = " << qlist[j] << std::endl;
        std::cout << "buslist[" << j << "] = " << buslist[j] << std::endl;
      }
      hadrec_app_sptr->scatterInjectionLoadNew(buslist, plist, qlist);
      // subvalues = {};
    }
    
    // Save observation file
    ob_vals.clear();
    ob_vals = hadrec_app_sptr->getObservations();
    if (boutputob && me == 0) {
      outbuf = "";
      sprintf(buf, "%.4f, ", time_step * isteps);
      outbuf += buf;
      for (idxtmp=0; idxtmp<ob_vals.size(); idxtmp++) {
        sprintf(buf, " %16.12f, ", ob_vals[idxtmp]);
        outbuf += buf;
      }
      sprintf(buf, " \n");
      outbuf += buf;
    }
    outFile << outbuf;

    isteps++;
  }

  outFile.close(); // close and save observation file

  if(me == 0){
    if(use_helics){
      // (*fed).finalize();
      helicsMsg.finalize();
    }
  }

  // Print if the system is secure
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

