/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
/*
 *  Created on: December 30, 2014
 *      Author: Kevin Glass, Bruce Palmer
 */

#ifndef BASEPTIPARSER_HPP_
#define BASEPTIPARSER_HPP_

#define OLD_MAP

#include <iostream>
#include <string>
#include <vector>
#include <map>
#ifndef OLD_MAP
#include <boost/unordered_map.hpp>
#endif


#include "gridpack/component/base_component.hpp"
#include "gridpack/timer/coarse_timer.hpp"
#include "gridpack/component/data_collection.hpp"
#include "gridpack/parser/dictionary.hpp"
#include "gridpack/utilities/exception.hpp"
#include "gridpack/utilities/string_utils.hpp"
#include "gridpack/network/base_network.hpp"
#include "gridpack/parser/base_parser.hpp"
#include "gridpack/parser/hash_distr.hpp"
#include "parser_classes/gencls.hpp"
#include "parser_classes/gensal.hpp"
#include "parser_classes/genrou.hpp"
#include "parser_classes/wsieg1.hpp"
#include "parser_classes/exdc1.hpp"
#include "parser_classes/esst1a.hpp"
#include "parser_classes/esst4b.hpp"
#include "parser_classes/ggov1.hpp"
#include "parser_classes/wshygp.hpp"
#include "parser_classes/lvshbl.hpp"
#include "parser_classes/frqtpat.hpp"
#include "parser_classes/distr1.hpp"
#include "parser_classes/cim6bl.hpp"
#include "parser_classes/acmtblu1.hpp"
#include "parser_classes/ieelbl.hpp"
#include "parser_classes/cmldblu1.hpp"

namespace gridpack {
namespace parser {

template <class _network>
class BasePTIParser : public BaseParser<_network>
{
  public:

    /**
     * Constructor
     */
    explicit BasePTIParser()
    {
      p_timer = gridpack::utility::CoarseTimer::instance();
    }


    /**
     * Destructor
     */
    virtual ~BasePTIParser(){}

    /**
     * Parse a second file after original network has been distributed. This
     * requires the data in the second file to be distributed to all network
     * objects that need the data
     * @param fileName name of file
     */
    void externalParse(const std::string &fileName)
    {
      std::string ext = getExtension(fileName);
      if (ext == "dyr") {
        getDSExternal(fileName);
      } else if (ext == "uc") {
        getUCExternal(fileName);
      }
    }

  protected:

    /* ************************************************************************
     **************************************************************************
     ***** PROTECTED SCOPE
     **************************************************************************
     *********************************************************************** */

    /**
     * Assign network to internal network pointer variable
     */
    void setNetwork(boost::shared_ptr<_network> network)
    {
      p_network = network;
      BaseParser<_network>::setNetwork(network);
    }

    /**
     * This routine opens up a .dyr file with parameters for dynamic
     * simulation. It assumes that a .raw file has already been parsed
     */
    void getDS(const std::string & fileName)
    {
      int t_ds = p_timer->createCategory("Parser:getDS");
      p_timer->start(t_ds);
      int me(p_network->communicator().rank());

      if (me == 0) {
        std::ifstream input;
        input.open(fileName.c_str());
        if (!input.is_open()) {
          p_timer->stop(t_ds);
          return;
        }
        find_ds_par(input);
        input.close();
      }
      p_timer->stop(t_ds);
#if 0
      int i;
      printf("BUS data size: %d\n",p_busData.size());
      for (i=0; i<p_network->numBuses(); i++) {
        printf("Dumping bus: %d\n",i);
        p_network->getBusData(i)->dump();
      }
#endif
    }

    // Data structure to hold generator params
    struct gen_params{
      // Generator parameters
      int bus_id; // ID of bus that owns device
      char model[8];  // Model represented by data
      char gen_id[3]; // Generator ID
      double inertia;  // Inertia constant 0
      double damping;  // Damping coefficient
      double reactance; // Transient reactance
      double tdop;
      double tdopp;
      double tqop;
      double tqopp;
      double gn_xd;
      double gn_xq;
      double xdp;
      double xqp;
      double xdpp;
      double gn_xl;
      double gn_s1;
      double s12;
      // Exciter parameters
      bool has_exciter;
      double ex_tr;
      double ex_ka;
      double ex_ta;
      double ex_tb;
      double ex_tc;
      double vrmax;
      double vrmin;
      double ex_ke;
      double ex_te;
      double ex_kf;
      double tf1;
      double rswitch;
      double ex_e1;
      double se1;
      double ex_e2;
      double se2;
      double uel;
      double vos;
      double vimax;
      double vimin;
      double tc1;
      double tb1;
      double vamax;
      double vamin;
      double ex_kc;
      double ex_tf;
      double klr;
      double ilr;
      double kpr;
      double kir;
      double kpm;
      double kim;
      double vmmax;
      double vmmin;
      double ex_kg;
      double ex_kp;
      double ex_ki;
      double vbmax;
      double ex_xl;
      double thetap;
      // Governor parameters
      bool has_governor;
      int jbus;
      int gv_m;
      double gv_k;
      double gv_t1;
      double gv_t2;
      double gv_t3;
      double gv_uo;
      double gv_uc;
      double pmax;
      double pmin;
      double gv_t4;
      double gv_k1;
      double gv_k2;
      double gv_t5;
      double gv_k3;
      double gv_k4;
      double gv_t6;
      double gv_k5;
      double gv_k6;
      double gv_t7;
      double gv_k7;
      double gv_k8;
      double db1;
      double err;
      double db2;
      double gv1;
      double pgv1;
      double gv2;
      double pgv2;
      double gv3;
      double pgv3;
      double gv4;
      double pgv4;
      double gv5;
      double pgv5;
      int iblock;
      double rselect;
      double flagswitch;
      double gv_r;
      double tpelec;
      double maxerr;
      double minerr;
      double kpgov;
      double kigov;
      double kdgov;
      double tdgov;
      double vmax;
      double vmin;
      double tact;
      double kturb;
      double wfnl;
      double gv_tb;
      double gv_tc;
      double teng;
      double tfload;
      double kpload;
      double kiload;
      double ldref;
      double gv_dm;
      double ropen;
      double rclose;
      double kimw;
      double aset;
      double gv_ka;
      double gv_ta;
      double trate;
      double gv_db;
      double tsa;
      double tsb;
      double rup;
      double rdown;
      double gv_td;
      double gv_ki;
      double gv_tf;
      double gv_kd;
      double gv_kp;
      double gv_tt;
      double gv_kg;
      double gv_tp;
      double velopen;
      double velclose;
      double aturb;
      double bturb;
      double tturb;
    };

    // Data structure to hold relay parameters on buses
    struct bus_relay_params{
      int bus_id; // ID of bus that owns device
      char model[8];  // Model represented by data
      char tag[3]; //
      int jbus;
      int mins;
      int frebus;
      double v1;
      double t1;
      double f1;
      double v2;
      double t2;
      double f2;
      double v3;
      double t3;
      double f3;
      double tb;
      double fl;
      double fu;
      double tp;
    };

    // Data structure to hold relay parameters on branches
    struct branch_relay_params{
      int from_bus; // ID of from bus
      int to_bus; // ID of to bus
      char model[8];  // Model represented by data
      char id[3];
      char id1[3];
      char id2[3];
      char id3[3];
      int ibus;
      int jbus;
      int ibus1;
      int jbus1;
      int ibus2;
      int jbus2;
      int ibus3;
      int jbus3;
      int rs;
      int mtype;
      int bmon;
      int bltype1;
      int bltype2;
      double zone1_time;
      double zone1_reach;
      double zone1_cenang;
      double zone1_cendis;
      double zone2_time;
      double zone2_reach;
      double zone2_cenang;
      double zone2_cendis;
      double zone3_time;
      double zone3_reach;
      double zone3_cenang;
      double zone3_cendis;
      double dirang;
      double thcur;
      double sebtime;
      double serctime;
      double trbtime;
      double trrctime;
      double blint1;
      double blro1;
      double blint2;
      double blro2;
    };

    // Data structure to hold parameters for conventional loads
    struct load_params{
      int bus_id; // ID of bus that owns load
      char model[9];  // Model represented by data
      char id[3]; // Load ID
      int it;
      double ra;
      double xa;
      double xm;
      double r1;
      double x1;
      double r2;
      double x2;
      double e1;
      double se1;
      double e2;
      double se2;
      double mbase;
      double pmult;
      double h;
      double vi;
      double ti;
      double tb;
      double a;
      double b;
      double d;
      double e;
      double c0;
      double tnom;

      double tstall;
      double trestart;
      double tv;
      double tf;
      double complf;
      double comppf;
      double vstall;
      double rstall;
      double xstall;
      double lfadj;
      double kp1;
      double np1;
      double kq1;
      double nq1;
      double kp2;
      double np2;
      double kq2;
      double nq2;
      double vbrk;
      double frst;
      double vrst;
      double cmpkpf;
      double cmpkqf;
      double vc1off;
      double vc2off;
      double vc1on;
      double vc2on;
      double tth;
      double th1t;
      double th2t;
      double fuvr;
      double uvtr1;
      double ttr1;
      double uvtr2;
      double ttr2;

      double a1;
      double a2;
      double a3;
      double a4;
      double a5;
      double a6;
      double a7;
      double a8;
      double n1;
      double n2;
      double n3;
      double n4;
      double n5;
      double n6;

      double mva;
      double bss;
      double rfdr;
      double xfdr;
      double fb;
      double xxf;
      double tfixhs;
      double tfixls;
      double ltc;
      double tmin;
      double tmax;
      double step;
      double vmin;
      double vmax;
      double tdel;
      double ttap;
      double rcomp;
      double xcomp;
      double fma;
      double fmb;
      double fmc;
      double fmd;
      double fel;
      double pfel;
      double vd1;
      double vd2;
      double frcel;
      double pfs;
      double p1e;
      double p1c;
      double p2e;
      double p2c;
      double pfreq;
      double q1e;
      double q1c;
      double q2e;
      double q2c;
      double qfreq;

      int mtpa;
      double lfma;
      double rsa;
      double lsa;
      double lpa;
      double lppa;
      double tpoa;
      double tppoa;
      double ha;
      double etrqa;
      double vtr1a;
      double ttr1a;
      double ftr1a;
      double vrc1a;
      double trc1a;
      double vtr2a;
      double ttr2a;
      double ftr2a;
      double vrc2a;
      double trc2a;

      int mtpb;
      double lfmb;
      double rsb;
      double lsb;
      double lpb;
      double lppb;
      double tpob;
      double tppob;
      double hb;
      double etrqb;
      double vtr1b;
      double ttr1b;
      double ftr1b;
      double vrc1b;
      double trc1b;
      double vtr2b;
      double ttr2b;
      double ftr2b;
      double vrc2b;
      double trc2b;

      int mtpc;
      double lfmc;
      double rsc;
      double lsc;
      double lpc;
      double lppc;
      double tpoc;
      double tppoc;
      double hc;
      double etrqc;
      double vtr1c;
      double ttr1c;
      double ftr1c;
      double vrc1c;
      double trc1c;
      double vtr2c;
      double ttr2c;
      double ftr2c;
      double vrc2c;
      double trc2c;

      int mtpd;
      double lfmd;
      double trst;
      double vtr1;
      double vtr2;
    };

    /**
     * This routine opens up a .dyr file with parameters for dynamic
     * simulation and distributes the parameters to whatever processor holds the
     * corresponding buses. It assumes that a .raw file has already been parsed
     */
    void getDSExternal(const std::string & fileName)
    {

      //      int t_ds = p_timer->createCategory("Parser:getDS");
      //      p_timer->start(t_ds);
      int me(p_network->communicator().rank());

      std::vector<gen_params> gen_data;
      std::vector<bus_relay_params> bus_relay_data;
      std::vector<branch_relay_params> branch_relay_data;
      std::vector<load_params> load_data;
      if (me == 0) {
        std::ifstream            input;
        input.open(fileName.c_str());
        if (!input.is_open()) {
          // p_timer->stop(t_ds);
          return;
        }
        find_ds_vector(input, &gen_data, &bus_relay_data,
            &branch_relay_data, &load_data);
        input.close();
      }
      int nsize = gen_data.size();
      std::vector<int> buses;
      int i;
      for (i=0; i<nsize; i++) {
        buses.push_back(gen_data[i].bus_id);
      }
      gridpack::hash_distr::HashDistribution<_network,gen_params,gen_params>
        distr(p_network);
      distr.distributeBusValues(buses,gen_data);
      // Now match data with corresponding data collection objects
      gridpack::component::DataCollection *data;
      nsize = buses.size();
      for (i=0; i<nsize; i++) {
        int l_idx = buses[i];
        data = dynamic_cast<gridpack::component::DataCollection*>
          (p_network->getBusData(l_idx).get());

        // Find out how many generators are already on bus
        int ngen = 0;
        data->getValue(GENERATOR_NUMBER, &ngen);
        // Identify index of generator to which this data applies
        int g_id = -1;
        if (ngen > 0) {
          // Clean up 2 character tag for generator ID
          std::string tag = gen_data[i].gen_id;
          int j;
          for (j=0; j<ngen; j++) {
            std::string t_id;
            data->getValue(GENERATOR_ID,&t_id,j);
            if (tag == t_id) {
              g_id = j;
              break;
            }
          }
        }

        // Check to see parameters can be assigned to a generator
        if (g_id > -1) {
          if (!strcmp(gen_data[i].model,"GENCLS")) {
            GenclsParser<gen_params> parser;
            parser.extract(gen_data[i], data, g_id);
          } else if (!strcmp(gen_data[i].model,"GENSAL")) {
            GensalParser<gen_params> parser;
            parser.extract(gen_data[i], data, g_id);
          } else if (!strcmp(gen_data[i].model,"GENROU")) {
            GenrouParser<gen_params> parser;
            parser.extract(gen_data[i], data, g_id);
          } else if (!strcmp(gen_data[i].model,"WSIEG1")) {
            Wsieg1Parser<gen_params> parser;
            parser.extract(gen_data[i], data, g_id);
          } else if (!strcmp(gen_data[i].model,"EXDC1") ||
              !strcmp(gen_data[i].model,"EXDC2")) {
            Exdc1Parser<gen_params> parser;
            parser.extract(gen_data[i], data, g_id);
          } else if (!strcmp(gen_data[i].model,"ESST1A")) {
            Esst1aParser<gen_params> parser;
            parser.extract(gen_data[i], data, g_id);
          } else if (!strcmp(gen_data[i].model,"ESST4B")) {
            Esst4bParser<gen_params> parser;
            parser.extract(gen_data[i], data, g_id);
          } else if (!strcmp(gen_data[i].model,"GGOV1")) {
            Ggov1Parser<gen_params> parser;
            parser.extract(gen_data[i], data, g_id);
          } else if (!strcmp(gen_data[i].model,"WSHYGP")) {
            WshygpParser<gen_params> parser;
            parser.extract(gen_data[i], data, g_id);
          }
        }
      }
      // Add parameters for a bus relay
      nsize = bus_relay_data.size();
      buses.clear();
      for (i=0; i<nsize; i++) {
        buses.push_back(bus_relay_data[i].bus_id);
      }
      gridpack::hash_distr::HashDistribution<_network,bus_relay_params,
        bus_relay_params> distr_bus(p_network);
      distr_bus.distributeBusValues(buses,bus_relay_data);
      // Now match data with corresponding data collection objects
      nsize = buses.size();
      for (i=0; i<nsize; i++) {
        int l_idx = buses[i];
        data = dynamic_cast<gridpack::component::DataCollection*>
          (p_network->getBusData(l_idx).get());
        // Add relay data to bus
        if (!strcmp(bus_relay_data[i].model,"LVSHBL")) {
          LvshblParser<bus_relay_params> parser;
          parser.extract(bus_relay_data[i], data);
        } else if (!strcmp(bus_relay_data[i].model,"FRQTPAT")) {
          FrqtpatParser<bus_relay_params> parser;
          parser.extract(bus_relay_data[i], data);
        }
      }
      buses.clear();
      // Add parameters for a branch relay
      nsize = branch_relay_data.size();
      std::vector<std::pair<int,int> > branches;
      std::vector<int> lbranch;
      for (i=0; i<nsize; i++) {
        branches.push_back(std::pair<int,int>(branch_relay_data[i].from_bus,
              branch_relay_data[i].to_bus));
      }
      gridpack::hash_distr::HashDistribution<_network,branch_relay_params,
        branch_relay_params> distr_branch(p_network);
      distr_branch.distributeBranchValues(branches,lbranch,branch_relay_data);
      // Now match data with corresponding data collection objects
      nsize = lbranch.size();
      for (i=0; i<nsize; i++) {
        int l_idx = lbranch[i];
        data = dynamic_cast<gridpack::component::DataCollection*>
          (p_network->getBranchData(l_idx).get());
        // Add relay data to branch
        if (!strcmp(branch_relay_data[i].model,"DISTR1")) {
          Distr1Parser<branch_relay_params> parser;
          parser.extract(branch_relay_data[i], data);
        }
      }

      // Add parameters for a load
      nsize = load_data.size();
      buses.clear();
      for (i=0; i<nsize; i++) {
        buses.push_back(load_data[i].bus_id);
      }
      gridpack::hash_distr::HashDistribution<_network,load_params,
        load_params> distr_load(p_network);
      distr_load.distributeBusValues(buses,load_data);
      // Now match data with corresponding data collection objects
      nsize = buses.size();
      for (i=0; i<nsize; i++) {
        int l_idx = buses[i];
        data = dynamic_cast<gridpack::component::DataCollection*>
          (p_network->getBusData(l_idx).get());

        // Find out how many generators are already on bus
        int nload = 0;
        data->getValue(LOAD_NUMBER, &nload);
        // Identify index of generator to which this data applies
        int l_id = -1;
        if (nload > 0) {
          // Clean up 2 character tag for generator ID
          std::string tag = load_data[i].id;
          int j;
          for (j=0; j<nload; j++) {
            std::string t_id;
            data->getValue(LOAD_ID,&t_id,j);
            if (tag == t_id) {
              l_id = j;
              break;
            }
          }
        }

        // Assign parameters to a load
        if (l_id > -1) {
          if (!strcmp(bus_relay_data[i].model,"CIM6BL")) {
            Cim6blParser<load_params> parser;
            parser.extract(load_data[i], data, l_id);
          } else if (!strcmp(bus_relay_data[i].model,"IEELBL")) {
            IeelblParser<load_params> parser;
            parser.extract(load_data[i], data, l_id);
          } else if (!strcmp(bus_relay_data[i].model,"ACMTBLU1")) {
            Acmtblu1Parser<load_params> parser;
            parser.extract(load_data[i], data, l_id);
          } else if (!strcmp(bus_relay_data[i].model,"CMLDBLU1")) {
            Cmldblu1Parser<load_params> parser;
            parser.extract(load_data[i], data, l_id);
          }
        }
      }
    }

    struct uc_params{
      int bus_id; // ID of bus that owns generator
      char gen_id[3]; // Generator ID
      int type;
      double init_level; // Initial production level
      double min_gen; // Minimum generation
      double max_gen; // Maximum generation
      double max_oper; // Maximum operating generation
      double min_up;
      double min_down;
      double ramp_up;
      double ramp_down;
      double start_up; // Start up cost
      double const_cost; // Constant cost
      double lin_cost; // Linear cost
      double co_2_cost;
      double init_prd; // Init periods
      double start_cap; // Startup cap
      double shut_cap; // Shutdown cap
    };

    /**
     * This routine opens up a .uc file with parameters for a unit commitment
     * calculationand distributes the parameters to whatever processor holds the
     * corresponding buses. It assumes that a .raw file has already been parsed
     */
    void getUCExternal(const std::string & fileName)
    {

      int me(p_network->communicator().rank());

      std::vector<uc_params> uc_data;
      if (me == 0) {
        std::ifstream            input;
        input.open(fileName.c_str());
        if (!input.is_open()) {
          return;
        }
        find_uc_vector(input, &uc_data);
        input.close();
      }
      int nsize = uc_data.size();
      std::vector<int> buses;
      int i;
      for (i=0; i<nsize; i++) {
        buses.push_back(uc_data[i].bus_id);
      }
      gridpack::hash_distr::HashDistribution<_network,uc_params,uc_params>
        distr(p_network);
      distr.distributeBusValues(buses,uc_data);
      // Now match data with corresponding data collection objects
      gridpack::component::DataCollection *data;
      nsize = buses.size();
      for (i=0; i<nsize; i++) {
        int l_idx = buses[i];
        data = dynamic_cast<gridpack::component::DataCollection*>
          (p_network->getBusData(l_idx).get());

        // Find out how many generators are already on bus
        int ngen;
        if (!data->getValue(GENERATOR_NUMBER, &ngen)) continue;
        // Identify index of generator to which this data applies
        int g_id = -1;
        // Clean up 2 character tag for generator ID
        std::string tag = uc_data[i].gen_id;
        int j;
        for (j=0; j<ngen; j++) {
          std::string t_id;
          data->getValue(GENERATOR_ID,&t_id,j);
          if (tag == t_id) {
            g_id = j;
            break;
          }
        }
        if (g_id == -1) continue;

        double rval;
        int ival;
        if (!data->getValue("GENERATOR_TYPE",&ival,g_id)) {
          data->addValue("GENERATOR_TYPE", uc_data[i].type, g_id);
        } else {
          data->setValue("GENERATOR_TYPE", uc_data[i].type, g_id);
        }

        if (!data->getValue("GENERATOR_INIT_LEVEL",&rval,g_id)) {
          data->addValue("GENERATOR_INIT_LEVEL", uc_data[i].init_level, g_id);
        } else {
          data->setValue("GENERATOR_INIT_LEVEL", uc_data[i].init_level, g_id);
        }

        if (!data->getValue("GENERATOR_MIN_GEN",&rval,g_id)) {
          data->addValue("GENERATOR_MIN_GEN", uc_data[i].min_gen, g_id);
        } else {
          data->setValue("GENERATOR_MIN_GEN", uc_data[i].min_gen, g_id);
        }

        if (!data->getValue("GENERATOR_MAX_GEN",&rval,g_id)) {
          data->addValue("GENERATOR_MAX_GEN", uc_data[i].max_gen, g_id);
        } else {
          data->setValue("GENERATOR_MAX_GEN", uc_data[i].max_gen, g_id);
        }

        if (!data->getValue("GENERATOR_MAX_OP_GEN",&rval,g_id)) {
          data->addValue("GENERATOR_MAX_OP_GEN", uc_data[i].max_gen, g_id);
        } else {
          data->setValue("GENERATOR_MAX_OP_GEN", uc_data[i].max_gen, g_id);
        }
        if (!data->getValue("GENERATOR_MIN_UP",&rval,g_id)) {
          data->addValue("GENERATOR_MIN_UP", uc_data[i].min_up, g_id);
        } else {
          data->setValue("GENERATOR_MIN_UP", uc_data[i].min_up, g_id);
        }

        if (!data->getValue("GENERATOR_MIN_DOWN",&rval,g_id)) {
          data->addValue("GENERATOR_MIN_DOWN", uc_data[i].min_down, g_id);
        } else {
          data->setValue("GENERATOR_MIN_DOWN", uc_data[i].min_down, g_id);
        }

        if (!data->getValue("GENERATOR_RAMP_UP",&rval,g_id)) {
          data->addValue("GENERATOR_RAMP_UP", uc_data[i].ramp_up, g_id);
        } else {
          data->setValue("GENERATOR_RAMP_UP", uc_data[i].ramp_up, g_id);
        }

        if (!data->getValue("GENERATOR_RAMP_DOWN",&rval,g_id)) {
          data->addValue("GENERATOR_RAMP_DOWN", uc_data[i].ramp_down, g_id);
        } else {
          data->setValue("GENERATOR_RAMP_DOWN", uc_data[i].ramp_down, g_id);
        }

        if (!data->getValue("GENERATOR_START_UP",&rval,g_id)) {
          data->addValue("GENERATOR_START_UP", uc_data[i].start_up, g_id);
        } else {
          data->setValue("GENERATOR_START_UP", uc_data[i].start_up, g_id);
        }

        if (!data->getValue("GENERATOR_CONST_COST",&rval,g_id)) {
          data->addValue("GENERATOR_CONST_COST", uc_data[i].const_cost, g_id);
        } else {
          data->setValue("GENERATOR_CONST_COST", uc_data[i].const_cost, g_id);
        }
        if (!data->getValue("GENERATOR_LIN_COST",&rval,g_id)) {
          data->addValue("GENERATOR_LIN_COST", uc_data[i].lin_cost, g_id);
        } else {
          data->setValue("GENERATOR_LIN_COST", uc_data[i].lin_cost, g_id);
        }
        if (!data->getValue("GENERATOR_CO_2_COST",&rval,g_id)) {
          data->addValue("GENERATOR_CO_2_COST", uc_data[i].co_2_cost, g_id);
        } else {
          data->setValue("GENERATOR_CO_2_COST", uc_data[i].co_2_cost, g_id);
        }

        if (!data->getValue("GENERATOR_INIT_PRD",&rval,g_id)) {
          data->addValue("GENERATOR_INIT_PRD", uc_data[i].init_prd, g_id);
        } else {
          data->setValue("GENERATOR_INIT_PRD", uc_data[i].init_prd, g_id);
        }

        if (!data->getValue("GENERATOR_START_CAP",&rval,g_id)) {
          data->addValue("GENERATOR_START_CAP", uc_data[i].start_cap, g_id);
        } else {
          data->setValue("GENERATOR_START_CAP", uc_data[i].start_cap, g_id);
        }

        if (!data->getValue("GENERATOR_SHUT_CAP",&rval,g_id)) {
          data->addValue("GENERATOR_SHUT_CAP", uc_data[i].shut_cap, g_id);
        } else {
          data->setValue("GENERATOR_SHUT_CAP", uc_data[i].shut_cap, g_id);
        }
      }
    }

    // Utility function to check if device is on a generator
    bool onGenerator(std::string &device) {
      bool ret = false;
      if (device == "GENCLS" || device == "GENSAL" || device == "GENROU" ||
          device == "WSIEG1" || device == "EXDC1" || device == "EXDC2" ||
          device == "ESST1A" || device == "ESST4B" || device == "GGOV1" ||
          device == "WSHYGP") {
        ret = true;
      }
      return ret;
    }

    // Utility function to check if device is on a bus
    bool onBus(std::string &device) {
      bool ret = false;
      if (device == "LVSHBL" || device == "FRQTPAT") {
        ret = true;
      }
      return ret;
    }

    // Utility function to check if device is on a branch
    bool onBranch(std::string &device) {
      bool ret = false;
      if (device == "DISTR1") {
        ret = true;
      }
      return ret;
    }

    // Utility function to check if parameters describe a load
    bool onLoad(std::string &device) {
      bool ret = false;
      if (device == "CIM6BL" || device == "USRLOD" ||
          device == "IEELBL") {
        ret = true;
      }
      return ret;
    }


    // Extract extension from file name and convert it to lower case
    std::string getExtension(const std::string file)
    {
      std::string ret;
      std::string line = file;
      int ntok1 = line.find('.',0);
      if (ntok1 == std::string::npos) return
        ret;
      ntok1++;
      int ntok2 = line.find(' ',ntok1);
      if (ntok2 == std::string::npos)
        ntok2 = line.size();
      // get extension
      ret = line.substr(ntok1,ntok2-ntok1);
      // convert all characters to lower case 
      int size = ret.size();
      int i;
      for (i=0; i<size; i++) {
        if (isalpha(ret[i])) {
          ret[i] = tolower(ret[i]);
        }
      }
      return ret;
    }

    void find_ds_par(std::ifstream & input)
    {
      std::string          line;
      gridpack::component::DataCollection *data;
      while(std::getline(input,line)) {
        // Check to see if line is blank
        int idx = line.find_first_not_of(' ');
        if (idx == std::string::npos) continue;

        std::string record = line;
        idx = line.find('/');
        while (idx == std::string::npos) {
          std::getline(input,line);
          idx = line.find('/');
          record.append(line);
        }
        idx = record.find('/');
        if (idx != std::string::npos) record.erase(idx,record.length()-idx);
        std::vector<std::string>  split_line;
        boost::split(split_line, record, boost::algorithm::is_any_of(","),
            boost::token_compress_on);

        std::string sval;
        // MODEL TYPE              "MODEL"                  string
        gridpack::utility::StringUtils util;
        sval = util.trimQuotes(split_line[1]);
        util.toUpper(sval);

        if (onGenerator(sval)) {
          // GENERATOR_BUSNUMBER               "I"                   integer
          int l_idx, o_idx;
          o_idx = atoi(split_line[0].c_str());
#ifdef OLD_MAP
          std::map<int, int>::iterator it;
#else
          boost::unordered_map<int, int>::iterator it;
#endif
          it = p_busMap->find(o_idx);
          if (it != p_busMap->end()) {
            l_idx = it->second;
          } else {
            continue;
          }
          data = dynamic_cast<gridpack::component::DataCollection*>
            (p_network->getBusData(l_idx).get());

          // Find out how many generators are already on bus
          int ngen = 0;
          data->getValue(GENERATOR_NUMBER, &ngen);
          // Identify index of generator to which this data applies
          int g_id = -1;
          if (ngen > 0) {
            // Clean up 2 character tag for generator ID
            std::string tag = util.clean2Char(split_line[2]);
            int i;
            for (i=0; i<ngen; i++) {
              std::string t_id;
              data->getValue(GENERATOR_ID,&t_id,i);
              if (tag == t_id) {
                g_id = i;
                break;
              }
            }
          }

          double rval;
          int ival;
          bool bval;

          if (g_id > -1) {
            if (sval == "GENCLS") {
              GenclsParser<gen_params> parser;
              parser.parse(split_line, data, g_id);
            } else if (sval == "GENSAL") {
              GensalParser<gen_params> parser;
              parser.parse(split_line, data, g_id);
            } else if (sval == "GENROU") {
              GenrouParser<gen_params> parser;
              parser.parse(split_line, data, g_id);
            } else if (sval == "WSIEG1") {
              Wsieg1Parser<gen_params> parser;
              parser.parse(split_line, data, g_id);
            } else if (sval == "EXDC1" || sval == "EXDC2") {
              Exdc1Parser<gen_params> parser;
              parser.parse(split_line, data, g_id);
            } else if (sval == "ESST1A") {
              Esst1aParser<gen_params> parser;
              parser.parse(split_line, data, g_id);
            } else if (sval == "ESST4A") {
              Esst4bParser<gen_params> parser;
              parser.parse(split_line, data, g_id);
            } else if (sval == "GGOV1") {
              Ggov1Parser<gen_params> parser;
              parser.parse(split_line, data, g_id);
            } else if (sval == "WSHYGP") {
              WshygpParser<gen_params> parser;
              parser.parse(split_line, data, g_id);
            }
          }
        } else if (onBus(sval)) {
          int l_idx, o_idx;
          if (sval == "LVSHBL") {
            o_idx = atoi(split_line[0].c_str());
          } else if (sval == "FRQTPAT") {
            o_idx = atoi(split_line[3].c_str());
          }
#ifdef OLD_MAP
          std::map<int, int>::iterator it;
#else
          boost::unordered_map<int, int>::iterator it;
#endif
          it = p_busMap->find(o_idx);
          if (it != p_busMap->end()) {
            l_idx = it->second;
          } else {
            continue;
          }
          data = dynamic_cast<gridpack::component::DataCollection*>
            (p_network->getBusData(l_idx).get());
          if (sval == "LVSHBL") {
            LvshblParser<gen_params> parser;
            parser.parse(split_line, data);
          } else if (sval == "FRQTPAT") {
            FrqtpatParser<gen_params> parser;
            parser.parse(split_line, data);
          }
        } else if (onLoad(sval)) {
          // Load bus number
          int l_idx, o_idx;
          o_idx = atoi(split_line[0].c_str());
#ifdef OLD_MAP
          std::map<int, int>::iterator it;
#else
          boost::unordered_map<int, int>::iterator it;
#endif
          it = p_busMap->find(o_idx);
          if (it != p_busMap->end()) {
            l_idx = it->second;
          } else {
            continue;
          }
          data = dynamic_cast<gridpack::component::DataCollection*>
            (p_network->getBusData(l_idx).get());

          // Find out how many generators are already on bus
          int nload = 0;
          data->getValue(LOAD_NUMBER, &nload);
          // Identify index of generator to which this data applies
          int l_id = -1;
          if (nload > 0) {
            // Clean up 2 character tag for generator ID
            std::string tag = util.clean2Char(split_line[2]);
            int i;
            for (i=0; i<nload; i++) {
              std::string t_id;
              data->getValue(LOAD_ID,&t_id,i);
              if (tag == t_id) {
                l_id = i;
                break;
              }
            }
          }
          if (sval == "CIM6BL") {
            Cim6blParser<load_params> parser;
            parser.parse(split_line, data, l_id);
          } else if (sval == "IEELBL") {
            IeelblParser<load_params> parser;
            parser.parse(split_line, data, l_id);
          } else if (sval == "USRLOD") {
            std::string sdev;
            sdev = util.trimQuotes(split_line[3]);
            if (sdev == "ACMTBLU1") {
              Acmtblu1Parser<load_params> parser;
              parser.parse(split_line, data, l_id);
            } else if (sdev == "CMLDBLU1") {
              Cmldblu1Parser<load_params> parser;
              parser.parse(split_line, data, l_id);
            }
          }
        } else if (onBranch(sval)) {
          int l_idx, from_idx, to_idx;
          from_idx = atoi(split_line[0].c_str());
          to_idx = atoi(split_line[2].c_str());
          std::map<std::pair<int, int>, int>::iterator it;
          it = p_branchMap->find(std::pair<int,int>(from_idx,to_idx));
          if (it != p_branchMap->end()) {
            l_idx = it->second;
          } else {
            continue;
          }
          if (sval == "DISTR1") {
            Distr1Parser<gen_params> parser;
            parser.parse(split_line, data);
          }
        }
      }
    }

    // Parse file to construct lists of structs representing different devices.
    void find_ds_vector(std::ifstream & input, std::vector<gen_params> *gen_vector,
        std::vector<bus_relay_params> *bus_relay_vector,
        std::vector<branch_relay_params> *branch_relay_vector,
        std::vector<load_params> *load_vector)
    {
      std::string          line;
      gen_vector->clear();
      while(std::getline(input,line)) {
        // Check to see if line is blank
        int idx = line.find_first_not_of(' ');
        if (idx == std::string::npos) continue;

        std::string record = line;
        idx = line.find('/');
        while (idx == std::string::npos) {
          std::getline(input,line);
          idx = line.find('/');
          record.append(line);
        }
        idx = record.find('/');
        if (idx != std::string::npos) record.erase(idx,record.length()-idx);
        std::vector<std::string>  split_line;
        boost::split(split_line, record, boost::algorithm::is_any_of(","),
            boost::token_compress_on);
        std::string sval;
        gridpack::utility::StringUtils util;
        sval = util.trimQuotes(split_line[1]);
        util.toUpper(sval);

        if (onGenerator(sval)) {
          gen_params data;

          // GENERATOR_BUSNUMBER               "I"                   integer
          int o_idx;
          o_idx = atoi(split_line[0].c_str());
          data.bus_id = o_idx;

          // Clean up 2 character tag for generator ID
          std::string tag = util.clean2Char(split_line[2]);
          strcpy(data.gen_id, tag.c_str());

          double rval;
          int ival;


          // GENERATOR_MODEL              "MODEL"                  integer
          strcpy(data.model, sval.c_str());

          if (sval == "GENCLS") {
            GenclsParser<gen_params> parser;
            parser.store(split_line,data);
          } else if (sval == "GENSAL") {
            GensalParser<gen_params> parser;
            parser.store(split_line,data);
          } else if (sval == "GENROU") {
            GenrouParser<gen_params> parser;
            parser.store(split_line,data);
          } else if (sval == "WSIEG1") {
            Wsieg1Parser<gen_params> parser;
            parser.store(split_line,data);
          } else if (sval == "EXDC1" || sval == "EXDC2") {
            Exdc1Parser<gen_params> parser;
            parser.store(split_line,data);
          } else if (sval == "ESST1A") {
            Esst1aParser<gen_params> parser;
            parser.store(split_line,data);
          } else if (sval == "ESST4B") {
            Esst4bParser<gen_params> parser;
            parser.store(split_line,data);
          } else if (sval == "GGOV1") {
            Ggov1Parser<gen_params> parser;
            parser.store(split_line,data);
          } else if (sval == "WSHYGP") {
            WshygpParser<gen_params> parser;
            parser.store(split_line,data);
          }
          gen_vector->push_back(data);
        } else if (onBus(sval)) {

          // RELAY_BUSNUMBER               "I"                   integer
          int o_idx;
          if (sval == "LVSHBL") {
            bus_relay_params data;
            o_idx = atoi(split_line[0].c_str());
            data.bus_id = o_idx;
            LvshblParser<bus_relay_params> parser;
            parser.store(split_line,data);
            bus_relay_vector->push_back(data);
          } else if (sval == "FRQTPAT") {
            bus_relay_params data;
            o_idx = atoi(split_line[3].c_str());
            data.bus_id = o_idx;
            FrqtpatParser<bus_relay_params> parser;
            parser.store(split_line,data);
            bus_relay_vector->push_back(data);
          }
        } else if (onLoad(sval)) {
          // ID of bus that owns load
          load_params data;
          int o_idx = atoi(split_line[0].c_str());
          data.bus_id = o_idx;

          // Clean up 2 character tag for generator ID
          std::string tag = util.clean2Char(split_line[2]);
          strcpy(data.id, tag.c_str());
          if (sval == "CIM6BL") {
            Cim6blParser<load_params> parser;
            parser.store(split_line,data);
          } else if (sval == "IEELBL") {
            IeelblParser<load_params> parser;
            parser.store(split_line,data);
          } else if (sval == "USRLOD") {
            std::string sdev;
            sdev = util.trimQuotes(split_line[3]);
            if (sdev == "ACMTBLU1") {
              Acmtblu1Parser<load_params> parser;
              parser.store(split_line,data);
            } else if (sdev == "CMLDBLU1") {
              Cmldblu1Parser<load_params> parser;
              parser.store(split_line,data);
            }
          }
          load_vector->push_back(data);
        } else if (onBranch(sval)) {
          branch_relay_params data;

          int from_idx, to_idx;
          if (sval == "DISTR1") {
            from_idx = atoi(split_line[0].c_str());
            to_idx = atoi(split_line[3].c_str());
            data.from_bus = from_idx;
            data.to_bus = to_idx;
            Distr1Parser<branch_relay_params> parser;
            parser.store(split_line,data);
          }
          branch_relay_vector->push_back(data);
        }
      }
    }

    void find_uc_vector(std::ifstream & input, std::vector<uc_params> *uc_vector)
    {
      std::string          line;
      uc_vector->clear();
      // Ignore first line containing header information
      std::getline(input,line);
      while(std::getline(input,line)) {
        std::vector<std::string>  split_line;
        boost::split(split_line, line, boost::algorithm::is_any_of(","),
            boost::token_compress_on);

        uc_params data;

        int nstr = split_line.size();
        if (nstr > 1) {
          data.type = atoi(split_line[1].c_str());
        }
        if (nstr > 2) {
          data.init_level = atof(split_line[2].c_str());
        }
        if (nstr > 3) {
          data.min_gen = atof(split_line[3].c_str());
        }
        if (nstr > 4) {
          data.max_gen = atof(split_line[4].c_str());
        }
        if (nstr > 5) {
          data.max_oper = atof(split_line[5].c_str());
        }
        if (nstr > 6) {
          data.min_up = atof(split_line[6].c_str());
        }
        if (nstr > 7) {
          data.min_down = atof(split_line[7].c_str());
        }
        if (nstr > 8) {
          data.ramp_up = atof(split_line[8].c_str());
        }
        if (nstr > 9) {
          data.ramp_down = atof(split_line[9].c_str());
        }
        if (nstr > 10) {
          data.start_up = atof(split_line[10].c_str());
        }
        if (nstr > 11) {
          data.const_cost = atof(split_line[11].c_str());
        }
        if (nstr > 12) {
          data.lin_cost = atof(split_line[12].c_str());
        }
        if (nstr > 13) {
          data.co_2_cost = atof(split_line[13].c_str());
        }
        if (nstr > 14) {
          data.init_prd = atof(split_line[14].c_str());
        }
        if (nstr > 15) {
          data.start_cap = atof(split_line[15].c_str());
        }
        if (nstr > 16) {
          data.shut_cap = atof(split_line[16].c_str());
        }
        if (nstr > 17) {
          data.bus_id = atoi(split_line[17].c_str());
        }
        if (nstr > 18) {
          // Clean up 2 character tag for generator ID
          gridpack::utility::StringUtils util;
          std::string tag = util.clean2Char(split_line[18]);
          strcpy(data.gen_id, tag.c_str());
        }
        uc_vector->push_back(data);
      }
    }

    /** Store pointer to bus and branch maps
     * @param busMap pointer to map between PTI and local indices
     * @param branchMap pointer to map between PTI bus pairs and local branch
     *        indices
     */
    void setMaps(std::map<int,int> *busMap,
                 std::map<std::pair<int, int>, int> *branchMap)
    {
      p_busMap = busMap;
      p_branchMap = branchMap;
    }

    /**
     * Expand any compound bus models that may need to be generated based on
     * parameters in the .dyr files. This function needs to be called after
     * calling the parser for the .dyr file
     */
    void expandBusModels(void)
    {
      // Find maximum original bus index. This value can be used to assign
      // indices to new buses generated by expanding the model
      int max_idx = -1;
      int nbus = p_network->numBuses();
      int totalBuses = p_network->totalBuses();
      int totalBranches = p_network->totalBranches();
      int i, j;
      for (i=0; i<nbus; i++) {
        if (p_network->activeBus(i)) {
          if (max_idx < p_network->getOriginalBusIndex(i)) {
            max_idx = p_network->getOriginalBusIndex(i);
          }
        }
      }
      p_network->communicator()->max(&max_idx,1);

      // Create new data collection objects corresponding to new buses and
      // branches for the expanded composite load
      std::vector<std::vector<gridpack::component::DataCollection*> >
        new_buses;
      std::vector<std::vector<gridpack::component::DataCollection*> >
        new_branches;
      gridpack::component::DataCollection *data;
      std::string model;
      int nload;
      std::vector<int> comp_buses;
      for (i=0; i<nbus; i++) {
        // Only expand bus on process that owns bus. New buses and branches will
        // not have connections to other processors
        if (p_network->activeBus(i)) {
         data = p_network->getBusData(i).get();
         if (data->getValue("LOAD_NUMBER",&nload)) {
           for (j=0; j<nload; j++) {
             if (data->getValue("LOAD_MODEL",&model, j)) {
               if (model == "CMLDBLU1") {
                 std::vector<gridpack::component::DataCollection*> buses;
                 std::vector<gridpack::component::DataCollection*> branches;
                 Cmldblu1Parser<load_params> parser;
                 parser.expandModel(data, buses, branches, j);
                 new_buses.push_back(buses);
                 new_buses.push_back(branches);
                 comp_buses.push_back(i);
               }
             }
           }
         }
        }
      }

      // Find new indices for the buses. Start by finding out how many new buses
      // have been added to the system
      int nprocs = p_network->communicator()->size();
      int me = p_network->communicator()->rank();
      int added_buses[nprocs];
      for (i=0; i<nprocs; i++) {
        added_buses[i] = 0;
      }
      int numLoads = new_buses.size();
      // Subtract 1 because bus 0 corresponds to the original bus and
      // already has an index
      for (i=0; i<numLoads; i++) {
        added_buses[me] += new_buses[i].size()-1;
      }
      p_network->component()->sum(added_buses,nprocs);
      int offset[nprocs];
      offset[0] = 0;
      for (i=1; i<nprocs; i++) {
        offset[i] = offset[i-1] + added_buses[i-1];
      }
      // Assign new indices to the new buses
      int icnt = 0;
      int idx;
      std::vector<std::vector<int> > local_bus;
      for (i=0; i<numLoads; i++) {
        std::vector<int> lidx;
        lidx.push_back(comp_buses[i]);
        for (j=1; j<new_buses.size(); j++) {
          data = (new_buses[i])[j];
          idx = offset[me]+icnt+max_idx+1;
          data->addValue(BUS_NUMBER, idx);
          p_network->addBus(idx);
          *(p_network->getBusData(nbus+icnt))
            = *data;
          p_network->setActiveBus(nbus+icnt,true);
          lidx.push_back(nbus+icnt);
          icnt++;
        }
        local_bus.push_back(lidx);
      }
      // New buses now have original IDs. Need to set "from" and "to"
      // indices for new branches as well as fixing up neighor lists
      int from, to, i1, i2;
      icnt = 0;
      int nbranch = p_network->numBranches();
      for (i=0; i<numLoads; i++) {
        for (j=0; j<new_branches.size(); j++) {
          // Assign new bus IDs to new branches
          data = (new_branches[i])[j];
          data->getValue(BRANCH_FROMBUS,&from);
          data->getValue(BRANCH_TOBUS,&to);
          (new_buses[i])[from]->getValue(BUS_NUMBER, &i1);
          (new_buses[i])[to]->getValue(BUS_NUMBER, &i2);
          data->setValue(BRANCH_FROMBUS,i1);
          data->setValue(BRANCH_TOBUS,i2);
          p_network->addBranch(i1,i2);
          *(p_network->getBranchData(nbranch+icnt))
            = *data;
          p_network->setActiveBranch(nbranch+icnt,true);
          // Modify neighbor lists etc to account for new buses
          p_network->setLocalBusIndex1(nbranch+icnt,(local_bus[i])[from]);
          p_network->setLocalBusIndex2(nbranch+icnt,(local_bus[i])[to]);
          p_network->addNeighbor((local_bus[i])[from],nbranch+icnt);
          p_network->addNeighbor((local_bus[i])[to],nbranch+icnt);
          icnt++;
        }
      }
      // Assign global indices to new buses.
      icnt = 0;
      int lcnt = 0;
      // loop over original buses
      for (i=0; i<nbus; i++) {
        if (p_network->getActiveBus(i)) {
          if (i == comp_buses[icnt]) {
            // if original bus is composite bus, loop over new buses;
            for (j=1; j<local_bus[i].size(); j++) {
              p_network->setGlobalBusIndex((local_bus[i])[j],
                  totalBuses+offset[me]+lcnt);
              lcnt++;
            }
            icnt++;
          }
        }
      }

      // Add global indices to new local branches.
      int nbranch_new = p_network->numBranches();
      int branch_buf[nprocs];
      for (i=0; i<nprocs; i++) {
        branch_buf[i] = 0;
      }
      branch_buf[me] = nbranch_new-nbranch;
      p_network->component()->sum(branch_buf,nprocs);
      offset[0] = 0;
      for (i=1; i<nprocs; i++) {
        offset[i] = offset[i-1]+branch_buf[i-1];
      }
      icnt = offset[me];
      for (i=nbranch; i<nbranch_new; i++) {
        if (p_network->getActiveBranch(i)) {
          p_network->setGlobalBranchIndex(i,totalBranches+icnt);
          icnt++;
        }
      }

      // Reset remaining indices
      p_network->resetGlobalIndices(false);

    }

    /* ************************************************************************
     **************************************************************************
     ***** OBJECT DATA
     **************************************************************************
     *********************************************************************** */
    boost::shared_ptr<_network> p_network;

    gridpack::utility::CoarseTimer *p_timer;

    // Map of PTI indices to index in p_busData
    std::map<int,int> *p_busMap;
    // Map of PTI index pair to index in p_branchData
    std::map<std::pair<int, int>, int> *p_branchMap;

}; /* end of PTI base parser */

} /* namespace parser */
} /* namespace gridpack */

#endif /* BASEPTIPARSER_HPP_ */
