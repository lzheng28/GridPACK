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
#include "gridpack/network/base_network.hpp"
#include "gridpack/parser/base_parser.hpp"
#include "gridpack/parser/hash_distr.hpp"

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

    struct ds_params{
      int bus_id; // ID of bus that owns generator
      char gen_id[3]; // Generator ID
      char gen_model[8];  // Generator model
      double inertia;  // Inertia constant 0
      double damping;  // Damping coefficient
      double reactance; // Transient reactance
    } ;

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

      std::vector<ds_params> ds_data;
      if (me == 0) {
        std::ifstream            input;
        input.open(fileName.c_str());
        if (!input.is_open()) {
          // p_timer->stop(t_ds);
          return;
        }
        find_ds_vector(input, &ds_data);
        input.close();
      }
      int nsize = ds_data.size();
      std::vector<int> buses;
      int i;
      for (i=0; i<nsize; i++) {
        buses.push_back(ds_data[i].bus_id);
      }
      gridpack::hash_distr::HashDistribution<_network,ds_params,ds_params>
        distr(p_network);
      distr.distributeBusValues(buses,ds_data);
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
        std::string tag = ds_data[i].gen_id;
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

        std::string sval;
        double rval;
        // GENERATOR_MODEL              "MODEL"        string
        if (!data->getValue(GENERATOR_MODEL,&sval,g_id)) {
          data->addValue(GENERATOR_MODEL, ds_data[i].gen_model, g_id);
        } else {
          data->setValue(GENERATOR_MODEL, ds_data[i].gen_model, g_id);
        }

        // GENERATOR_INERTIA_CONSTANT_H                float
        if (!data->getValue(GENERATOR_INERTIA_CONSTANT_H,&rval,g_id)) {
          data->addValue(GENERATOR_INERTIA_CONSTANT_H,
              ds_data[i].inertia, g_id);
        } else {
          data->setValue(GENERATOR_INERTIA_CONSTANT_H,
              ds_data[i].inertia, g_id);
        }

        // GENERATOR_DAMPING_COEFFICIENT_0             float
        if (!data->getValue(GENERATOR_DAMPING_COEFFICIENT_0,&rval,g_id)) {
          data->addValue(GENERATOR_DAMPING_COEFFICIENT_0,
              ds_data[i].damping, g_id);
        } else {
          data->setValue(GENERATOR_DAMPING_COEFFICIENT_0,
              ds_data[i].damping, g_id);
        }

        // GENERATOR_TRANSIENT_REACTANCE               float
        if (!data->getValue(GENERATOR_TRANSIENT_REACTANCE,&rval,g_id)) {
          data->addValue(GENERATOR_TRANSIENT_REACTANCE,
              ds_data[i].reactance, g_id);
        } else {
          data->setValue(GENERATOR_TRANSIENT_REACTANCE,
              ds_data[i].reactance, g_id);
        }
      }
      //      p_timer->stop(t_ds);
    }

    // Clean up 2 character tags so that single quotes are removed and single
    // character tags are right-justified. These tags can be delimited by a
    // pair of single quotes, a pair of double quotes, or no quotes
    std::string clean2Char(std::string string)
    {
      std::string tag = string;
      // Find and remove single or double quotes
      int ntok1 = tag.find('\'',0);
      bool sngl_qt = true;
      bool no_qt = false;
      // if no single quote found, then assume double quote or no quote
      if (ntok1 == std::string::npos) {
        ntok1 = tag.find('\"',0);
        // if no double quote found then assume no quote
        if (ntok1 == std::string::npos) {
          ntok1 = tag.find_first_not_of(' ',0);
          no_qt = true;
        } else {
          sngl_qt = false;
        }
      }
      int ntok2;
      if (sngl_qt) {
        ntok1 = tag.find_first_not_of('\'',ntok1);
        ntok2 = tag.find('\'',ntok1);
      } else if (no_qt) {
        ntok2 = tag.find(' ',ntok1);
      } else {
        ntok1 = tag.find_first_not_of('\"',ntok1);
        ntok2 = tag.find('\"',ntok1);
      }
      if (ntok2 == std::string::npos) ntok2 = tag.length();
      std::string clean_tag = tag.substr(ntok1,ntok2-ntok1);
      //get rid of white space
      ntok1 = clean_tag.find_first_not_of(' ',0);
      ntok2 = clean_tag.find(' ',ntok1);
      if (ntok2 == std::string::npos) ntok2 = clean_tag.length();
      tag = clean_tag.substr(ntok1,ntok2-ntok1);
      if (tag.length() == 1) {
        clean_tag = " ";
        clean_tag.append(tag);
      } else {
        clean_tag = tag;
      }
      return clean_tag;
    }

    // Tokenize a string on blanks, but ignore blanks within a text string
    // delimited by single quotes
    std::vector<std::string> blankTokenizer(std::string input)
    {
      std::vector<std::string> ret;
      std::string line = input;
      int ntok1 = line.find_first_not_of(' ',0);
      int ntok2 = ntok1;
      while (ntok1 != std::string::npos) {
        bool quote = false;
        if (line[ntok1] != '\'') {
          ntok2 = line.find(' ',ntok1);
        } else {
          bool quote = true;
          ntok2 = line.find('\'',ntok1);
        }
        if (ntok2 == std::string::npos) ntok2 = line.length();
        if (quote) {
          if (line[ntok2-1] == '\'') {
            ret.push_back(line.substr(ntok1,ntok2-ntok1));
          }
        } else {
          ret.push_back(line.substr(ntok1,ntok2-ntok1));
        }
        ntok1 = line.find_first_not_of(' ',0);
        ntok2 = ntok1;
      }
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
        std::string record = line;
        int idx = line.find('/');
        while (idx == std::string::npos) {
          std::getline(input,line);
          idx = line.find('/');
          record.append(line);
        }
        idx = record.find('/');
        if (idx != std::string::npos) record.erase(idx,record.length()-idx);
        std::vector<std::string>  split_line;
        boost::split(split_line, record, boost::algorithm::is_any_of(","), boost::token_compress_on);

        // GENERATOR_BUSNUMBER               "I"                   integer
        int l_idx, o_idx;
        o_idx = atoi(split_line[0].c_str());
#ifdef OLD_MAP
        std::map<int, int>::iterator it;
#else
        boost::unordered_map<int, int>::iterator it;
#endif
        int nstr = split_line.size();
        it = p_busMap->find(o_idx);
        if (it != p_busMap->end()) {
          l_idx = it->second;
        } else {
          continue;
        }
        data = dynamic_cast<gridpack::component::DataCollection*>
          (p_network->getBusData(l_idx).get());

        // Find out how many generators are already on bus
        int ngen;
        if (!data->getValue(GENERATOR_NUMBER, &ngen)) continue;
        // Identify index of generator to which this data applies
        int g_id = -1;
        // Clean up 2 character tag for generator ID
        std::string tag = clean2Char(split_line[2]);
        int i;
        for (i=0; i<ngen; i++) {
          std::string t_id;
          data->getValue(GENERATOR_ID,&t_id,i);
          if (tag == t_id) {
            g_id = i;
            break;
          }
        }
        if (g_id == -1) continue;

        std::string sval;
        double rval;
        // GENERATOR_MODEL              "MODEL"                  string
        if (!data->getValue(GENERATOR_MODEL,&sval,g_id)) {
          data->addValue(GENERATOR_MODEL, split_line[1].c_str(), g_id);
        } else {
          data->setValue(GENERATOR_MODEL, split_line[1].c_str(), g_id);
        }

        // GENERATOR_INERTIA_CONSTANT_H                           float
        if (nstr > 3) {
          if (!data->getValue(GENERATOR_INERTIA_CONSTANT_H,&rval,g_id)) {
            data->addValue(GENERATOR_INERTIA_CONSTANT_H,
                atof(split_line[3].c_str()), g_id);
          } else {
            data->setValue(GENERATOR_INERTIA_CONSTANT_H,
                atof(split_line[3].c_str()), g_id);
          }
        } 

        // GENERATOR_DAMPING_COEFFICIENT_0                           float
        if (nstr > 4) {
          if (!data->getValue(GENERATOR_DAMPING_COEFFICIENT_0,&rval,g_id)) {
            data->addValue(GENERATOR_DAMPING_COEFFICIENT_0,
                atof(split_line[4].c_str()), g_id);
          } else {
            data->setValue(GENERATOR_DAMPING_COEFFICIENT_0,
                atof(split_line[4].c_str()), g_id);
          }
        }
      }
    }

    void find_ds_vector(std::ifstream & input, std::vector<ds_params> *ds_vector)
    {
      std::string          line;
      ds_vector->clear();
      while(std::getline(input,line)) {
        std::string record = line;
        int idx = line.find('/');
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

        ds_params data;

        // GENERATOR_BUSNUMBER               "I"                   integer
        int o_idx;
        o_idx = atoi(split_line[0].c_str());
        data.bus_id = o_idx;

        int nstr = split_line.size();

        // Clean up 2 character tag for generator ID
        std::string tag = clean2Char(split_line[2]);
        strcpy(data.gen_id, tag.c_str());

        std::string sval;
        double rval;
        // GENERATOR_MODEL              "MODEL"                  integer
        strcpy(data.gen_model, split_line[1].c_str());

        // GENERATOR_INERTIA_CONSTANT_H                           float
        if (nstr > 3) {
          data.inertia = atof(split_line[3].c_str());
        } 

        // GENERATOR_DAMPING_COEFFICIENT_0                           float
        if (nstr > 4) {
          data.damping = atof(split_line[4].c_str());
        }

        // GENERATOR_TRANSIENT_REACTANCE                             float
        if (nstr > 5) {
          data.reactance = atof(split_line[5].c_str());
        }
        ds_vector->push_back(data);
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
