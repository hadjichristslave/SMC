#include <iostream>
#include <string>
#include <stdlib.h>
#include <stdio.h>
#include <sstream>
#include <libconfig.h>
#include <algorithm>
#include <vector>
#include "SMC.h"
#include "Utilities.h"

#include <stdexcept>
#include <vector>
#include <random>
#include "Landmark.h"
#include "Landmarks.h"

#include "DBWrapper.h"
using namespace std;


template <typename T>
vector<size_t> sort_indexes(const vector<T> &v) {

  // initialize original index locations
  vector<size_t> idx(v.size());
  for (size_t i = 0; i != idx.size(); ++i) idx[i] = i;

  // sort indexes based on comparing values in v
  sort(idx.begin(), idx.end(),
       [&v](size_t i1, size_t i2) {return v[i1] > v[i2];});

  return idx;
}
int main(int argc, char* argv[])
{

    // Config file parsing

    int numOfParticles , numOfSamples;
    const char *filepath = NULL;
    const char *database = NULL;
    config_t cfg, *cf;
    cf = &cfg;
    config_init(cf);
    if (!config_read_file(cf, "config.cfg"))
        return(EXIT_FAILURE);
    config_lookup_int(cf, "particles" , &numOfParticles);
    config_lookup_int(cf, "samples" , &numOfSamples);
    config_lookup_string(cf, "filepath" , &filepath);
    config_lookup_string(cf, "database" , &database);

    // Variable declaration
    DBWrapper dbwr(database);
    SMC smc;
    Utilities ut;
    vector< vector< vector< double > > > dataPoints  = ut.readFile("-----", filepath);
    int timeStates = dataPoints.size(); // All different points in time of our pointclouds.
    SMC::Params Baseparams;
    Baseparams.cloudInstances = dataPoints.size();
    vector < SMC::StateProgression > particles(numOfParticles, timeStates);
    // Initialize the method
    smc.init();
    //Get the clusters
    smc.infer( &particles, & dataPoints, Baseparams, numOfParticles , numOfSamples);
    //Compare and decide on whether you have new landmarks or not.
    SMC::StateProgression temp = particles[0];
    Landmarks observations;
    for(unsigned int i = 0;i< temp.stateProg[0].size();i++){
        Landmark land(smc.numOfLandmarks , temp.stateProg[0][i]);
        smc.numOfLandmarks++;
        observations.addLandMark(land);
    }
    Landmarks landmarks =  dbwr.getCurrentLandmarks();
    //Save landmarks to db
    vector<double> current_observations;
    for(unsigned int i=1;i<observations.size();i++){
        vector< vector< double > > distanceFeatures  = landmarks.extractDistances(& observations.landmarks[i],  & ut );
        vector<double> probabilities = ut.observationProbabilities(& distanceFeatures);
        for(auto i: sort_indexes(probabilities)){
            //if probability is larger than .9 then we have a match
            if( probabilities[i]>.9){
                //do stuff so that the landmark is registered
                current_observations.push_back(i);
                break;
            }else{
                dbwr.insertLandmark(& observations.landmarks[i].distribution);
                current_observations.push_back(i);
                landmarks =  dbwr.getCurrentLandmarks();
                break;
            }
        }
        //;

    }


    return 0;
}
