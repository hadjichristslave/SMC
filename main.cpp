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
#include "decision-trees.hxx"
#include "marray.hxx"
#include "DBWrapper.h"
using namespace std;

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
    smc.infer( &particles, dataPoints, Baseparams, numOfParticles , numOfSamples);
    //Compare and decide on whether you have new landmarks or not.
    SMC::StateProgression temp = particles[0];
    Landmarks lands;
    for(unsigned int i = 0;i< temp.stateProg[0].size();i++){
        Landmark land(smc.numOfLandmarks , temp.stateProg[0][i]);
        smc.numOfLandmarks++;
        lands.addLandMark(land);
    }
    //Save landmarks to db
    for(unsigned int i=1;i<lands.size();i++){
        //cout << " new distance pair " << endl;
        SMC::SufficientStatistics dist = lands.landmarks[i].distribution;
        dbwr.insertLandmark(dist);
        //extractDistances(lands.landmarks[i], lands.landmarks[i-1], ut);
    }




/*    const size_t numberOfSamples = 100;
    const size_t numberOfFeatures = 2;

    // define random feature matrix
    std::default_random_engine RandomNumberGenerator;
    typedef double Feature;
    std::uniform_real_distribution<double> randomDistribution(0.0, 1.0);
    const size_t shape[] = {numberOfSamples, numberOfFeatures};
    andres::Marray<Feature> features(shape, shape + 2);
    for(size_t sample = 0; sample < numberOfSamples; ++sample)
    for(size_t feature = 0; feature < numberOfFeatures; ++feature) {
        features(sample, feature) = randomDistribution(RandomNumberGenerator);
    }

    // define labels
    typedef unsigned char Label;
    andres::Marray<Label> labels(shape, shape + 1);
    for(size_t sample = 0; sample < numberOfSamples; ++sample) {
        if((features(sample, 0) <= 0.5 && features(sample, 1) <= 0.5)
        || (features(sample, 0) > 0.5 && features(sample, 1) > 0.5)) {
            labels(sample) = 0;
        }
        else {
            labels(sample) = 1;
        }
    }

    // learn decision forest
    typedef double Probability;
    andres::ml::DecisionForest<Feature, Label, Probability> decisionForest;
    const size_t numberOfDecisionTrees = 10;
    decisionForest.learn(features, labels, numberOfDecisionTrees);

    // predict probabilities for every label and every training sample
    andres::Marray<Probability> probabilities(shape, shape + 2);
    decisionForest.predict(features, probabilities);
    // TODO: test formally
*/
    return 0;
}
