#include <iostream>
#include <string>
#include <stdlib.h>
#include <stdio.h>
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

using namespace std;

void extractDistances(Landmark landmark1 , Landmark landmark2 , Utilities ut){
        cout << "wasser " << endl;
        cout << ut.Wasserstein(landmark1.distribution.mean, landmark1.distribution.covar, landmark2.distribution.mean, landmark2.distribution.covar) << endl;
        cout << "gauskld" << endl;
        cout << ut.GaussKLDivergence(landmark1.distribution.mean, landmark1.distribution.covar, landmark2.distribution.mean, landmark2.distribution.covar) << endl;
        cout << "expkld " << endl;
        cout << ut.ExpKLDivergence(landmark1.distribution.exponential, landmark2.distribution.exponential) << endl;
        cout << " exphellinger " << endl;
        cout << ut.Expsquaredhellinger(landmark1.distribution.exponential, landmark2.distribution.exponential) << endl;
        float histogram1[ landmark1.distribution.categorical.size() ];
        float histogram2[ landmark2.distribution.categorical.size() ];
        for (unsigned int o=0;o<landmark2.distribution.categorical.size();o++){
            histogram1[o] = landmark1.distribution.categorical[o]==0?.000001:landmark1.distribution.categorical[o];
            histogram2[o] = landmark2.distribution.categorical[o]==0?.000001:landmark1.distribution.categorical[o];
        }
        vector<float> distances = ut.categoricalhistogramCompare(histogram1, histogram2, sizeof(histogram1)/sizeof(float));
        for_each(distances.begin(), distances.end(), [] (float y) { cout << y << ",";});
        cout << endl;
}

int main(int argc, char* argv[])
{

    // Config file parsing

    int numOfParticles , numOfSamples;
    const char *filepath = NULL;
    config_t cfg, *cf;
    cf = &cfg;
    config_init(cf);
    if (!config_read_file(cf, "config.cfg"))
        return(EXIT_FAILURE);
    config_lookup_int(cf, "particles" , &numOfParticles);
    config_lookup_int(cf, "samples" , &numOfSamples);
    config_lookup_string(cf, "filepath" , &filepath);
    // Variable declaration

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
    for(unsigned int i=1;i<lands.size();i++){
        cout << " new distance pair " << endl;
        extractDistances(lands.landmarks[i], lands.landmarks[i-1], ut);
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
