#include <iostream>
#include <stdlib.h>
#include <libconfig.h>
#include "Utilities.h"
#include "Landmarks.h"
#include "DBWrapper.h"
#include <omp.h>

using namespace std;
using namespace Structures;

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
int main(int argc, char* argv[]){

    cout << "Starting, reading config" << endl;
    double start = omp_get_wtime();

    // Config file parsing
    int numOfParticles , numOfSamples;
    double landmarkThreshold;
    const char *filepath = NULL;
    const char *database = NULL;
    config_t cfg, *cf;
    cf = &cfg;
    config_init(cf);
    if (!config_read_file(cf, "config.cfg"))
        return(EXIT_FAILURE);
    config_lookup_int(cf, "particles" , &numOfParticles);
    config_lookup_int(cf, "samples" , &numOfSamples);
    config_lookup_float(cf, "landmarkThreshold" , &landmarkThreshold);
    config_lookup_string(cf, "filepath" , &filepath);
    config_lookup_string(cf, "database" , &database);

    // Variable declaration
    //The database wrapper object

    double end = omp_get_wtime();
    cout << "Config file read in "  << end-start << " s" <<endl;
    cout << "Starting, parameter init" << endl;
    start = omp_get_wtime();


    DBWrapper dbwr(database);
    // The smc sampler object
    SMC smc;
    // Utilities object. Everything that does not fit something else just goes there
    Utilities ut;
    //Data points, our cloud information goes here
    vector< vector< vector< double > > > dataPoints  = ut.readFile("-----", filepath);
    // All different points in time of our pointclouds.
    int timeStates = dataPoints.size();
    // Define the base params size
    Params baseparams;
    baseparams.cloudInstances = dataPoints.size();
    // Create as many particles as config file defines
    vector < StateProgression > particles(numOfParticles, timeStates);
    // Initialize the method
    smc.init();
    // Cluster the elements

    end = omp_get_wtime();
    cout << "Everything initialized in " << end-start << " s " << endl;

    cout << "Clustering starting"  << endl;
    start = omp_get_wtime();

    smc.infer( &particles, & dataPoints, baseparams, numOfParticles , numOfSamples);

    end= omp_get_wtime();
    cout << " clustering finished in  " << end-start <<  " s" << endl;
    //Get the clusters output. Single particle or a mixture can be used. Only particle 1 is used here

    Landmarks observations;
    for(unsigned int i = 0;i< particles[0].stateProg[0].size();i++){
            Landmark land(smc.numOfLandmarks , particles[0].stateProg[0][i]);
            observations.addLandMark(land);
    }
    // Get landmarks currently in the database
    cout << "Retrieving landmarks and trainign set"  << endl;
    start = omp_get_wtime();

    Landmarks                   landmarks   =  dbwr.getCurrentLandmarks();
    vector< vector< double > >  trainingSet =  dbwr.getTrainingSet();
    ut.decisionTrain(& trainingSet);

    //for(unsigned int i=0;i<landmarks.size();i++){
        //vector< vector< double > > distanceFeatures  = landmarks.extractDistances(& landmarks.landmarks[i],  & ut );
        //dbwr.insertLabeledDistances(distanceFeatures, i);
    //}
    //return 0;

    end = omp_get_wtime();
    cout << "Landmarks and db retireved and trained in "  << end-start << " s" << endl;

    int initialDbSize   = landmarks.size();
    vector<double> current_observations;

    cout << "Clasifying observations to landmarks;"  << endl;
    start = omp_get_wtime();

    for(unsigned int i=0;i<observations.size();i++){

        if(landmarks.size()==0){
            dbwr.insertLandmark(& observations.landmarks[i].distribution);
            landmarks  =   dbwr.getCurrentLandmarks();
            current_observations.push_back(landmarks.landmarks.size()-1);
            continue;
        }
        // For every landmark calculate its distances with stored landmarks
        vector< vector< double > > distanceFeatures  = landmarks.extractDistances(& observations.landmarks[i],  & ut );
        // Get the probability of being the same instance as that given landmark

        vector<double> probabilities = ut.observationProbabilities(& distanceFeatures);

        //Normalized quantities reduce the confidene interval < than .2 due to the large number of landmarks in the training set.
        //ut.normalizeVec(&probabilities);
        for(auto ij: sort_indexes(probabilities)){
            cout << "Top probability" << probabilities[ij]<< endl;
            //if(probabilities[ij]<.1)
                //cout << observations.landmarks[i].distribution;

            //if probability is larger than .9 then we have a matchs
            if( probabilities[ij]>landmarkThreshold){
                //Landmark is registered as currently detected
                current_observations.push_back(ij);
                break;
            }else{
                //for_each(probabilities.begin(),probabilities.end(),[](double y){ cout << y <<","; });
                //cout << endl;
                //Insert the landmark, update landmark db, add new landmark id to the currently detected list
                dbwr.insertLandmark(& observations.landmarks[i].distribution);
                current_observations.push_back(observations.size()-1);
                landmarks =  dbwr.getCurrentLandmarks();
                break;
            }
        }
        //;
    }
    end = omp_get_wtime();
    cout << "Landmarks classified in "  << end-start <<" s"<<  endl;

    // if initial db size is zero create an initial training sample for the random forest
    // For now, I assume the training db was populated
    /*if(initialDbSize==0){
        landmarks =  dbwr.getCurrentLandmarks();
        for(unsigned int i=0;i<landmarks.size();i++){
            vector< vector< double > > distanceFeatures  = landmarks.extractDistances(& landmarks.landmarks[i],  & ut );
            dbwr.insertLabeledDistances(distanceFeatures, i);
        }
    }*/
    return 0;
}
