#include <iostream>
#include <stdlib.h>
#include <libconfig.h>
#include "Utilities.h"
#include "Landmarks.h"
#include "DBWrapper.h"


using namespace std;
using namespace Structures;
typedef vector< vector< vector< double > > > dataBuffer;

int main(int argc, char* argv[]){

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

    DBWrapper dbwr(database);
    // The Sequential monte carlo sampler object
    SMC smc;
    Utilities ut;
    dataBuffer dataPoints  = ut.readFile("-----", filepath);

    int timeStates = dataPoints.size();
    // Define the base params size

    Params baseparams;
    // Create as many particles as config file defines
    vector < StateProgression > particles(numOfParticles, timeStates);
    // Initialize the method
    smc.init();
    // Cluster the elements

    smc.infer( &particles, & dataPoints, baseparams, numOfParticles , numOfSamples);
    Landmarks observations;
    for(unsigned int i = 0;i< particles[0].stateProg[0].size();i++){
            Landmark land(smc.numOfLandmarks , particles[0].stateProg[0][i],i);
            observations.addLandMark(land);
    }
    //for(unsigned int i = 0;i< particles[0].assignments.size();i++)
      //      cout  << dataPoints[0][i][0] << " " << dataPoints[0][i][1] << " " << dataPoints[0][i][2] << " " << particles[0].assignments[i] <<endl;
    //cout << endl;

    Landmarks landmarks   =  dbwr.getCurrentLandmarks();
    vector<double> current_observations;
    std::map<int, int >  landmark_associations;

    for(unsigned int i=0;i<observations.size();i++){
        if(landmarks.size()==0){
            dbwr.insertLandmark(& observations.landmarks[i].distribution);
            landmarks  =   dbwr.getCurrentLandmarks();
            current_observations.push_back(0);
            continue;
        }
        vector< vector< double > > distanceFeatures  = landmarks.extractDistances(& observations.landmarks[i],  & ut);
        bool found = false;
        for(unsigned int ijk = 0 ; ijk<distanceFeatures.size();ijk++){
            //cout << distanceFeatures[ijk][0] <<  " "  <<distanceFeatures[ijk][1] << endl;
            if( distanceFeatures[ijk][0]<.3 && distanceFeatures[ijk][1]<20 && found ==false){
                //cout << " dding known landmark for cluster" << i <<endl;
                landmark_associations.insert(std::make_pair(i, landmarks.landmarks[ijk].getId()));
                //Landmark is registered as currently detected
                current_observations.push_back(ijk);
                found  =true;
                break;
            }
        }
        if(found==false){
                dbwr.insertLandmark(& observations.landmarks[i].distribution);
                current_observations.push_back(observations.size());
                landmarks = dbwr.getCurrentLandmarks();
                landmark_associations.insert(std::make_pair(current_observations.back(), landmarks.landmarks.back().getId()));
        }
    }
    ofstream myfile;
    myfile.open ("/home/panos/Desktop/aggregated.txt",std::ofstream::out | std::ofstream::app);
    myfile << "[";
    for(unsigned int i =0 ;i< dataPoints[0].size() ;i++)
        myfile <<  dataPoints[0][i][0] << " " << dataPoints[0][i][1] << " " << dataPoints[0][i][2] << " " << landmark_associations[particles[0].assignments[i]]<< endl;
    myfile<< "]";
    myfile.close();
    return 0;
}
