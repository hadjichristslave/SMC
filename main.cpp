#include <iostream>
#include <stdlib.h>
#include <libconfig.h>
#include "Utilities.h"
#include "Landmarks.h"
#include "DBWrapper.h"
#include <map>
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
    //cout << "Config file read in "  << end-start << " s" <<endl;
    //cout << "Starting, parameter init" << endl;
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
    //cout << "Everything initialized in " << end-start << " s " << endl;

    start = omp_get_wtime();
    smc.infer( &particles, & dataPoints, baseparams, numOfParticles , numOfSamples);
    cout << " Gost out in  " << omp_get_wtime() - start<<endl;
    Landmarks observations;
    for(unsigned int i = 0;i< particles[0].stateProg[0].size();i++){
            //cout <<particles[0].stateProg[0][i] <<endl;
            Landmark land(smc.numOfLandmarks , particles[0].stateProg[0][i],i);
            observations.addLandMark(land);
    }
    //cout << endl<< endl<< endl<< endl<< endl<< " assignemnts" <<endl;
    for(unsigned int i = 0;i< particles[0].assignments.size();i++)
            cout  << dataPoints[0][i][0] << " " << dataPoints[0][i][1] << " " << dataPoints[0][i][2] << " " << particles[0].assignments[i] <<endl;
    cout << endl;

    Landmarks                   landmarks   =  dbwr.getCurrentLandmarks();


    //sort( particles[0].assignments.begin(), particles[0].assignments.end() );
    //particles[0].assignments.erase( unique( particles[0].assignments.begin(), particles[0].assignments.end() ), particles[0].assignments.end() );
    //for(unsigned int p = 0;p<particles[0].assignments.size();p++)
    //    cout<< particles[0].assignments[p] << ",";
    //cout<<endl;
    //for(unsigned int i=0;i<landmarks.size();i++){
        //vector< vector< double > > distanceFeatures  = landmarks.extractDistances(& landmarks.landmarks[i],  & ut );
        //cout <<"distance fieatureas size" << distanceFeatures.size() << endl;
        //dbwr.insertLabeledDistances(distanceFeatures, i);
    //}
    //return 0;


    //cout << "Landmarks and db retireved and trained in "  << end-start << " s" << endl;

    int initialDbSize   = landmarks.size();
    vector<double> current_observations;

    //cout << "Clasifying observations to landmarks;"  << endl;
    std::map<int, int >  landmark_associations;


    start = omp_get_wtime();
    for(unsigned int i=0;i<observations.size();i++){
        if(landmarks.size()==0){
            dbwr.insertLandmark(& observations.landmarks[i].distribution);
            landmarks  =   dbwr.getCurrentLandmarks();
            current_observations.push_back(0);
            continue;
        }
        // For every landmark calculate its distances with stored landmarks
        vector< vector< double > > distanceFeatures  = landmarks.extractDistances(& observations.landmarks[i],  & ut,0 , 40 );
        // Get the probability of being the same instance as that given landmark
        bool found = false;
        for(int ijk = 0 ; ijk<distanceFeatures.size();ijk++){
            //cout << distanceFeatures[ijk][0] <<  " "  <<distanceFeatures[ijk][1] << endl;
            if( distanceFeatures[ijk][0]<.3 && distanceFeatures[ijk][1]<20 && found ==false){
                //cout << " dding known landmark for cluster" << i <<endl;
                landmark_associations.insert(std::make_pair(i, landmarks.landmarks[ijk].uuid));
                //Landmark is registered as currently detected
                current_observations.push_back(ijk);
                found  =true;
                break;
            }
        }
        cout << "gotoutp" <<endl;
        if(found==false){
                //cout << " adding a new landmark" <<endl;
                dbwr.insertLandmark(& observations.landmarks[i].distribution);
                current_observations.push_back(observations.size());
                landmarks =  dbwr.getCurrentLandmarks();
                landmark_associations.insert(std::make_pair(current_observations.back(), landmarks.landmarks.back().uuid));
        }
    }
    cout << " landmark matching in   " << omp_get_wtime()-start<<endl;
    ofstream myfile;
    myfile.open ("/home/panos/Desktop/aggregated.txt",std::ofstream::out | std::ofstream::app);
    myfile << "[";


    for(int i =0 ;i< dataPoints[0].size() ;i++){
        //myfile <<  dataPoints[0][i][0] << " " << dataPoints[0][i][1] << " " << dataPoints[0][i][2] << " " << landmark_associations[particles[0].assignments[i]]<< endl;
    }
    myfile<< "]";
    myfile.close();
    // if initial db size is zero create an initial training sample for the random forest
    // For now, I assume the training db was populated
    landmarks =  dbwr.getCurrentLandmarks();
    for(int i =0;i<landmarks.size();i++){
        //if(landmarks.landmarks[i].uuid>2480) break;
        vector< vector< double > > distanceFeatures = landmarks.extractDistances(& landmarks.landmarks[i],  & ut,0 , 40 );
        //cout << "-------------- " <<endl;
        //cout << " Distances for landmark with id"  << landmarks.landmarks[i].uuid <<endl;
        for(int j = 0 ; i < distanceFeatures.size(); j ++){
            if(landmarks.landmarks[j].uuid<1 ||landmarks.landmarks[j].uuid>10000  ) break;
            //cout << landmarks.landmarks[j].uuid  << ":" << distanceFeatures[j][0] << "," << distanceFeatures[j][1] << ","<< \
            distanceFeatures[j][2] << ","<< distanceFeatures[j][3] << ","<< distanceFeatures[j][4] << ","<<endl;
        }
    }
    return 0;
}
