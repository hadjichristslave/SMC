#include "SMC.h"

SMC::SMC()
{
    //ctor

}

SMC::~SMC()
{
    //dtor
}

const void SMC::infer(vector< State > particles, \
                      vector < vector < vector<double> > > cloudData, \
                      SMC::Params params, int numOfParticles,\
                      int numOfSamples){
//for all time intervals

    for(int t=0;t<cloudData.size();t++){
    //for all particles
        for(int j=0;j< numOfParticles; j++){
            SMC::State currState = particles[j];
            //for all samples in time t
            for(int s=0;s< numOfSamples;s++)
                smc_sample(currState, cloudData[t], params, t,  s , cloudData[t].size());
        }
    }


}
// Sample  cluster parameters for cloud data at a given time
const void SMC::smc_sample(SMC::State currState, \
                vector< vector<double> >cloudInstance, \
                SMC::Params params, \
                int currentTime, \
                int currentSample , \
                int dataSize){


    //Get the current clusters.
    int currentClusters = currState.assignments.size();

    // Update cluster parameterse given the points within them

    for(int i=0;i<currentClusters;i++){

    }


    //sample an assignment for every datapoint in the dataset
    for(int i=0;i<cloudInstance.size();i++){
        vector<double> pointInstance = cloudInstance[i];
        vector<double> clusterLogProb;


        std::vector<double> sizes;
        if(currState.clusterSizes.size() >0){
            // take the sizes of the cluster for every item in the dataset.
            //log_prob_obs(k) = log()
            sizes = currState.clusterSizes;
        }
        sizes.push_back(params.crp);
        double sum = 0;
        for_each(sizes.begin(), sizes.end(), [&sum] (double y) mutable { sum +=y; });
        for_each(sizes.begin(), sizes.end(), [&sum] (double &y) mutable { y = y/sum; });


        for( int j =0 ; j<sizes.size(); j++){
             if(sizes(j)>0){
                    clusterLogProb.push_back(log( ut.sampleMultivariateNormal( currState.clusterParams ) ));
             }else{
                    clusterLogProb.push_back( -INFINITY);
             }
        }


    }



}
