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
    for(int t=0;t<cloudData.size();t++)
        for(int j=0;j< numOfParticles; j++)
            particles[j].assignments.resize(cloudData[t].size());



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
    int currentClusters = currState.clusterParams.size();

    // Update cluster parameterse given the points within them

    for(int i=0;i<currentClusters;i++){

    }


    //sample an assignment for every datapoint in the dataset
    for(int i=0;i<cloudInstance.size();i++){
        vector<double> pointInstance = cloudInstance[i];
        vector<double> clusterLogProb(currentClusters);


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
        vector<double> prob_assig = sizes;

        RowVector3d instance(3);
        instance << pointInstance[0], pointInstance[1], pointInstance[2];


        for( int j =0 ; j<currentClusters; j++){
             if(sizes[j]>0){
                    clusterLogProb[j] = log(ut.sampleMultivariateNormal(instance, params.mu , params.tau0, 3));
             }else{
                    clusterLogProb[j] = -INFINITY;
             }
        }
        RowVector3d       BaseMu = RowVector3d::Zero(1,3);
        Matrix3d          BaseTau = MatrixXd::Identity(3,3);
        BaseTau           *= 1.0000e+20f;
        clusterLogProb.push_back(log(ut.sampleMultivariateNormal(instance, BaseMu , BaseTau, 3)));

        sum =0;
        for(int ik=0;ik< clusterLogProb.size();ik++) {
            prob_assig[ik] = prob_assig[ik]* exp(clusterLogProb[ik]);
            sum +=prob_assig[ik];
        }
        for_each(prob_assig.begin(), prob_assig.end(), [&sum] (double &y) mutable { y = y/sum;});
        int sample_k = -1;
        if(sum>0)
            sample_k = ut.randcat( & prob_assig);
        else
            sample_k =0;

        cout << currState.assignments.size() << " - " << i << endl;
        currState.assignments[i] = sample_k;

        if( sample_k == currentClusters+1 ){
            currentClusters++;
            SMC::State tempState = newCluster(currState, params, pointInstance, currentTime);

        }


    }

}

SMC::State SMC::newCluster(SMC::State currState, \
                      SMC::Params params,\
                      vector<double> pointInstance,\
                      int currentTime){
    int K = currState.clusterParams.size() + 1;

    const gsl_rng_type * T;
    gsl_rng * r;
    gsl_rng_env_setup();

    T = gsl_rng_default;
    r = gsl_rng_alloc(T);




    SMC::State stet(10);
    return stet;
}
