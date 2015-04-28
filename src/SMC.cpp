#include "SMC.h"

SMC::SMC()
{
    //ctor

}

SMC::~SMC()
{
    //dtor
}

const void SMC::infer(vector< StateProgression >  * particles, \
                      vector < vector < vector<double> > > cloudData, \
                      SMC::Params params, int numOfParticles,\
                      int numOfSamples){
//for all time intervals
    for(int t=0;t<cloudData.size();t++)
        for(int j=0;j< numOfParticles; j++)
            particles->at(j).assignments.resize(cloudData[t].size() , -1);

    for(int t=0;t<cloudData.size();t++){
        //for all particles
        cout << "T ->" << t<< endl;
        if(t%2==1) continue;
        for(int j=0;j< numOfParticles; j++)
            for(int s=0;s< numOfSamples;s++)
                smc_sample( & particles->at(j), cloudData[t], params, t,  s , cloudData[t].size());

    }


}
// Sample  cluster parameters for cloud data at a given time
const void SMC::smc_sample(StateProgression  * currState, \
                vector< vector<double> >cloudInstance, \
                SMC::Params params, \
                int currentTime, \
                int currentSample , \
                int dataSize){

    //Get the current clusters.
    // Time - 2 due to the error in my ut file reader. Will make it -1 as soon as i fix the eerror
    int currentClusters = currState->stateProg[currentTime>0?(currentTime-2):0].size();
    int tempTime        = currentTime>0?(currentTime-1):0;
    // Update cluster parameterse given the points within them
    for(int i=0;i<currentClusters;i++){
        if( currentTime>0 && currentSample==0 ){
        // Time - 2 due to the error in my ut file reader. Will make it -1 as soon as i fix the eerror
            cout <<" got in if" << endl;
            currState->stateProg[currentTime] = currState->stateProg[currentTime-2];
        }else{
            cout <<" got in else" << endl;
            if( currState->clusterSizes.back()[i] > 0 ){
                if( currentTime == 0){
                    vector< vector< double> > dataOfCluster;
                    for( int jkl = 0 ;jkl< currState->assignments.size();jkl++)
                        if( assignments[jkl] == i)
                            dataOfCluster.push_back( cloudInstance[jkl]);
                    updateParams()

                }else{



                }




            }
        }

    }


    //sample an assignment for every datapoint in the dataset
    for(int i=0;i<dataSize;i++){
        vector<double> pointInstance = cloudInstance[i];
        vector<double> clusterLogProb(currentClusters);
        int old_k = -1;
        //second time within the algorithm
        if( currState->assignments.size()>=i+1)
            old_k = currState->assignments[i];

        std::vector<double> sizes;
        if(currState->stateProg[currentTime].size() >0)
            for( int kk =0; kk< currState->clusterSizes.size(); kk++)
                sizes.push_back(currState->clusterSizes[kk].back());

        double sum = 0;
        sizes.push_back(params.crp);
        for_each(sizes.begin(), sizes.end(), [&sum] (double y) mutable { sum +=y; });
        for_each(sizes.begin(), sizes.end(), [&sum] (double &y) mutable { y = y/sum; });
        vector<double> prob_assig = sizes;
        RowVector3d instance(3);
        instance << pointInstance[0], pointInstance[1], pointInstance[2];

        for( int j =0 ; j<currentClusters; j++){
             if(sizes[j]>0)
                    clusterLogProb[j] = log(ut.sampleMultivariateNormal(instance, params.mu , params.tau0, 3));
             else
                    clusterLogProb[j] = -INFINITY;

        }
        RowVector3d       BaseMu = RowVector3d::Zero(1,3);
        Matrix3d          BaseTau = MatrixXd::Identity(3,3);
        BaseTau           *= 1.0000e+20f;
        clusterLogProb.push_back(log(ut.sampleMultivariateNormal(instance, BaseMu , BaseTau, 3)));
        sum =0;
        for(int ik=0;ik<clusterLogProb.size();ik++) {
            prob_assig[ik] = prob_assig[ik]* exp(clusterLogProb[ik]);
            sum +=prob_assig[ik];
        }
        for_each(prob_assig.begin(), prob_assig.end(), [&sum] (double &y) mutable { y = y/sum;});
        int sample_k = -1;
        if(sum>0)
            sample_k = ut.randcat( & prob_assig);
        else
            sample_k =0;

        currState->assignments[i] = sample_k;
        if( sample_k +1  >= currentClusters+1){
            currentClusters++;
            newCluster(currState, params, pointInstance, currentTime);
            vector<int> sizeVec(cloudInstance.size(),0);
            currState->clusterSizes.push_back(sizeVec);
        }
        if(old_k!=-1)
            for(int k=i; k<dataSize;k++)
                currState->clusterSizes[old_k][k]--;
        for(int k=i; k<dataSize;k++)
            currState->clusterSizes[sample_k][k]++;
    }

}

const void SMC::newCluster(SMC::StateProgression * currState, \
                      SMC::Params params,\
                      vector<double> pointInstance,\
                      int currentTime){
    SMC::SufficientStatistics pr(pointInstance, 6);
    pr.mean[0]  = pointInstance[0];
    pr.mean[1]  = pointInstance[1];
    pr.mean[2]  = pointInstance[2];
    pr.covar = ut.iwishrnd(params.tau0, params.nu0, params.tau0.rows());
    pr.exponential = pointInstance[3];
    currState->stateProg[currentTime].push_back(pr);
}

const void SMC::updateParams(MatrixXd position, MatrixXd colourInformation, MatrixXd angleInformation, SMC::Params params){


}
