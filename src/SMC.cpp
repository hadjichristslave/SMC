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
    for(unsigned int t=0;t<cloudData.size();t++)
        for(int j=0;j< numOfParticles; j++)
            particles->at(j).assignments.resize(cloudData[t].size() , -1);

    for(unsigned int t=0;t<cloudData.size();t++){
        //for all particles
        for(int j=0;j< numOfParticles; j++){
            for(int s=0;s< numOfSamples;s++)
                smc_sample( & particles->at(j), cloudData[t], params, t,  s , cloudData[t].size());
            cout << "---------------------" << endl;
            for(int k=0; k<particles->at(j).clusterSizes.size();k++)
                cout << particles->at(j).clusterSizes[k].back() << ",";
            cout<< endl;
            removeEmptyStates( & particles->at(j), t);
            for(int k=0; k<particles->at(j).clusterSizes.size();k++)
                cout << particles->at(j).clusterSizes[k].back() << ",";
            cout<< endl;
        }
        resample( particles, cloudData[t],  params, t , numOfParticles);
    }


}
const void SMC::removeEmptyStates(SMC::StateProgression * state, int currTime){
    int currentObs = state->assignments.size();
    vector< vector<int> > data = state->clusterSizes;
    vector< int > sums(data.size(),0);
    for(unsigned int i = 0 ;i< data.size(); i++)
        for_each(data[i].begin(), data[i].end(), [&] (double y) mutable { sums[i] +=y; });
    int counter = 0;
    for(int i =0 ;i<sums.size();i++){
        if(sums[i]>0){
            for(int j =0;j<state->assignments.size();j++)
                state->assignments[j]        =  (state->assignments[j]==i)?counter:state->assignments[j];
            // change the index of a cluster to the reduced one.
            for(int jj =0 ; jj< state->stateProg.size(); jj++){
                if(state->stateProg[jj].size()==0) continue;
                state->stateProg[jj][counter]    = state->stateProg[jj][i];
            }
            state->clusterSizes[counter]           = state->clusterSizes[i];
            counter++;
        }
    }
    // Fix remove the redundant ones. Only removing the states of curr time due dynamic number of clusters on every step
    // CLuster that survive can be traced back though
    state->stateProg[currTime].erase (state->stateProg[currTime].begin()+counter,state->stateProg[currTime].end());
    state->clusterSizes.erase (state->clusterSizes.begin()+counter,state->clusterSizes.end());

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
    int currentClusters = currState->stateProg[currentTime>0?(currentTime-timeOffset):0].size();
    // Update cluster parameterse given the points within them
    for(int i=0;i<currentClusters;i++){
        if( currentTime>0 && currentSample==0 ){
        // Time - 2 due to the error in my ut file reader. Will make it -1 as soon as i fix the eerror
            currState->stateProg[currentTime] = currState->stateProg[currentTime-timeOffset];
        }else{
            Params par(CRP);
            if( currState->clusterSizes.back()[i] > 0 ){
                if( currentTime == 0){
                    Eigen::MatrixXd clusteredData(getDataOfCluster(i, & currState->assignments, &cloudInstance));
                    // Update the parameters given the data that belong to the cluster
                    par = updateParams(clusteredData, params, 3);

                }else{
                    par = calculatePosteriorParams(currentTime, currState,  params, &cloudInstance,i);
                    //Eigen::MatrixXd clusteredData(getDataOfCluster(i, & currState->assignments, &cloudInstance));
                    //par = updateParams(clusteredData, params, 3);
                    SMC::SufficientStatistics pr;
                    pr.covar      = ut.iwishrnd(par.tau0, par.nu0, 3, par.nu0);
                    Vector3d tempMean(ut.sampleMultivariateNormal(par.mu0, par.tau0.array()/par.kappa0,1,3));
                    pr.mean[0]     = tempMean(0);
                    pr.mean[1]     = tempMean(1);
                    pr.mean[2]     = tempMean(2);
                    pr.exponential = ut.gammarnd(params.gamma_alpha0  , params.gamma_beta0);
                    RowVectorXd temp(ut.dirrnd(params.q0));
                    pr.categorical.resize(temp.cols());
                    for(int ii=0;ii<temp.cols(); ii++) pr.categorical[ii] = temp(ii);
                    currState->stateProg[currentTime][i] = pr;
                }
            }
        }
    }
    //sample an assignment for every datapoint in the dataset
    for(int i=0;i<dataSize;i++){

        vector<double> pointInstance = cloudInstance[i];
        vector<double> clusterLogProb(currentClusters);
        int old_k = -1;
        //second time within the algorith
        if( currState->assignments.size()>(unsigned int)i)
            old_k = currState->assignments[i];
        std::vector<double> sizes;
        if(currState->clusterSizes.size() >0){
            for( unsigned int kk =0; kk< currState->clusterSizes.size(); kk++)
                sizes.push_back(currState->clusterSizes[kk].back());
        }
        double sum = 0;
        sizes.push_back(params.crp);
        for_each(sizes.begin(), sizes.end(), [&sum] (double y) mutable { sum +=y; });
        for_each(sizes.begin(), sizes.end(), [&sum] (double &y) mutable { y = y/sum; });

        vector<double> prob_assig = sizes;
        RowVector3d instance(3);
        instance << pointInstance[0], pointInstance[1], pointInstance[2];
        for( int j =0 ; j<currentClusters; j++)
             if(sizes[j]>0)
                    clusterLogProb[j] = log(ut.multivariateNormalPDF(instance, params.mu0 , params.tau0, 3));
             else   clusterLogProb[j] = -INFINITY;

        clusterLogProb.push_back(log(ut.multivariateNormalPDF(instance, params.mu0 , params.tau0, 3)));
        //TODO Must add exponent here
        sum = 0;
        for(unsigned int ik=0;ik<clusterLogProb.size();ik++) {
            prob_assig[ik] = prob_assig[ik]* exp(clusterLogProb[ik]) * ut.exppdf(pointInstance[5] , params.exp_lambda0);
            sum +=prob_assig[ik];
        }for_each(prob_assig.begin(), prob_assig.end(), [&sum] (double &y) mutable { y = y/sum;});

        int sample_k = -1;
        if(sum>0)  sample_k = ut.randcat( & prob_assig);
        else       sample_k = 0;
        currState->assignments[i] = sample_k;
        // K+ 1 because k is a cpp 0-index variable
        if( sample_k  == currentClusters){
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
// Get all the data that are assigned to cluster CLUSTER
Eigen::MatrixXd SMC::getDataOfCluster(int cluster, vector<int> * assignments , vector< vector<double> > * cloudInstance){
    int numOfAssign= 0 ;
    for(unsigned int jkl = 0 ;jkl< assignments->size();jkl++)
        if( assignments->at(jkl) == cluster)  numOfAssign++;
    Eigen::MatrixXd clusteredData(numOfAssign,cloudInstance->at(0).size());
    numOfAssign = 0;
    for(unsigned int jkl = 0 ;jkl< assignments->size();jkl++)
        if( assignments->at(jkl) == cluster){
            std::vector<double> v = cloudInstance->at(jkl);
            double* ptr = &v[0];
            Eigen::Map<Eigen::RowVectorXd> my_vect(ptr, v.size());
            clusteredData.row(numOfAssign) = my_vect;
            numOfAssign++;
        }
    return clusteredData;
}

const void SMC::newCluster(SMC::StateProgression * currState, \
                      SMC::Params params,\
                      vector<double> pointInstance,\
                      int currentTime){
    SMC::SufficientStatistics pr(pointInstance, 6);
    pr.mean[0]  = pointInstance[0];
    pr.mean[1]  = pointInstance[1];
    pr.mean[2]  = pointInstance[2];
    pr.covar = ut.iwishrnd(params.tau0, params.nu0, 3,0);
    pr.exponential = ut.gammarnd(params.gamma_alpha0  , params.gamma_beta0);
    currState->stateProg[currentTime].push_back(pr);
}

SMC::Params SMC::updateParams(MatrixXd data, SMC::Params params , int colourbins){
    SMC::Params newParams;
    if( data.rows()>0){
        MatrixXd position = data.leftCols(3);
        MatrixXd colours  = data.middleCols(6,colourbins*colourbins*colourbins);
        MatrixXd angles   = data.block(0, 3, data.rows(), 3);
        int obsSize  = position.rows();
        MatrixXd mean   = position.colwise().mean();
        MatrixXd S = ((position.row(0) - mean).transpose())* (position.row(0) - mean);
        for(int i=1;i< position.rows();i++)
            S =  S.array() + (( (position.row(i) - mean).transpose())*((position.row(i) - mean))).array();
        newParams.kappa0   = (double)params.kappa0 + obsSize;
        newParams.nu0      = (double)params.nu0    + obsSize;
        newParams.mu0 = (mean.array()* obsSize / newParams.kappa0) + (params.mu0.array() * params.kappa0/newParams.kappa0);
        newParams.tau0 = params.tau0.array() + S.array(); + ( params.kappa0 * obsSize * ( params.mu0 - mean )*(params.mu0 - mean).transpose() );
        newParams.q0 = params.q0.array()  + colours.colwise().sum().array();
        newParams.gamma_alpha0 = params.gamma_alpha0 + obsSize , newParams.gamma_beta0 =params.gamma_beta0 + angles.sum();
    }
    return newParams;
}
const void SMC::resample( vector< SMC::StateProgression > * particles, \
                          vector < vector<double> >  cloudData, \
                          SMC::Params params, \
                          int currTime,\
                          int numOfParticles){
    vector<double> weights(particles->size());

    for(int i =0 ; i< numOfParticles; i ++ ){
        weights[i] = computeWeights(& particles->at(i) , currTime , &  cloudData , params);
    }
}
double SMC::computeWeights( SMC::StateProgression * stuff, \
                            int currTime , \
                            vector<vector <double > > * cloudData , \
                            SMC::Params params){
    int clusters = stuff->stateProg[currTime].size();
    return exp(getWeightNumerator( stuff , currTime, cloudData, params)- getWeightDenominator( stuff , currTime, cloudData, params ));
}

double SMC::getWeightNumerator(SMC::StateProgression * stuff , int currTime, vector< vector<double> > * cloudData,  SMC::Params params){
    double term1 = getJointProbTheta(stuff,  currTime, cloudData, params);
    double term2 = getJointProbAssig(stuff, currTime, cloudData, params);
    double term3 = getJointProbData(stuff,  currTime, cloudData, params);
    return term1 + term2 + term3;
}
double SMC::getWeightDenominator(SMC::StateProgression * stuff , int currTime, vector< vector<double> > * cloudData, SMC::Params params){
    double term1 = getPosteriorTheta(stuff , currTime, cloudData, params);
    double term2 = getPosteriorAssignments(stuff , currTime, cloudData);
    return term1  + term2;
}

double SMC::getJointProbData(SMC::StateProgression * currState , int currTime, vector< vector<double> > * cloudData,  SMC::Params params){
    int dataSize = cloudData->size();
    RowVectorXd probAssig_i(dataSize);
    for(int i=0;i<dataSize;i++){
        Vector3d instance(3);
        instance(0) = cloudData->at(i)[0];
        instance(1) = cloudData->at(i)[1];
        instance(2) = cloudData->at(i)[2];
        int assignment = currState->assignments[i];

        Vector3d mean(3);
        mean(0) = currState->stateProg[currTime][assignment].mean[0];
        mean(1) = currState->stateProg[currTime][assignment].mean[1];
        mean(2) = currState->stateProg[currTime][assignment].mean[2];
        probAssig_i(i) = log(ut.multivariateNormalPDF( instance,\
                                                       mean ,\
                                                       currState->stateProg[currTime][assignment].covar,\
                                                       3));
    }
    return probAssig_i.sum();
}
double SMC::getJointProbAssig(SMC::StateProgression * currState , int currTime, vector< vector<double> > * cloudData,  SMC::Params params){
    int dataSize = cloudData->size();
    RowVectorXd probAssig_i(dataSize);
    for(int i=0;i<dataSize;i++){
        vector<double> pointInstance = cloudData->at(i);
        std::vector<double> sizes(0);
        if(currState->stateProg[currTime].size() >0)
            for( int kk =0; kk< currState->clusterSizes.size(); kk++)
                sizes.push_back(currState->clusterSizes[kk].back());
        double sum = 0;
        sizes.push_back(CRP);
        for_each(sizes.begin(), sizes.end(), [&sum] (double y) { sum +=y; });
        for_each(sizes.begin(), sizes.end(), [&sum] (double &y) mutable { y = y/sum; });
        vector<double> prob_assig = sizes;
        probAssig_i(i) = log(ut.catpdf(currState->assignments[i] , prob_assig));
    }
    return probAssig_i.sum();
}
double SMC::getJointProbTheta(SMC::StateProgression * currState,\
                              int currTime,\
                              vector< vector<double> > * cloudData,\
                              SMC::Params params){

 int currentClusters = currState->stateProg[currTime>0?(currTime-1):0].size();
    VectorXd postProbTheta(currentClusters);
    // Update cluster parameterse given the points within them
    for(int i=0;i<currentClusters;i++){
        Params par(CRP);
        if( currState->clusterSizes.back()[i] > 0 ){
            if( currTime == 0){
                Eigen::MatrixXd clusteredData(getDataOfCluster(i, & currState->assignments, cloudData));
                par = updateParams(clusteredData, params, 3);
            }else
                par = calculateJointParams(currTime, currState,  params, cloudData,i);
                vector<double> pr = currState->stateProg[currTime][i].mean;
                Vector3d tempMean(3);
                tempMean(0) = pr[0];
                tempMean(1) = pr[1];
                tempMean(2) = pr[2];
                postProbTheta(i) = log(ut.multivariateNormalPDF( tempMean,\
                                                                par.mu0,\
                                                                currState->stateProg[currTime][i].covar.array()/par.kappa0,\
                                                                3))+\
                                       log(ut.exppdf(currState->stateProg[currTime][i].exponential,\
                                                 par.exp_lambda0));
        }
    }
    return postProbTheta.sum();
}

double SMC::getPosteriorTheta(SMC::StateProgression * currState,\
                              int currTime,\
                              vector< vector<double> > * cloudData,\
                              SMC::Params params){
    int currentClusters = currState->stateProg[currTime>0?(currTime-1):0].size();
    VectorXd postProbTheta(currentClusters);
    // Update cluster parameterse given the points within them
    for(int i=0;i<currentClusters;i++){
        Params par(CRP);
        if( currState->clusterSizes.back()[i] > 0 ){
            if( currTime == 0){
                Eigen::MatrixXd clusteredData(getDataOfCluster(i, & currState->assignments, cloudData));
                // Update the parameters given the data that belong to the cluster
                par = updateParams(clusteredData, params, 3);
            }else
                par = calculatePosteriorParams(currTime, currState,  params, cloudData,i);
                vector<double> pr = currState->stateProg[currTime][i].mean;
                Vector3d tempMean(3);
                tempMean(0) = pr[0];
                tempMean(1) = pr[1];
                tempMean(2) = pr[2];
                postProbTheta(i) = log(ut.multivariateNormalPDF( tempMean,\
                                                                par.mu0,\
                                                                currState->stateProg[currTime][i].covar.array()/par.kappa0,\
                                                                3))+\
                                       log(ut.exppdf(currState->stateProg[currTime][i].exponential,\
                                                 par.exp_lambda0));
        }
    }

    return postProbTheta.sum();
};

double SMC::getPosteriorAssignments(SMC::StateProgression * currState , int currTime, vector< vector<double> > * cloudData){
    //sample an assignment for every datapoint in the dataset
    int dataSize = cloudData->size();
    RowVectorXd probAssig_i(dataSize);
    for(int i=0;i<dataSize;i++){
        vector<double> pointInstance = cloudData->at(i);
        vector<double> clusterLogProb(currState->stateProg[currTime].size());
        int old_k = -1;
        if( currState->assignments.size()>(unsigned int)i) old_k = currState->assignments[i];

        std::vector<double> sizes(0);
        if(currState->stateProg[currTime].size() >0)
            for( int kk =0; kk< currState->clusterSizes.size(); kk++)
                sizes.push_back(currState->clusterSizes[kk].back());
        double sum = 0;
        sizes.push_back(CRP);
        for_each(sizes.begin(), sizes.end(), [&sum] (double y) { sum +=y; });
        for_each(sizes.begin(), sizes.end(), [&sum] (double &y) mutable { y = y/sum; });
        vector<double> prob_assig = sizes;
        RowVector3d instance(3);
        instance << pointInstance[0], pointInstance[1], pointInstance[2];
        for( int j =0 ; j<currState->stateProg[currTime].size(); j++){
            RowVector3d mean(3);
            instance << currState->stateProg[currTime][j].mean[0], \
                        currState->stateProg[currTime][j].mean[1], \
                        currState->stateProg[currTime][j].mean[2];
             if(sizes[j]>0)
                   clusterLogProb[j] = log(ut.multivariateNormalPDF(instance, \
                                                                     mean, \
                                                                     currState->stateProg[currTime][j].covar, 3));
             else   clusterLogProb[j] = -INFINITY;
        }
        Eigen::RowVectorXd mean    = Eigen::RowVectorXd::Zero(3);
        Eigen::MatrixXd covar      = Eigen::MatrixXd::Identity(3,3);
        clusterLogProb.push_back(log(ut.multivariateNormalPDF(instance, mean, covar*1.0000e+20f, 3)));
        sum = 0;
        for(unsigned int ik=0;ik<currState->stateProg[currTime].size();ik++) {
            prob_assig[ik] = prob_assig[ik]* exp(clusterLogProb[ik]) \
             * ut.exppdf(pointInstance[5], currState->stateProg[currTime][ik].exponential);
        }
        probAssig_i(i) = log(ut.catpdf(currState->assignments[i] , prob_assig));
    }
    return probAssig_i.sum();
}

SMC::Params SMC::calculatePosteriorParams( int currTime,\
                                      SMC::StateProgression * currState,\
                                      SMC::Params params,\
                                      vector< vector<double> > * cloudData,
                                      int curDataPoint){

    MatrixXd data_t= getDataOfCluster(curDataPoint, & currState->assignments, cloudData);
    Eigen::MatrixXd clusteredData(data_t.rows()+ params.auxiliaryNum, data_t.cols() );
    if(clusteredData.rows()>0 &&  currState->stateProg[currTime-timeOffset].size()>0){
        SufficientStatistics ss = currState->stateProg[currTime-timeOffset][curDataPoint];
        VectorXd mean(ss.mean.size());
        mean(0) = ss.mean[0];
        mean(1) = ss.mean[1];
        mean(2) = ss.mean[2];
        //Format: params = {crp, del, #aux, tau0, v0, mu0, k0, q0, _,_,_<-#colorbins?}
        Eigen::MatrixXd auxGausSamples =  ut.sampleMultivariateNormal(mean,ss.covar, params.auxiliaryNum, 3);
        Eigen::MatrixXd auxMultSamples = ut.sampleMultinomial( ss.categorical, params.auxiliaryNum);
        Eigen::VectorXd auxExpSamples  = ut.exprnd(ss.exponential , params.auxiliaryNum);

        MatrixXd C(auxGausSamples.transpose().rows(), auxExpSamples.cols()*3+auxGausSamples.transpose().cols() +\
         auxMultSamples.cols());
        C << auxGausSamples.transpose() , auxExpSamples,auxExpSamples,auxExpSamples, auxMultSamples;
        clusteredData << data_t,  C;
        return updateParams(clusteredData, params, 3);
    }else if(currState->stateProg[currTime-timeOffset].size()>0){
        SufficientStatistics ss = currState->stateProg[currTime-timeOffset][curDataPoint];
        VectorXd mean(ss.mean.size());
        mean(0) = ss.mean[0];
        mean(1) = ss.mean[1];
        mean(2) = ss.mean[2];
        Eigen::MatrixXd auxGausSamples =  ut.sampleMultivariateNormal(mean,ss.covar, params.auxiliaryNum, 3);
        Eigen::MatrixXd auxMultSamples = ut.sampleMultinomial( ss.categorical, params.auxiliaryNum);
        Eigen::VectorXd auxExpSamples  = ut.exprnd(ss.exponential , params.auxiliaryNum);

        MatrixXd C( auxGausSamples.transpose().rows(), auxExpSamples.cols()*3 +\
                    auxGausSamples.transpose().cols() + auxMultSamples.cols());
        C << auxGausSamples.transpose() , auxExpSamples,auxExpSamples,auxExpSamples, auxMultSamples;
        return updateParams(C, params, 3);
    }else
        return  updateParams(clusteredData, params, 3);

}
SMC::Params SMC::calculateJointParams( int currTime,\
                                      SMC::StateProgression * currState,\
                                      SMC::Params params,\
                                      vector< vector<double> > * cloudData,
                                      int curDataPoint){

    MatrixXd data_t= getDataOfCluster(curDataPoint, & currState->assignments, cloudData);
    Eigen::MatrixXd clusteredData(data_t.rows()+ params.auxiliaryNum, data_t.cols() );
    if(currState->stateProg[currTime-timeOffset].size()>0){
        SufficientStatistics ss = currState->stateProg[currTime-timeOffset][curDataPoint];
        VectorXd mean(ss.mean.size());
        mean(0) = ss.mean[0];
        mean(1) = ss.mean[1];
        mean(2) = ss.mean[2];
        Eigen::MatrixXd auxGausSamples =  ut.sampleMultivariateNormal(mean,ss.covar, params.auxiliaryNum, 3);
        Eigen::MatrixXd auxMultSamples = ut.sampleMultinomial( ss.categorical, params.auxiliaryNum);
        Eigen::VectorXd auxExpSamples  = ut.exprnd(ss.exponential , params.auxiliaryNum);

        MatrixXd C( auxGausSamples.transpose().rows(), auxExpSamples.cols()*3 +\
                    auxGausSamples.transpose().cols() + auxMultSamples.cols());
        C << auxGausSamples.transpose() , auxExpSamples,auxExpSamples,auxExpSamples, auxMultSamples;
        return updateParams(C, params, 3);
    }else
        return  updateParams(clusteredData, params, 3);
}
