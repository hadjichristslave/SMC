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
    int currentClusters = currState->stateProg[currentTime>0?(currentTime-timeOffset):0].size();
    // Update cluster parameterse given the points within them
    for(int i=0;i<currentClusters;i++){
        if( currentTime>0 && currentSample==0 ){
        // Time - 2 due to the error in my ut file reader. Will make it -1 as soon as i fix the eerror
            currState->stateProg[currentTime] = currState->stateProg[currentTime-timeOffset];
        }else{
            Params par(.1);
            if( currState->clusterSizes.back()[i] > 0 ){
                if( currentTime == 0){
                    Eigen::MatrixXd clusteredData(getDataOfCluster(i, & currState->assignments, &cloudInstance));
                    // Update the parameters given the data that belong to the cluster
                    par = updateParams(clusteredData, params, 3);

                }else{
                    MatrixXd data_t= getDataOfCluster(i, & currState->assignments, &cloudInstance);
                    Eigen::MatrixXd clusteredData(data_t.rows()+ params.auxiliaryNum, data_t.cols() );
                    //Eigen::MatrixXd clusteredData(getDataOfCluster(i, & currState->assignments, &cloudInstance));
                    //par = updateParams(clusteredData, params, 3);

                    if(clusteredData.rows()>0 &&  currState->stateProg[currentTime-timeOffset].size()>0){
                        SufficientStatistics ss = currState->stateProg[currentTime-timeOffset][i];
                        // This is an ugly way to initialize an Eigen vector to an std::vector.
                        // Must search alternatives
                        VectorXd mean(ss.mean.size());
                        mean(0) = ss.mean[0];
                        mean(1) = ss.mean[1];
                        mean(2) = ss.mean[2];
                        //Format: params = {crp, del, #aux, tau0, v0, mu0, k0, q0, _,_,_<-#colorbins?}
                        Eigen::MatrixXd auxGausSamples =  ut.sampleMultivariateNormal(mean,ss.covar, params.auxiliaryNum, 3);
                        Eigen::MatrixXd auxMultSamples = ut.sampleMultinomial( ss.categorical, params.auxiliaryNum);
                        Eigen::VectorXd auxExpSamples  = ut.exprnd(ss.exponential , params.auxiliaryNum);

                        MatrixXd C(auxGausSamples.transpose().rows(), auxExpSamples.cols()*3+auxGausSamples.transpose().cols() + auxMultSamples.cols());
                        C << auxGausSamples.transpose() , auxExpSamples,auxExpSamples,auxExpSamples, auxMultSamples;
                        clusteredData << data_t,  C;
                        par = updateParams(clusteredData, params, 3);

                    }else if(currState->stateProg[currentTime-timeOffset].size()>0){
                        SufficientStatistics ss = currState->stateProg[currentTime-timeOffset][i];
                        VectorXd mean(ss.mean.size());
                        mean(0) = ss.mean[0];
                        mean(1) = ss.mean[1];
                        mean(2) = ss.mean[2];
                        Eigen::MatrixXd auxGausSamples =  ut.sampleMultivariateNormal(mean,ss.covar, params.auxiliaryNum, 3);
                        Eigen::MatrixXd auxMultSamples = ut.sampleMultinomial( ss.categorical, params.auxiliaryNum);
                        Eigen::VectorXd auxExpSamples  = ut.exprnd(ss.exponential , params.auxiliaryNum);

                        MatrixXd C(auxGausSamples.transpose().rows(), auxExpSamples.cols()*3+auxGausSamples.transpose().cols() + auxMultSamples.cols());
                        C << auxGausSamples.transpose() , auxExpSamples,auxExpSamples,auxExpSamples, auxMultSamples;
                        par = updateParams(C, params, 3);

                    }else{
                        par = updateParams(clusteredData, params, 3);

                    }
//                  Format: params = {crp, del, #aux, tau0, v0, mu0, k0, q0, _,_,_<-#colorbins?}
//                  state{2}{t,k}{2} = iwishrnd(params_post{4},params_post{5});
//                  state{2}{t,k}{1} = mvnrnd(params_post{6},state{2}{t,k}{2}/params_post{7});
//                  state{2}{t,k}{3} = dirrnd(params_post{8},1);
//                  state{2}{t,k}{3} = gammarnd(alpha, beta);
                    SMC::SufficientStatistics pr;
                    pr.covar = ut.iwishrnd(par.tau0, par.nu0, 3, par.nu0);
                    pr.mean  = ut.sampleMultivariateNormal(par.mu0, par.tau0,1,3);
                    //pr.exponential = ut.gammarnd(params.gamma_alpha0  , params.gamma_beta0);
                    //pr.covar = ut.iwishrnd(params.tau0, params.nu0, params.tau0.rows());
                    //currState->stateProg[currentTime] = iwishrnd();


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

        for( int j =0 ; j<currentClusters; j++)
             if(sizes[j]>0)
                    clusterLogProb[j] = log(ut.multivariateNormalPDF(instance, params.mu0 , params.tau0, 3));
             else   clusterLogProb[j] = -INFINITY;

        clusterLogProb.push_back(log(ut.multivariateNormalPDF(instance, params.mu0 , params.tau0*1.0000e+20f, 3)));
        //TODO Must add exponent here
        sum = 0;
        for(int ik=0;ik<clusterLogProb.size();ik++) {
            prob_assig[ik] = prob_assig[ik]* exp(clusterLogProb[ik]) * ut.exppdf(pointInstance[6] , params.exp_lambda0);
            sum +=prob_assig[ik];
        } for_each(prob_assig.begin(), prob_assig.end(), [&sum] (double &y) mutable { y = y/sum;});

        int sample_k = -1;
        if(sum>0)  sample_k = ut.randcat( & prob_assig);
        else       sample_k =0;

        currState->assignments[i] = sample_k;
        // K+ 1 because k is a cpp 0-index variable
        if( sample_k +1  > currentClusters){
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
    for( int jkl = 0 ;jkl< assignments->size();jkl++)
        if( assignments->at(jkl) == cluster)  numOfAssign++;
    Eigen::MatrixXd clusteredData(numOfAssign,cloudInstance->at(0).size());
    numOfAssign = 0;
    for( int jkl = 0 ;jkl< assignments->size();jkl++)
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
    SMC::Params newParams(.1);
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
    return newParams;
}
