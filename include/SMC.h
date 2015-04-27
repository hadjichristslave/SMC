#ifndef SMC_H
#define SMC_H
#include <eigen3/Eigen/Dense>
#include "Utilities.h"
#include "Point.h"
#include "PointCloud.h"
#include <iostream>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>


using namespace std;
using namespace Eigen;

class SMC
{
    public:
        SMC();
        virtual ~SMC();


        struct SufficientStatistics{
            std::vector<double> mean        = decltype(mean)(3,0);
            Matrix3d            covar       = MatrixXd::Identity(3,3);
            vector<double>      categorical;
            double              exponential;
            SufficientStatistics(vector<double> pointInstance , int begin){
                double sum = 0;
                for(int i =begin-1; i < pointInstance.size(); i++){
                    categorical.push_back(pointInstance[i]);
                    sum += pointInstance[i];
                }
                if(sum!=0)
                    for_each(categorical.begin(), categorical.end(), [&sum] (double &y) mutable { y = y/sum;});
            }
            friend ostream& operator<<(ostream& out ,const SufficientStatistics& rhs){
                out << "========Start of SS print======"<< endl;
                out <<" Mean statistics";
                out << rhs.mean[0] << "," << rhs.mean[1] << "," << rhs.mean[2] << endl;
                out << "Covariance statistics" << endl;
                out << rhs.covar;
                out << endl;
                for( int i =0 ; i< rhs.categorical.size();i++)
                    out << rhs.categorical[i] << ",";
                out << endl;
                out << " Exponent" << rhs.exponential<< endl;
                out << "==========End of SS print====="<< endl;
                return out;
            }
        };
        struct Params{
            std::vector<double> position = decltype(position)(3,0);
            int   cloudInstances, auxiliaryNum , colourBin , cBin1, cBin2, cBin3, colourBins = 3;
            double crp, del, nu0, kappa0, KullbackDistance, EMDDistance;
            RowVector3d  mu = RowVector3d::Zero(1,3);
            Matrix3d tau0 = MatrixXd::Identity(3,3);
            RowVectorXd q0 = RowVectorXd::Ones(1,colourBins*colourBins*colourBins);

            //Format: params = {crp, del, #aux,  tau0, v0, mu0, k0, q0, _,_,_<-#colorbin
            //params = {0.1, 0.7, 10, 3*eye(2), 60, [0,0], 0.05, 10*ones(1,10), 1, 1, 1};
            Params(void)
            {
                crp = .1;
                del = .7;
                auxiliaryNum = 10;
                nu0 = 60;
                kappa0 = .05;
                cBin1  = 1,cBin2 =1, cBin3  = 1;
              // Initialize Foo
            }
        };
        struct State{
            Params clusterParams;
        };
        struct StateProgression{
            vector< vector< SMC::SufficientStatistics > > stateProg;// Outside vector defines time index
                                                      // Inside vector defines cluster index
                                                      // Params define cluster statistics at time t at cluster k
            vector<int>                    assignments;
            vector < vector<double> >      clusterSizes; // Double despite them being integers because I'm bored to do typecasting on every iteration :)
            StateProgression(int timeWindow){
                stateProg.resize(timeWindow);
            }
        };

        Utilities ut;

        const void infer( vector< SMC::StateProgression > particles, \
                          vector < vector < vector<double> > > cloudData, \
                          SMC::Params params, int numOfParticles,\
                          int numOfSamples);
        const void smc_sample(  SMC::StateProgression  * currState, \
                                vector< vector<double> >cloudInstance, \
                                SMC::Params params, \
                                int currentTime,
                                int currentSample, \
                                int dataSize);
        const void newCluster(SMC::StateProgression * currState,\
                                   SMC::Params params, \
                                   vector<double> pointInstance,\
                                   int currentTime);

    protected:
    private:
};

#endif // SMC_H
