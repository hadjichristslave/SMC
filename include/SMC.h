#ifndef SMC_H
#define SMC_H
#include <eigen3/Eigen/Dense>
#include "Utilities.h"
#include "Point.h"
#include "PointCloud.h"

using namespace std;
using namespace Eigen;

class SMC
{
    public:
        SMC();
        virtual ~SMC();

        struct Params{
            int   cloudInstances, auxiliaryNum , colourBin , cBin1, cBin2, cBin3, colourBins = 3;
            double crp, del, nu0, kappa0, KullbackDistance, EMDDistance;
            Vector3d  mu = Vector3d(.0,.0,.0);
            MatrixXd tau0 = MatrixXd::Identity(3,3);
            RowVectorXd q0 = RowVectorXd::Ones(1,colourBins*colourBins*colourBins);

            //Format: params = {crp, del, #aux,  tau0, v0, mu0, k0, q0, _,_,_<-#colorbin
            //params = {0.1, 0.7, 10, 3*eye(2), 60, [0,0], 0.05, 10*ones(1,10), 1, 1, 1};
            Params(void)
            {
                crp = .1;
                del = .7;
                auxiliaryNum = 10;
                nu0 = 60;
                kappa0 - .05;
                cBin1  = 1,cBin2 =1, cBin3  = 1;
              // Initialize Foo
            }
        };
        struct State{
            vector<int>    assignments;
            vector<Params> clusterParams;
            vector<double>    clusterSizes; // Double despite them being integers because I'm bored to do typecasting on every iteration :)

            State(int numsamples){
                assignments.resize(numsamples,0);
            }
            State(void){

            }
        };

        Utilities ut;

        const void infer( vector< SMC::State > particles, \
                          vector < vector < vector<double> > > cloudData, \
                          SMC::Params params, int numOfParticles,\
                          int numOfSamples);
        const void smc_sample(SMC::State currState, \
                                vector< vector<double> >cloudInstance, \
                                SMC::Params params, \
                                int currentTime,
                                int currentSample, \
                                int dataSize);


    protected:
    private:
};

#endif // SMC_H
