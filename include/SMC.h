#ifndef SMC_H
#define SMC_H
#include <eigen3/Eigen/Dense>
#include "Utilities.h"
#include <iostream>
#include <fstream>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

using namespace std;
using namespace Eigen;

#if defined(__WIN32__) || defined(_WIN32) || defined(WIN32) || defined(__WINDOWS__) || defined(__TOS_WIN__)

  #include <windows.h>

  inline void delay( unsigned long ms )
    {
    Sleep( ms );
    }

#else  /* presume POSIX */

  #include <unistd.h>

  inline void delay( unsigned long ms )
    {
    usleep( ms * 1000 );
    }

#endif
class SMC
{
    private:
        const static int CRP = 1;
        const static int timeOffset = 1;
        const static unsigned long delaytime = 1;
    public:
        SMC();
        virtual ~SMC();
        const void init();
        struct SufficientStatistics{
            std::vector<double> mean        = decltype(mean)(3,0);
            Matrix3d            covar       = MatrixXd::Identity(3,3);
            vector<double>      categorical;
            double              exponential;
            SufficientStatistics(vector<double> pointInstance , int begin){
                double sum = 0;
                for(unsigned int i =begin; i < pointInstance.size(); i++){
                    categorical.push_back(pointInstance[i]);
                    sum += pointInstance[i];
                }
                if(sum!=0)
                    for_each(categorical.begin(), categorical.end(), [&sum] (double &y) mutable { y = y/sum;});
            }
            SufficientStatistics(void){}
            void copy(SufficientStatistics stat){
                mean        = stat.mean;
                covar       = stat.covar;
                categorical = stat.categorical;
                exponential = stat.exponential;
            }
            friend ostream& operator<<(ostream& out ,const SufficientStatistics& rhs){
                out << "====="<< endl;
                //out <<" Mean statistics " << endl;
                out << rhs.mean[0] << "," << rhs.mean[1] << "," << rhs.mean[2] << endl;
                //out << "-" << endl;
                //out << "Covariance statistics" << endl;
                for( int i = 0 ; i < rhs.covar.rows(); i ++){
                    for( int j = 0 ; j < rhs.covar.cols(); j ++){
                        out << rhs.covar(i,j) ;
                        out << "," ;
                    }
                }
                out << endl;
                //out << endl;
                //for(unsigned int i =0 ; i< rhs.categorical.size();i++)
                    //out << rhs.categorical[i] << ",";
                //out << endl;
                //out << " Exponent" << rhs.exponential<< endl;
                //out << "==========End of SS print====="<< endl;
                return out;
            }
        };
        struct Params{
            std::vector<double> position = decltype(position)(3,0);
            int   cloudInstances, auxiliaryNum , colourBin , cBin1, cBin2, cBin3, colourBins = 3;
            double  crp, del, nu0, kappa0, gamma_alpha0, gamma_beta0 , exp_lambda0 =1;
            RowVector3d  mu0 = RowVector3d::Zero(1,3);
            Matrix3d tau0 = MatrixXd::Identity(3,3);
            RowVectorXd q0 = RowVectorXd::Ones(1,colourBins*colourBins*colourBins);
            Params(void){
                crp = CRP ; del = .7;auxiliaryNum = 10;
                nu0 = 60; kappa0 = .05;
                cBin1  = 1,cBin2 =1, cBin3  = 1;
                // For my exponential alpha and beta are the parameters of the prior distribution Gamma
                // alpha is updated by the numver of observations whereas  beta by their sum
                gamma_alpha0 = 1;  gamma_beta0 = 1;
            }
            Params(double a) : crp(a){}
            friend ostream& operator<<(ostream& out ,const Params& rhs){
                out << "----------------Params printing-------------------" << endl;
                out << " mu " << rhs.mu0 << endl;
                out << " nu " << rhs.nu0 << endl;
                out << " tau " << rhs.tau0 << endl;
                out << " q0 " << rhs.q0 << endl;
                out << " kappa0 " << rhs.kappa0 << endl;
                out << " gamma alpha " << rhs.gamma_alpha0 << endl;
                out << " gamma beta  " << rhs.gamma_beta0 << endl;
                out << "----------------Params printing end-------------------" << endl;
                return out;
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
            vector < vector<int> >         clusterSizes;
            StateProgression(int timeWindow){
                stateProg.resize(timeWindow);
            }
        };
        Utilities ut;

        const void infer( vector< SMC::StateProgression > * particles, \
                          vector < vector < vector<double> > > cloudData, \
                          SMC::Params params, int numOfParticles,\
                          int numOfSamples);
        const void sample(  SMC::StateProgression  * currState, \
                                vector< vector<double> >cloudInstance, \
                                SMC::Params params, \
                                int currentTime,
                                int currentSample, \
                                int dataSize);
        const void newCluster(SMC::StateProgression * currState,\
                                   SMC::Params params, \
                                   vector<double> pointInstance,\
                                   int currentTime);
        SMC::Params updateParams(MatrixXd data, \
                                SMC::Params params , \
                                int colourbins);
        Eigen::MatrixXd getDataOfCluster(int cluster, vector<int> * assignments,\
                                         vector< vector<double> > * cloudInstance);
        const void removeEmptyStates(SMC::StateProgression * state, int currTime);

        const void resample(         vector< SMC::StateProgression > * particles, \
                                     vector < vector<double> > cloudData, \
                                     SMC::Params params,
                                     int currTime,
                                     int numOfParticles);
        double computeWeights(       SMC::StateProgression * stuff,\
                                     int currTime,\
                                     vector< vector<double> > * cloudData,\
                                     SMC::Params params);
        double getWeightNumerator(   SMC::StateProgression * stuff,\
                                     int currTime,\
                                     vector< vector<double> > * cloudData,\
                                     SMC::Params params);
        double getWeightDenominator( SMC::StateProgression * stuff ,\
                                     int currTime,\
                                     vector< vector<double> > * cloudData,\
                                     SMC::Params params);
        double getJointProbTheta(    SMC::StateProgression * currState,\
                                     int currTime,\
                                     vector< vector<double> > * cloudData,\
                                     SMC::Params params);
        double getJointProbAssig(    SMC::StateProgression * currState,\
                                     int currTime,\
                                     vector< vector<double> > * cloudData,\
                                     SMC::Params params);
        double getJointProbData(     SMC::StateProgression * currState,\
                                     int currTime,\
                                     vector< vector<double> > * cloudData,\
                                     SMC::Params params);
        double getPosteriorTheta(    SMC::StateProgression * currState,\
                                     int currTime,\
                                     vector< vector<double> > * cloudData,\
                                       SMC::Params params);
        double getPosteriorAssignments(SMC::StateProgression * currState,\
                                       int currTime,\
                                       vector< vector<double> > * cloudData,\
                                       SMC::Params par);

        SMC::Params calculatePosteriorParams( int currTime,\
                                      SMC::StateProgression * currState,\
                                      SMC::Params params,\
                                      vector< vector<double> > * cloudData,
                                      int curDatapoint);
        SMC::Params calculateJointParams( int currTime,\
                                      SMC::StateProgression * currState,\
                                      SMC::Params params,\
                                      vector< vector<double> > * cloudData,
                                      int curDataPoint);
        int numOfLandmarks;
    protected:

};

#endif // SMC_H
