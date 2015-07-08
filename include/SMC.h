#ifndef SMC_H
#define SMC_H
#include "Utilities.h"

using namespace std;
using namespace Eigen;
using namespace Structures;
class SMC
{
    private:
        const static int timeOffset = 1;
        const static unsigned long delaytime = 1;
    public:
        SMC();
        virtual ~SMC();
        Utilities ut;

        const void infer( vector< StateProgression > * particles, \
                          vector < vector < vector<double> > > * cloudData, \
                          Params params, int numOfParticles,\
                          int numOfSamples);
        const void sample(  StateProgression  * currState, \
                                vector< vector<double> > * cloudInstance, \
                                Params params, \
                                int currentTime,
                                int currentSample, \
                                int dataSize);
        const void newCluster(StateProgression * currState,\
                                   Params params, \
                                   vector<double> pointInstance,\
                                   int currentTime);
        Params updateParams(MatrixXd data, \
                                Params params , \
                                int colourbins);
        Eigen::MatrixXd getDataOfCluster(int cluster, vector<int> * assignments,\
                                         vector< vector<double> > * cloudInstance);
        const void removeEmptyStates(StateProgression * state, int currTime);

        const void resample(         vector< StateProgression > * particles, \
                                     vector < vector<double> > * cloudData, \
                                     Params params,
                                     int currTime,
                                     int numOfParticles);
        double computeWeights(       StateProgression * stuff,\
                                     int currTime,\
                                     vector< vector<double> > * cloudData,\
                                     Params params);
        double getWeightNumerator(   StateProgression * stuff,\
                                     int currTime,\
                                     vector< vector<double> > * cloudData,\
                                     Params params);
        double getWeightDenominator( StateProgression * stuff ,\
                                     int currTime,\
                                     vector< vector<double> > * cloudData,\
                                     Params params);
        double getJointProbTheta(    StateProgression * currState,\
                                     int currTime,\
                                     vector< vector<double> > * cloudData,\
                                     Params params);
        double getJointProbAssig(    StateProgression * currState,\
                                     int currTime,\
                                     vector< vector<double> > * cloudData,\
                                     Params params);
        double getJointProbData(     StateProgression * currState,\
                                     int currTime,\
                                     vector< vector<double> > * cloudData,\
                                     Params params);
        double getPosteriorTheta(    StateProgression * currState,\
                                     int currTime,\
                                     vector< vector<double> > * cloudData,\
                                       Params params);
        double getPosteriorAssignments(StateProgression * currState,\
                                       int currTime,\
                                       vector< vector<double> > * cloudData,\
                                       Params par);

        Params calculatePosteriorParams( int currTime,\
                                      StateProgression * currState,\
                                      Params params,\
                                      vector< vector<double> > * cloudData,
                                      int curDatapoint);
        Params calculateJointParams( int currTime,\
                                      StateProgression * currState,\
                                      Params params,\
                                      vector< vector<double> > * cloudData,
                                      int curDataPoint);
        int numOfLandmarks;
    protected:

};
#endif // SMC_H
