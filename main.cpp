#include <iostream>
#include <string>
#include <stdlib.h>
#include <stdio.h>
#include <libconfig.h>
#include <algorithm>
#include <vector>

#include "SMC.h"
#include "Utilities.h"

using namespace std;

int main(int argc, char* argv[])
{
    // Config file parsing
    int numOfParticles , numOfSamples;
    config_t cfg, *cf;
    const char *base = NULL;
    cf = &cfg;
    config_init(cf);
    if (!config_read_file(cf, "config.cfg"))
        return(EXIT_FAILURE);
    config_lookup_int(cf, "particles" , &numOfParticles);
    config_lookup_int(cf, "samples" , &numOfSamples);
    // Variable declaration
    SMC smc;
    Utilities ut;
    cout << "Stuff;" << endl;
    vector< vector< vector< double > > > dataPoints  = ut.readFile("-----");
    int timeStates = dataPoints.size(); // All different points in time of our pointclouds.

    //Format: params = {crp, del, #aux, tau0, v0, mu0, k0, q0, _,_,_<-#colorbins?,lambda0(angle distance measure)}
    //State format: = { assignments, cluster parameters, clusterSizes(Number of elements) }
    //RowVector3d x;
    //x << -.886363 ,-.808702, -.0398121;
    //RowVector3d mu;
    //mu <<  1.20956e-316, 1.4822e-323, 1.4822e-323;
    //Matrix3d covar;
    //covar  << .00560623, -.00302697, .000888741,
             //-.00302697, .0017546, -.000225909,
              //.000888741, -.000225909,  .00351697;
    //cout << ut.iwishrnd( covar,60 ,3,1)<< endl;
    //return 0;
    SMC::Params Baseparams;
    Baseparams.cloudInstances = dataPoints.size();
    vector < SMC::StateProgression > particles(numOfParticles, timeStates);
    smc.infer( &particles, dataPoints, Baseparams, numOfParticles , numOfSamples);
    return 0;
}
