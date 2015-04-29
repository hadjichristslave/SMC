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
    vector< vector< vector< double > > > dataPoints  = ut.readFile("-----");
    int timeStates = dataPoints.size(); // All different points in time of our pointclouds.
    //Format: params = {crp, del, #aux, tau0, v0, mu0, k0, q0, _,_,_<-#colorbins?,lambda0(angle distance measure)}
    //State format: = { assignments, cluster parameters, clusterSizes(Number of elements) }
    SMC::Params Baseparams;
    Baseparams.cloudInstances = dataPoints.size();
    vector < SMC::StateProgression > particles(timeStates, timeStates);
    smc.infer( &particles, dataPoints, Baseparams, numOfParticles , numOfSamples);
    return 0;
}
