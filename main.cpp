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

    if (!config_read_file(cf, "config.cfg")) {
        fprintf(stderr, "%s:%d - %s\n",
        config_error_file(cf),
        config_error_line(cf),
        config_error_text(cf));
        config_destroy(cf);
        return(EXIT_FAILURE);
    }
    config_lookup_int(cf, "particles" , &numOfParticles);
    config_lookup_int(cf, "samples" , &numOfSamples);

    // Variable declaration
    SMC smc;
    Utilities ut;
    cout << "parsing file..." <<  endl;
    vector< vector< vector< double > > > dataPoints  = ut.readFile("-----");
    cout << "parsing of the file succesfull!"  << endl;
    cout << "Procceding with the clustering.." << endl;

    //Format: params = {crp, del, #aux, tau0, v0, mu0, k0, q0, _,_,_<-#colorbins?,lambda0(angle distance measure)}
    SMC::Params params;

    params.cloudInstances = dataPoints.size();
    //State format: = { assignments, cluster parameters, clusterSizes(Number of elements) }
    SMC::State state(dataPoints.size());

    vector< SMC::State > particles(numOfParticles);

    smc.infer(particles, dataPoints, params, numOfParticles , numOfSamples);



//    vector<double> vec(10,100);
//    for_each(vec.begin(), vec.end(), [] (int y) mutable {
//        y++;
//        cout << y << endl;
//    });




    return 0;
}
