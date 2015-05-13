#include "Landmark.h"

Landmark::Landmark(int LandId , SMC::SufficientStatistics stats)
{
    //ctor

    int uuid = LandId;
    SMC::SufficientStatistics distribution = stats;
}

Landmark::~Landmark()
{
    //dtor
}
int Landmark::getId(){
    return uuid;
}
void Landmark::setId(int id ){
    uuid = id;
}

void Landmark::print(){
    cout << "Landmark id " << endl;
    cout << " With distribution " << endl << distribution << endl;
}
