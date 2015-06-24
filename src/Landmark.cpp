#include "Landmark.h"

Landmark::Landmark(){
    uuid = 0;
}

Landmark::Landmark(int LandId , SufficientStatistics stats, int initialIdz){
    uuid = LandId;
    distribution.copy(stats);
    initialId = initialIdz;
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
    cout << "Landmark id " << uuid << endl;
    cout << " With distribution " << endl << distribution << endl;
}
