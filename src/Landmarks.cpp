#include "Landmarks.h"

Landmarks::Landmarks()
{
    //ctor
}

Landmarks::~Landmarks()
{
    //dtor
}

void Landmarks::addLandMark(Landmark land){
    landmarks.push_back( land );
}
unsigned int  Landmarks::size(){
    return landmarks.size();
}
void Landmarks::print(int landmarkIndex){
    landmarks[landmarkIndex].print();
}
