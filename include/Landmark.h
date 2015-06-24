#ifndef LANDMARK_H
#define LANDMARK_H
#include "Structures.h"

using namespace Structures;

class Landmark
{
    public:
        int uuid,initialId;
        SufficientStatistics distribution;
        Landmark(int LandId , SufficientStatistics stats,int initialId);
        Landmark();

        virtual ~Landmark();
        int getId();
        void setId(int id );
        void print();
    protected:
    private:
};

#endif // LANDMARK_H
