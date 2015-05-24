#ifndef LANDMARK_H
#define LANDMARK_H
#include "Structures.h"

using namespace Structures;

class Landmark
{
    public:
        int uuid;
        SufficientStatistics distribution;
        Landmark(int LandId , SufficientStatistics stats);
        Landmark();

        virtual ~Landmark();
        int getId();
        void setId(int id );
        void print();
    protected:
    private:
};

#endif // LANDMARK_H
