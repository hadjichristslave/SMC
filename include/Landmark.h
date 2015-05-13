#ifndef LANDMARK_H
#define LANDMARK_H
#include "SMC.h"

class Landmark
{
    int uuid;
    SMC::SufficientStatistics distribution;
    public:
        Landmark(int LandId , SMC::SufficientStatistics stats);
        virtual ~Landmark();
        int getId();
        void setId(int id );
        void print();
    protected:
    private:
};

#endif // LANDMARK_H
