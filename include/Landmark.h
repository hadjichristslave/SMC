#ifndef LANDMARK_H
#define LANDMARK_H
#include "SMC.h"

class Landmark
{
    public:
        int uuid;
        SMC::SufficientStatistics distribution;
        Landmark(int LandId , SMC::SufficientStatistics stats);

        virtual ~Landmark();
        int getId();
        void setId(int id );
        void print();
    protected:
    private:
};

#endif // LANDMARK_H
