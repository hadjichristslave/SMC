#ifndef LANDMARKS_H
#define LANDMARKS_H
#include "Landmark.h"
#include "Utilities.h"



using namespace std;
class Landmarks
{
    public:
        Landmarks();
        virtual ~Landmarks();
        vector<Landmark> landmarks;
        void addLandMark( Landmark land);
        vector< vector< double > > extractDistances(Landmark * observation , Utilities * ut);
        unsigned int size();
        void print(int landmarkIndex);
    protected:
    private:
};

#endif // LANDMARKS_H
