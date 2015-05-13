#ifndef LANDMARKS_H
#define LANDMARKS_H
#include "Landmark.h"

using namespace std;
class Landmarks
{
    public:
        Landmarks();
        virtual ~Landmarks();
        vector<Landmark> landmarks;
        void addLandMark( Landmark land);
        unsigned int size();
        void print(int landmarkIndex);
    protected:
    private:
};

#endif // LANDMARKS_H
