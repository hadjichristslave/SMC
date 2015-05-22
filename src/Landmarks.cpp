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
void Landmarks::extractDistances(Landmark landmark1 , Landmark landmark2 , Utilities ut){
            cout << "wasser " << endl;
        cout << ut.Wasserstein(landmark1.distribution.mean, landmark1.distribution.covar, landmark2.distribution.mean, landmark2.distribution.covar) << endl;
        cout << "gauskld" << endl;
        cout << ut.GaussKLDivergence(landmark1.distribution.mean, landmark1.distribution.covar, landmark2.distribution.mean, landmark2.distribution.covar) << endl;
        cout << "expkld " << endl;
        cout << ut.ExpKLDivergence(landmark1.distribution.exponential, landmark2.distribution.exponential) << endl;
        cout << " exphellinger " << endl;
        cout << ut.Expsquaredhellinger(landmark1.distribution.exponential, landmark2.distribution.exponential) << endl;
        float histogram1[ landmark1.distribution.categorical.size() ];
        float histogram2[ landmark2.distribution.categorical.size() ];
        for (unsigned int o=0;o<landmark2.distribution.categorical.size();o++){
            histogram1[o] = landmark1.distribution.categorical[o]==0?.000001:landmark1.distribution.categorical[o];
            histogram2[o] = landmark2.distribution.categorical[o]==0?.000001:landmark1.distribution.categorical[o];
        }
        vector<float> distances = ut.categoricalhistogramCompare(histogram1, histogram2, sizeof(histogram1)/sizeof(float));
        for_each(distances.begin(), distances.end(), [] (float y) { cout << y << ",";});
        cout << endl;
}
