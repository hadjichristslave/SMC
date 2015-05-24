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
std::vector<double>  Landmarks::extractDistances(Landmark landmark1 , Landmark landmark2 , Utilities ut){
        vector<double> distances;
        // Gaussian distances
        distances.push_back(ut.Wasserstein(landmark1.distribution.mean, landmark1.distribution.covar, landmark2.distribution.mean, landmark2.distribution.covar));
        distances.push_back(ut.GaussKLDivergence(landmark1.distribution.mean, landmark1.distribution.covar, landmark2.distribution.mean, landmark2.distribution.covar));
        // exponential distances
        distances.push_back(ut.ExpKLDivergence(landmark1.distribution.exponential, landmark2.distribution.exponential));
        distances.push_back(ut.Expsquaredhellinger(landmark1.distribution.exponential, landmark2.distribution.exponential));
        //categorical distances
        float histogram1[ landmark1.distribution.categorical.size() ];
        float histogram2[ landmark2.distribution.categorical.size() ];
        for (unsigned int o=0;o<landmark2.distribution.categorical.size();o++){
            histogram1[o] = landmark1.distribution.categorical[o]==0?.000001:landmark1.distribution.categorical[o];
            histogram2[o] = landmark2.distribution.categorical[o]==0?.000001:landmark1.distribution.categorical[o];
        }
        vector<float> catfdistances = ut.categoricalhistogramCompare(histogram1, histogram2, sizeof(histogram1)/sizeof(float));
        vector<double> catddistances(catfdistances.begin() , catfdistances.end());
        distances.insert(distances.end(), catddistances.begin(), catddistances.end());
        return distances;
}
