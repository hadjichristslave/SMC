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
vector< vector< double > >  Landmarks::extractDistances(Landmark * observation , Utilities * ut){
        vector< double > currDist;
        vector< vector< double > > distances;
        if( landmarks.size()==0)
            return distances;
        for(unsigned int i=0;i< landmarks.size();i++){
            distances.push_back(currDist);
            //distances.push_back(currDist);
            Landmark landmark = landmarks[i];
            // Gaussian distances
            distances[i].push_back(ut->Wasserstein(observation->distribution.mean, observation->distribution.covar,\
             landmark.distribution.mean, landmark.distribution.covar));

            distances[i].push_back(ut->GaussKLDivergence(observation->distribution.mean, observation->distribution.covar,\
             landmark.distribution.mean, landmark.distribution.covar));
            // exponential distances
            distances[i].push_back(ut->ExpKLDivergence(observation->distribution.exponential,\
             landmark.distribution.exponential));
            distances[i].push_back(ut->Expsquaredhellinger(observation->distribution.exponential,\
             landmark.distribution.exponential));
            float histogram1[ observation->distribution.categorical.size() ];
            float histogram2[ landmark.distribution.categorical.size() ];
            for (unsigned int o=0;o<landmark.distribution.categorical.size();o++){
                histogram1[o] = observation->distribution.categorical[o]==0?.000001:observation->distribution.categorical[o];
                histogram2[o] = landmark.distribution.categorical[o]==0?.000001:landmark.distribution.categorical[o];
            }
            vector<float> catfdistances = ut->categoricalhistogramCompare(histogram1, histogram2, sizeof(histogram1)/sizeof(float));
            vector<double> catddistances(catfdistances.begin() , catfdistances.end());
            distances[i].insert(distances[i].end(), catddistances.begin(), catddistances.end());
        }
        return distances;
}
