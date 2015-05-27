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
            distances.push_back(currDist);
            Landmark landmark = landmarks[i];
            // Gaussian distances
            distances[2*i].push_back(ut->Wasserstein(observation->distribution.mean, observation->distribution.covar,\
             landmark.distribution.mean, landmark.distribution.covar));
            distances[2*i+1].push_back(ut->Wasserstein(observation->distribution.mean, observation->distribution.covar,\
             landmark.distribution.mean, landmark.distribution.covar)+40);
            distances[2*i].push_back(ut->GaussKLDivergence(observation->distribution.mean, observation->distribution.covar,\
             landmark.distribution.mean, landmark.distribution.covar));
            distances[2*i+1].push_back(ut->GaussKLDivergence(observation->distribution.mean, observation->distribution.covar,\
             landmark.distribution.mean, landmark.distribution.covar)+40);
            // exponential distances
            distances[2*i].push_back(ut->ExpKLDivergence(observation->distribution.exponential,\
             landmark.distribution.exponential));
            distances[2*i+1].push_back(ut->ExpKLDivergence(observation->distribution.exponential,\
             landmark.distribution.exponential));
            distances[2*i].push_back(ut->Expsquaredhellinger(observation->distribution.exponential,\
             landmark.distribution.exponential));
            distances[2*i+1].push_back(ut->Expsquaredhellinger(observation->distribution.exponential,\
             landmark.distribution.exponential));
            //categorical distances
            float histogram1[ observation->distribution.categorical.size() ];
            float histogram2[ landmark.distribution.categorical.size() ];
            for (unsigned int o=0;o<landmark.distribution.categorical.size();o++){
                histogram1[o] = observation->distribution.categorical[o]==0?.000001:observation->distribution.categorical[o];
                histogram2[o] = landmark.distribution.categorical[o]==0?.000001:landmark.distribution.categorical[o];
            }
            vector<float> catfdistances = ut->categoricalhistogramCompare(histogram1, histogram2, sizeof(histogram1)/sizeof(float));
            vector<double> catddistances(catfdistances.begin() , catfdistances.end());
            distances[2*i].insert(distances[2*i].end(), catddistances.begin(), catddistances.end());
            distances[2*i+1].insert(distances[2*i+1].end(), catddistances.begin(), catddistances.end());
        }
        normalizeFeatures(&distances);
        return distances;
}
//NOrmalize gaussian distances using formula
// feature = feature-dev/mean
void Landmarks::normalizeFeatures(vector <vector <double > > * features){
    double WassersteinSum =0;
    double gausKLSum      =0;

    for (int i = 0;i< features->size();i++){
        WassersteinSum += features->at(i)[0];
        gausKLSum      += features->at(i)[1];
    }
    double WassersteinAvg  = WassersteinSum/features->size();
    double gausKLAvg       = gausKLSum/features->size();

    double WassersteinStd = 0;
    double gausKLStd      = 0;

    for (int i =0;i< features->size();i++){
        WassersteinStd += pow(features->at(i)[0]-WassersteinAvg,2);
        gausKLStd      += pow(features->at(i)[1]-gausKLAvg     ,2);
    }
    for (int i =0;i< features->size();i++){
        features->at(i)[0] = (features->at(i)[0]-WassersteinAvg)/WassersteinStd;
        features->at(i)[1] = (features->at(i)[1]-gausKLAvg)/gausKLStd;
    }
}
