#include <iostream>
#include <vector>
#include <algorithm>

#ifndef POINT_H
#define POINT_H

using namespace std;
class Point{

    public:
        int x , y , z ,  colourBin ;        // Define x,y,z,r,g,b and point id
        double KLDist, EMD;
        inline Point();
        const inline void add(const double val);
        const inline void printDistances(const vector<double> & v);
    private:



};

inline Point::Point(){
    x = 0; y = 0 ; z = 0, colourBin =0;
}
const inline void Point::printDistances(const vector<double>& v){
    std::for_each(v.begin(), v.end(), [](double v) { cout << v << "," <<endl; });
}
#endif // POINT_H
