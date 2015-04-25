#include <iostream>
#include <vector>
#include <algorithm>
#include "Point.h"

#ifndef POINTCLOUD_H
#define POINTCLOUD_H

using namespace std;

class PointCloud{

    public:
        vector<Point> Points;  // Point cloud is a vector of points
        int timeId;            // Define time index

};

#endif // POINTCLOUD_H
