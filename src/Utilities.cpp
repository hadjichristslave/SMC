#include "Utilities.h"

Utilities::Utilities()
{
    //ctor
}

Utilities::~Utilities()
{
    //dtor
}

vector < vector < vector<double> > >  Utilities::readFile(string CloudSeperator){
    vector< vector< vector<double> > > clouds;
    std::string delimiter = ",";
    std::string line;
    int currentCloud   = -1, currPointIndex = -1;

    ifstream myfile ("/home/panos/Desktop/cloudData/aggregated.csv");

    if (myfile.is_open())
    {
        while ( getline (myfile,line) )
        {

            if(line == CloudSeperator){
                vector< vector< double> > dat;
                clouds.push_back(dat);
                currentCloud++;
                currPointIndex = -1;
                continue;
            }else if(currPointIndex == -1){
                vector< vector< double> > dat;
                clouds.push_back(dat);
                currentCloud++;
                currPointIndex = -1;
            }
            vector<double> newPoint;
            clouds[currentCloud].push_back(newPoint);
            currPointIndex++;
            size_t pos = 0;
            std::string token;
            while ((pos = line.find(delimiter)) != std::string::npos) {
                token = line.substr(0, pos);
                double dat = atof(token.c_str());
                clouds[currentCloud][currPointIndex].push_back(dat);
                line.erase(0, pos + delimiter.length());
            }
            clouds[currentCloud][currPointIndex].push_back( atof(line.c_str()));
        }
        myfile.close();
    }else{
        cout << " could not read data file";
    }


    return clouds;
}
