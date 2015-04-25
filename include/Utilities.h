#ifndef UTILITIES_H
#define UTILITIES_H
#include <vector>
#include <string>
#include <iostream>
#include <fstream>

using namespace std;

class Utilities
{
    public:
        Utilities();
        virtual ~Utilities();
        //Read a file that has a specific number of dimensions for each attribute
        virtual vector < vector < vector<double> > >  readFile(string CloudSeperator);
    protected:
    private:
};

#endif // UTILITIES_H
