#ifndef DBWRAPPER_H
#define DBWRAPPER_H
#include "sqlite3pp.h"
#include "sqlite3ppext.h"
#include "SMC.h"
#include "Landmark.h"
#include "Landmarks.h"
#include "Structures.h"

#define SSTR( x ) dynamic_cast< std::ostringstream & >( \
        ( std::ostringstream() << std::dec << x ) ).str()

using namespace std;
using namespace Structures;
class DBWrapper
{

    public:
        DBWrapper(string db);
        virtual ~DBWrapper();
        void ifNotExistscreateDB();
        void insertLandmark(SufficientStatistics * dist);
        Landmark getLandmark(int LandId);
        Landmarks getCurrentLandmarks();

        string database;
    protected:
    private:
};

#endif // DBWRAPPER_H
