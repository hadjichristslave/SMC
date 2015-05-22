#ifndef DBWRAPPER_H
#define DBWRAPPER_H
#include "../sqlite3pp.h"
#include "../sqlite3ppext.h"
#include "SMC.h"
#include "Landmark.h"

using namespace std;
class DBWrapper
{

    public:
        DBWrapper(string db);
        virtual ~DBWrapper();
        void ifNotExistscreateDB();
        void insertLandmark(SMC::SufficientStatistics dist);

        string database;
    protected:
    private:
};

#endif // DBWRAPPER_H
