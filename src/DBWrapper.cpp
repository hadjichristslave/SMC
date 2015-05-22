#include "DBWrapper.h"

DBWrapper::DBWrapper(string db){
    DBWrapper::ifNotExistscreateDB();
    database = db;
}
DBWrapper::~DBWrapper(){}

void DBWrapper::ifNotExistscreateDB(){
    sqlite3pp::database db("landmarks.db");
    db.execute("create table landmarks(id integer primary key autoincrement,\
    x double, y double, z double, \
    covar11 double, covar12 double, covar13 double,\
    covar21 double, covar22 double, covar23 double,\
    covar31 double, covar32 double, covar33 double,\
    exponentialLambda double,\
    bin0 double, bin2 double, bin3 double, bin4 double, bin5 double, bin6 \
    double, bin7 double, bin8 double, bin9 double, bin10 double, bin11 double, bin12 double, \
    bin13 double, bin14 double, bin15 double, bin16 double, bin17 double, bin18 double, bin19 \
    double, bin20 double, bin21 double, bin22 double, bin23 double, bin24 double, bin25 double, bin26 double);");
}
void DBWrapper::insertLandmark(SMC::SufficientStatistics dist){
    sqlite3pp::database db("landmarks.db");
    sqlite3pp::command cmd(
        db, "INSERT INTO landmarks(x , y , z , \
            covar11 , covar12 , covar13 ,\
            covar21 , covar22 , covar23 ,\
            covar31 , covar32 , covar33 ,\
            exponentialLambda ,\
            bin0 , bin2 , bin3 , bin4 , bin5 , \
            bin6, bin7 , bin8 , bin9 , bin10 , bin11 , bin12 , \
            bin13 , bin14 , bin15 , bin16 , bin17 , bin18 , bin19, \
            bin20 , bin21 , bin22 , bin23 , bin24 , bin25 , bin26)\
            VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?\
            , ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)");
        int place = 1;
        cmd.bind(place++, dist.mean[0]);
        cmd.bind(place++, dist.mean[1]);-
        cmd.bind(place++, dist.mean[2]);
        for( int jk = 0 ; jk< dist.covar.rows();jk++)
            for( int kl =0 ; kl<dist.covar.cols();kl++,place++)
                cmd.bind(place,dist.covar(jk,kl));
        //cmd.binder() << dist.exponential;
        for(int jk = 0;jk < dist.categorical.size(); jk ++, place++)
            cmd.bind(place,dist.categorical[jk]);
        cmd.execute();
}
