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
Landmark DBWrapper::getLandmark(int LandId){
    sqlite3pp::database db("landmarks.db");
    string query = "SELECT * FROM landmarks where id="+SSTR(LandId);
    sqlite3pp::query qry(db, query.c_str());
    if(qry.column_count() ==0){
        Landmark land;
        return  land;
    }
    for (int i = 0; i < qry.column_count(); ++i) {
        int id;
        double x,y,z, covar11 , covar12 , covar13 ,\
        covar21 , covar22 , covar23 ,\
        covar31 , covar32 , covar33 ,\
        exponentialLambda ,\
        bin0 , bin1, bin2 , bin3 , bin4 , bin5 , \
        bin6, bin7 , bin8 , bin9 , bin10 , bin11 , bin12 , \
        bin13 , bin14 , bin15 , bin16 , bin17 , bin18 , bin19, \
        bin20 , bin21 , bin22 , bin23 , bin24 , bin25 , bin26;

        SMC::SufficientStatistics stats;
        stats.categorical.resize(27);
        for(sqlite3pp::query::iterator ij = qry.begin();ij != qry.end();++ij){

            std::tie(id, stats.mean[0],stats.mean[1],stats.mean[2])                             =\
            (*ij).get_columns<int , double, double, double>(0,1,2,3);
            std::tie(id, stats.covar(0,0), stats.covar(0,1), stats.covar(0,2))                  =\
            (*ij).get_columns<int , double, double, double>(0,4,5,6);
            std::tie(id, stats.covar(1,0), stats.covar(1,1), stats.covar(2,2))                  =\
            (*ij).get_columns<int , double, double, double>(0,7,8,9);
            std::tie(id, stats.covar(2,0), stats.covar(2,1), stats.covar(2,2))                  =\
            (*ij).get_columns<int , double, double, double>(0,10,11,12);
            std::tie(id, stats.exponential, stats.categorical[0], stats.categorical[1])         =\
            (*ij).get_columns<int , double, double, double>(0,13,14,15);
            std::tie(id, stats.categorical[2], stats.categorical[3], stats.categorical[4])      =\
            (*ij).get_columns<int , double, double, double>(0,16,17,18);
            std::tie(id, stats.categorical[5], stats.categorical[6], stats.categorical[7])      =\
            (*ij).get_columns<int , double, double, double>(0,19,20,21);
            std::tie(id, stats.categorical[8], stats.categorical[9], stats.categorical[10])     =\
            (*ij).get_columns<int , double, double, double>(0,22,23,24);
            std::tie(id, stats.categorical[11], stats.categorical[12], stats.categorical[13])   =\
            (*ij).get_columns<int , double, double, double>(0,25,26,27);
            std::tie(id, stats.categorical[14], stats.categorical[15], stats.categorical[16])   =\
            (*ij).get_columns<int , double, double, double>(0,28,29,30);
            std::tie(id, stats.categorical[17], stats.categorical[18], stats.categorical[19])   =\
            (*ij).get_columns<int , double, double, double>(0,31,32,33);
            std::tie(id, stats.categorical[20], stats.categorical[21], stats.categorical[22])   =\
            (*ij).get_columns<int , double, double, double>(0,34,35,36);
            std::tie(id, stats.categorical[23], stats.categorical[24], stats.categorical[25])   =\
            (*ij).get_columns<int , double, double, double>(0,37,38,39);
            std::tie(id, stats.categorical[26])                                                 =\
            (*ij).get_columns<int , double>(0,40);
        }
        Landmark land(id, stats);
        return land;
    }
}
Landmarks DBWrapper::getCurrentLandmarks(){
    Landmarks landmarks;
    sqlite3pp::database db("landmarks.db");
    sqlite3pp::query qry(db, "SELECT * FROM landmarks");

    for (int i = 0; i < qry.column_count(); ++i) {
        int id;
        double x,y,z, covar11 , covar12 , covar13 ,\
        covar21 , covar22 , covar23 ,\
        covar31 , covar32 , covar33 ,\
        exponentialLambda ,\
        bin0 , bin1, bin2 , bin3 , bin4 , bin5 , \
        bin6, bin7 , bin8 , bin9 , bin10 , bin11 , bin12 , \
        bin13 , bin14 , bin15 , bin16 , bin17 , bin18 , bin19, \
        bin20 , bin21 , bin22 , bin23 , bin24 , bin25 , bin26;

        SMC::SufficientStatistics stats;
        stats.categorical.resize(27);
        for(sqlite3pp::query::iterator ij = qry.begin();ij != qry.end();++ij){

            std::tie(id, stats.mean[0],stats.mean[1],stats.mean[2])                             =\
            (*ij).get_columns<int , double, double, double>(0,1,2,3);
            std::tie(id, stats.covar(0,0), stats.covar(0,1), stats.covar(0,2))                  =\
            (*ij).get_columns<int , double, double, double>(0,4,5,6);
            std::tie(id, stats.covar(1,0), stats.covar(1,1), stats.covar(2,2))                  =\
            (*ij).get_columns<int , double, double, double>(0,7,8,9);
            std::tie(id, stats.covar(2,0), stats.covar(2,1), stats.covar(2,2))                  =\
            (*ij).get_columns<int , double, double, double>(0,10,11,12);
            std::tie(id, stats.exponential, stats.categorical[0], stats.categorical[1])         =\
            (*ij).get_columns<int , double, double, double>(0,13,14,15);
            std::tie(id, stats.categorical[2], stats.categorical[3], stats.categorical[4])      =\
            (*ij).get_columns<int , double, double, double>(0,16,17,18);
            std::tie(id, stats.categorical[5], stats.categorical[6], stats.categorical[7])      =\
            (*ij).get_columns<int , double, double, double>(0,19,20,21);
            std::tie(id, stats.categorical[8], stats.categorical[9], stats.categorical[10])     =\
            (*ij).get_columns<int , double, double, double>(0,22,23,24);
            std::tie(id, stats.categorical[11], stats.categorical[12], stats.categorical[13])   =\
            (*ij).get_columns<int , double, double, double>(0,25,26,27);
            std::tie(id, stats.categorical[14], stats.categorical[15], stats.categorical[16])   =\
            (*ij).get_columns<int , double, double, double>(0,28,29,30);
            std::tie(id, stats.categorical[17], stats.categorical[18], stats.categorical[19])   =\
            (*ij).get_columns<int , double, double, double>(0,31,32,33);
            std::tie(id, stats.categorical[20], stats.categorical[21], stats.categorical[22])   =\
            (*ij).get_columns<int , double, double, double>(0,34,35,36);
            std::tie(id, stats.categorical[23], stats.categorical[24], stats.categorical[25])   =\
            (*ij).get_columns<int , double, double, double>(0,37,38,39);
            std::tie(id, stats.categorical[26])                                                 =\
            (*ij).get_columns<int , double>(0,40);
        }
        cout << stats;
        Landmark land(id, stats);
        landmarks.addLandMark(land);
        cout << endl;
    }
    return landmarks;
}
