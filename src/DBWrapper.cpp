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
    db.execute("create table distances(id integer primary key autoincrement,\
     feature1 double, feature2 double ,feature3 double, feature4 double,\
     feature5 double, feature6 double, feature7 double, label integer)");
}
void DBWrapper::insertLandmark(SufficientStatistics * dist){
    sqlite3pp::database db(database.c_str());
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
        cmd.bind(place++, dist->mean[0]);
        cmd.bind(place++, dist->mean[1]);
        cmd.bind(place++, dist->mean[2]);
        for( int jk = 0 ; jk< dist->covar.rows();jk++)
            for( int kl =0 ; kl<dist->covar.cols();kl++,place++)
                cmd.bind(place,dist->covar(jk,kl));
        cmd.bind(place++ ,dist->exponential);
        for(unsigned int jk = 0;jk < dist->categorical.size(); jk ++, place++)
            cmd.bind(place,dist->categorical[jk]);
        cmd.execute();
}
//Retrieve landmark with ID id. Id's are considered unique
Landmark DBWrapper::getLandmark(int LandId){
    sqlite3pp::database db(database.c_str());
    string query = "SELECT * FROM landmarks where id="+SSTR(LandId);
    sqlite3pp::query qry(db, query.c_str());
    for(sqlite3pp::query::iterator ij = qry.begin();ij != qry.end();++ij){
        int id=-1;
        SufficientStatistics stats;
        stats.categorical.resize(27);
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
        Landmark land(id, stats);
        return land;
    }
    //Case landmark id does not exist
    Landmark land;
    return  land;
}
Landmarks DBWrapper::getCurrentLandmarks(){
    Landmarks landmarks;
    sqlite3pp::database db(database.c_str());
    sqlite3pp::query qry(db, "SELECT * FROM landmarks");
    for(sqlite3pp::query::iterator ij = qry.begin();ij != qry.end();++ij){
        int id=-1;
        SufficientStatistics stats;
        stats.categorical.resize(27);
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
        Landmark land(id, stats);
        landmarks.addLandMark(land);
    }
    return landmarks;
}
vector< vector< double > > DBWrapper::getTrainingSet(){
    //Retrieve 20 and 20 samples for the training
    vector< vector< double > >  trainingSet;
    Landmarks landmarks;
    sqlite3pp::database db(database.c_str());
    sqlite3pp::query qry(db,  "SELECT * FROM distances2 where label=1 limit 20");
    sqlite3pp::query qry2(db, "SELECT * FROM distances2 where label=0 limit 20");
    for(sqlite3pp::query::iterator ij = qry.begin();ij != qry.end();++ij){
        int id=-1;
        vector<double> currentDist(8);
        std::tie(id, currentDist[0], currentDist[1], currentDist[2])=\
        (*ij).get_columns<int , double, double, double>(0,1,2,3);
        std::tie(id, currentDist[3], currentDist[4], currentDist[5])=\
        (*ij).get_columns<int , double, double, double>(0,4,5,6);
        std::tie(id, currentDist[6], currentDist[7])=\
        (*ij).get_columns<int , double, double>(0,7,8);
        trainingSet.push_back(currentDist);
    }
    for(sqlite3pp::query::iterator ij = qry2.begin();ij != qry2.end();++ij){
        int id=-1;
        vector<double> currentDist(8);
        std::tie(id, currentDist[0], currentDist[1], currentDist[2])=\
        (*ij).get_columns<int , double, double, double>(0,1,2,3);
        std::tie(id, currentDist[3], currentDist[4], currentDist[5])=\
        (*ij).get_columns<int , double, double, double>(0,4,5,6);
        std::tie(id, currentDist[6], currentDist[7])=\
        (*ij).get_columns<int , double, double>(0,7,8);
        trainingSet.push_back(currentDist);
    }
    return trainingSet;
}


void DBWrapper::insertLabeledDistances(vector< vector< double > > distanceFeatures, int positiveLabel){
    sqlite3pp::database db(database.c_str());
    for( unsigned int i = 0; i < distanceFeatures.size();i++){
        sqlite3pp::command cmd(
            db, "INSERT INTO distances2(feature1,\
                feature2, feature3, feature4, feature5,\
                feature6, feature7, label)\
                VALUES (?, ?, ?, ?, ?, ?, ?, ?)");
            int place = 1;
            for(unsigned int j = 0; j < distanceFeatures[i].size();j++)
                cmd.bind(place++, distanceFeatures[i][j]);
            int label = (2*positiveLabel==i || 2*positiveLabel==i-1 )?1:0;
            cmd.bind(place++ , label);
            cmd.execute();
    }
}
