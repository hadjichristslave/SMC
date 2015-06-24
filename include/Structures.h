#ifndef STRUCTURES_H
#define STRUCTURES_H
#include <eigen3/Eigen/Dense>
#include <vector>
#include <algorithm>
#include <iostream>

using namespace std;
using namespace Eigen;
namespace Structures {
    const double CRP=.1;
    struct SufficientStatistics{
        std::vector<double> mean        = decltype(mean)(3,0);
        Matrix3d            covar       = MatrixXd::Identity(3,3);
        vector<double>      categorical;
        double              exponential;
        SufficientStatistics(vector<double> pointInstance , int begin){
            double sum = 0;
            for(unsigned int i =begin; i < pointInstance.size(); i++){
                categorical.push_back(pointInstance[i]);
                sum += pointInstance[i];
            }
            if(sum!=0)
                for_each(categorical.begin(), categorical.end(), [&sum] (double &y) mutable { y = y/sum;});
        }
        SufficientStatistics(void){}
        void copy(SufficientStatistics stat){
                mean        = stat.mean;
                covar       = stat.covar;
                categorical = stat.categorical;
                exponential = stat.exponential;
        }
            friend ostream& operator<<(ostream& out ,const SufficientStatistics& rhs){
                out << "====="<< endl;
                out <<" Mean statistics " << endl;
                //out << rhs.mean[0] << "," << rhs.mean[1] << "," << rhs.mean[2] << endl;
                //out << "-" << endl;
                //out << "Covariance statistics" << endl;
                for( int i = 0 ; i < rhs.covar.rows(); i ++){
                    for( int j = 0 ; j < rhs.covar.cols(); j ++){
                        //out << rhs.covar(i,j) ;
                        //out << "," ;
                    }
                }
                out << endl;
                //out << endl;
                for(unsigned int i =0 ; i< rhs.categorical.size();i++)
                    out << rhs.categorical[i] << ",";
                out << endl;
                //out << " Exponent" << rhs.exponential<< endl;
                //out<< rhs.exponential<<endl;
                //out << "==========End of SS print====="<< endl;
                return out;
            }
        };
        struct Params{
            std::vector<double> position = decltype(position)(3,0);
            int  auxiliaryNum , colourBin , colourBins = 3;
            double  crp, del, nu0, kappa0, gamma_alpha0, gamma_beta0 , exp_lambda0 =1;
            RowVector3d  mu0 = RowVector3d::Zero(1,3);
            Matrix3d tau0 = MatrixXd::Identity(3,3);
            RowVectorXd q0 = RowVectorXd::Ones(1,colourBins*colourBins*colourBins);
            Params(void){
                crp = CRP ; del = .7;auxiliaryNum = 20;
                nu0 = 60; kappa0 = .05;
                gamma_alpha0 = 1;  gamma_beta0 = 1;
            }
            Params(double a) : crp(a){
                crp = CRP ; del = .7;auxiliaryNum = 20;
                nu0 = 60; kappa0 = .05;
                gamma_alpha0 = 1;  gamma_beta0 = 1;
            }
            friend ostream& operator<<(ostream& out ,const Params& rhs){
                out << "----------------Params printing-------------------" << endl;
                out << " mu " << rhs.mu0 << endl;
                out << " nu " << rhs.nu0 << endl;
                out << " tau " << rhs.tau0 << endl;
                out << " q0 " << rhs.q0 << endl;
                out << " kappa0 " << rhs.kappa0 << endl;
                out << " gamma alpha " << rhs.gamma_alpha0 << endl;
                out << " gamma beta  " << rhs.gamma_beta0 << endl;
                out << "----------------Params printing end-------------------" << endl;
                return out;
            }
        };
        struct State{
            Params clusterParams;
        };
        struct StateProgression{
            vector< vector< Structures::SufficientStatistics > > stateProg;// Outside vector defines time index
                                                                    // Inside vector defines cluster index
                                                                    // Params define cluster statistics at time t at cluster k
            vector<int>                    assignments;
            vector < vector<int> >         clusterSizes;
            StateProgression(int timeWindow){
                stateProg.resize(timeWindow);
            }
        };
}
#endif // STRUCTURES_H
