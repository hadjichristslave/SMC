#ifndef UTILITIES_H
#define UTILITIES_H
#include <vector>
#include <math.h>
#include <random>
#include <iostream>
#include <fstream>
#include <eigen3/Eigen/Dense>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/math/distributions/normal.hpp>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <opencv2/opencv.hpp>
#include <opencv2/ml/ml.hpp>
#include  "Structures.h"
#include <map>

using namespace std;
using namespace boost::math;
using namespace Eigen;
using namespace cv;

namespace Eigen {
    namespace internal {
        template<typename Scalar>
        struct scalar_normal_dist_op
        {
          static boost::mt19937 rng;    // The uniform pseudo-random algorithm
          mutable boost::normal_distribution<Scalar> norm;  // The gaussian combinator
          EIGEN_EMPTY_STRUCT_CTOR(scalar_normal_dist_op)
          template<typename Index>
          inline const Scalar operator() (Index, Index = 0) const { return norm(rng); }
        };
        template<typename Scalar> boost::mt19937 scalar_normal_dist_op<Scalar>::rng;
        template<typename Scalar>
        struct functor_traits<scalar_normal_dist_op<Scalar> >{
            enum{
                Cost = 50 * NumTraits<Scalar>::MulCost, PacketAccess = false, IsRepeatable = false
            };
        };
    }
}
class Utilities{
    public:
        Utilities();
        virtual ~Utilities();
        //Read a file that has a specific number of dimensions for each attribute
        virtual vector < vector < vector<double> > >  readFile(string CloudSeperator, string filepath);
        double multivariateNormalPDF(Vector3d instance, Vector3d mu, Matrix3d covar, int dimensionality);
        Eigen::MatrixXd sampleMultivariateNormal( Eigen::VectorXd mean, Eigen::MatrixXd covar, int numOfSamples, int dimensionality);
        Eigen::MatrixXd sampleMultinomial(vector<double> probabilities, int samples);
        int randcat( vector<double> * vec);
        Eigen::Matrix3d iwishrnd( Matrix3d tau, double nu, int dimensionality, int df);
        Eigen::VectorXd exprnd(double rate, int samples);
        double exppdf(double x , double lambda);
        double gammarnd(double alpha , double beta);
        RowVectorXd dirrnd(RowVectorXd q0);
        double catpdf(int index , vector<double>  probabilities);
        // Exponential distribution distances
        double Expsquaredhellinger(double lambda1, double lambda2);
        double ExpKLDivergence(double lambda1, double lambda2);
        // Gaussian distribution distances
        double GaussKLDivergence(std::vector<double> mean1, Matrix3d covar1, std::vector<double> mean2, Matrix3d covar2 );
        double Wasserstein(std::vector<double> mean1, Matrix3d covar1, std::vector<double> mean2, Matrix3d covar2 );

        vector<float>  categoricalhistogramCompare( float histA[], float histB[] , int N);
        float categoricalKLDivergence( cv::Mat * mat1, cv::Mat * mat2);
    protected:
    private:
};
#endif // UTILITIES_H
