#ifndef UTILITIES_H
#define UTILITIES_H
#include <vector>
#include <math.h>
#include <cmath>
#include <random>
#include <string>
#include <iostream>
#include <fstream>
#include <eigen3/Eigen/Dense>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/math/distributions/normal.hpp>
#include <boost/math/distributions/exponential.hpp>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

using namespace std;
using namespace boost::math;
using namespace Eigen;

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

class Utilities
{
    public:
        Utilities();
        virtual ~Utilities();
        //Read a file that has a specific number of dimensions for each attribute
        virtual vector < vector < vector<double> > >  readFile(string CloudSeperator);
        double multivariateNormalPDF(Vector3d instance, Vector3d mu, Matrix3d covar, int dimensionality);
        Eigen::MatrixXd sampleMultivariateNormal( Eigen::RowVector3d mean, Eigen::Matrix3d covar, int numOfSamples, int dimensionality);
        int randcat( vector<double> * vec);
        Eigen::Matrix3d iwishrnd( Matrix3d tau, double nu, int dimensionality);

    protected:
    private:
};

#endif // UTILITIES_H
