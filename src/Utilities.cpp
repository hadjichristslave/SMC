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

Eigen::MatrixXd Utilities::sampleMultivariateNormal(Eigen::VectorXd mean , Eigen::MatrixXd covar , int dimensionality ){
  Eigen::internal::scalar_normal_dist_op<double> randN; // Gaussian functor
  Eigen::internal::scalar_normal_dist_op<double>::rng.seed(1); // Seed the rng


  Eigen::MatrixXd normTransform(dimensionality, dimensionality);
  Eigen::LLT<Eigen::MatrixXd> cholSolver(covar);

  // We can only use the cholesky decomposition if
  // the covariance matrix is symmetric, pos-definite.
  // But a covariance matrix might be pos-semi-definite.
  // In that case, we'll go to an EigenSolver
  if (cholSolver.info()==Eigen::Success) {
    // Use cholesky solver
    normTransform = cholSolver.matrixL();
  } else {
    // Use eigen solver
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigenSolver(covar);
    normTransform = eigenSolver.eigenvectors()
                   * eigenSolver.eigenvalues().cwiseSqrt().asDiagonal();
  }

  Eigen::MatrixXd samples = (normTransform
                           * Eigen::MatrixXd::NullaryExpr(size,1,randN)).colwise()
                           + mean;
    std::cout << "Mean\n" << mean << std::endl;
    std::cout << "Covar\n" << covar << std::endl;
    std::cout << "Samples\n" << samples << std::endl;
    return samples;
}
