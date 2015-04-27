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
    int lineCount = 0;
    if (myfile.is_open())
    {
        while ( getline (myfile,line) )
        {
            if(lineCount ==0){
                lineCount++;
                continue;
            }

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

double Utilities::sampleMultivariateNormal( Vector3d instance, Vector3d mu, Matrix3d covar , int dimensionality){

  Eigen::Matrix3d normTransform(dimensionality, dimensionality);
  Eigen::LLT<Eigen::Matrix3d> cholSolver(covar);

  // We can only use the cholesky decomposition if
  // the covariance matrix is symmetric, pos-definite.
  // But a covariance matrix might be pos-semi-definite.
  // In that case, we'll go to an EigenSolver

  if (cholSolver.info()==Eigen::Success) {
    // Use cholesky solver
    normTransform = cholSolver.matrixL();
  } else {
    // Use eigen solver
    Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> eigenSolver(covar);
    normTransform = eigenSolver.eigenvectors()
                   * eigenSolver.eigenvalues().cwiseSqrt().asDiagonal();
  }
    RowVector3d difference = instance - mu;
    RowVector3d tempDif  = difference * normTransform.inverse();
    double pdf = pow( (2*M_PI) , ((double)-dimensionality/2) ) * exp(tempDif.array().square().sum()/2)/(double)normTransform.eigenvalues().sum().real();
    return pdf;
}
int Utilities::randcat( vector<double> * vec){
    // GEt cumsum
    for(int i =0;i<vec->size();i++) vec->at(i) += i>0?vec->at(i-1):0;
    double r = ((double) rand() / (RAND_MAX)) + 1;
    for(int i=0;i< vec->size();i++) if( vec->at(i) < r) return i;
}
Matrix3d Utilities::iwishrnd( Matrix3d tau, double nu, int dimensionality){
    Matrix3d normTransform(dimensionality, dimensionality);
    LLT<Eigen::Matrix3d> cholSolver(tau);
    if (cholSolver.info()==Eigen::Success)
        normTransform = cholSolver.matrixL();
    else{
     SelfAdjointEigenSolver<Eigen::Matrix3d> eigenSolver(tau);
     normTransform = eigenSolver.eigenvectors()
                    * eigenSolver.eigenvalues().cwiseSqrt().asDiagonal();}

    int sizeOfMat = normTransform.rows();
    const gsl_rng_type * T;
    gsl_rng * r;
    gsl_rng_env_setup();

    T = gsl_rng_default;
    r = gsl_rng_alloc(T);

    RowVector3d chi2rand(dimensionality);
    chi2rand << gsl_ran_chisq(r , nu) , gsl_ran_chisq(r , nu-1) , gsl_ran_chisq(r , nu-2);
    Matrix3d tempMat = chi2rand.cwiseSqrt().asDiagonal();

    std::random_device rd;
    std::mt19937 gen(rd());
    std::normal_distribution<> d(0,1);
    for( int i=0;i<tempMat.rows();i++)
        for(int j=0;j<tempMat.cols();j++)
            if(i>j) tempMat(i,j) = d(gen);
    tempMat = normTransform * tempMat;
    Eigen::Matrix3d wishartstuff= tempMat * tempMat.transpose();
    wishartstuff = wishartstuff.array().inverse();
    return wishartstuff;
}
