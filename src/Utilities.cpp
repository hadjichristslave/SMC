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

    ifstream myfile ("/home/panos/Desktop/cloudData/subset.csv");
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
                currPointIndex = -1;
                continue;
            }if(currPointIndex == -1){
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
                if (line=="") break;
            }
            clouds[currentCloud][currPointIndex].push_back( atof(line.c_str()));
        }
        myfile.close();
    }else{
        cout << " could not read data file";
    }

    return clouds;
}

double Utilities::multivariateNormalPDF( Vector3d instance, Vector3d mu, Matrix3d covar , int dimensionality){
    Eigen::Matrix3d normTransform(dimensionality, dimensionality);
    Eigen::LLT<Eigen::Matrix3d> cholSolver(covar);
    normTransform = cholSolver.matrixL();
    RowVector3d difference = instance - mu;
    RowVector3d tempDif  = difference * normTransform.transpose().inverse();
    double pdf = pow( (2*M_PI) , ((double)-dimensionality/2) ) * \
    exp(-tempDif.array().square().sum()/2)/normTransform.transpose().diagonal().prod();
    return pdf==0?.00000000000000001:pdf;
}
int Utilities::randcat( vector<double> * vec){
    // GEt cumsum
    for(int i =0;i<vec->size();i++) vec->at(i) += i>0?vec->at(i-1):0;
    double r = ((double) rand() / (RAND_MAX));
    for(int i=0;i< vec->size();i++) if( vec->at(i) >= r) return i;
}
Eigen::MatrixXd Utilities::sampleMultinomial(vector<double> probabilities, int samples){
    // Return samples number of items from a multinomial with probabilities vector defining the probability of samplign each item
    MatrixXd mat(samples, probabilities.size());
    vector<double> cumsum =probabilities;
    for(int j=1; j<cumsum.size(); ++j)  cumsum[j] += cumsum[j-1];

    for(int i=0;i< samples; i++){
        double r = ((double) rand() / (RAND_MAX));
        int sample = randcat( & cumsum);
        Eigen::RowVectorXd tempVec;
        tempVec.resize(cumsum.size());
        tempVec.setZero();
        tempVec(sample)++;
        mat.row(i) = tempVec;
    }
    return mat;
}
Eigen::VectorXd Utilities::exprnd(double rate , int samples){
    VectorXd vec(samples);
    const gsl_rng_type * T;
    gsl_rng * r;
    gsl_rng_env_setup();
    T = gsl_rng_default;
    r = gsl_rng_alloc(T);
    for(int i =0 ; i<samples;i++)
        vec(i) = gsl_ran_exponential(r, rate);
    return vec;
}
double Utilities::exppdf(double x , double lambda){
    return gsl_ran_exponential_pdf(x , lambda);
}
double Utilities::catpdf(int index , vector<double>  probabilities){
    double sum = 0;
    for_each(probabilities.begin(), probabilities.end(), [&sum] (double y) mutable { sum +=y; });
    for_each(probabilities.begin(), probabilities.end(), [&sum] (double y) mutable { y /= sum; });
    return probabilities[index];
}
double Utilities::gammarnd(double alpha, double beta){
    const gsl_rng_type * T;
    gsl_rng * r;
    gsl_rng_env_setup();
    T = gsl_rng_default;
    r = gsl_rng_alloc(T);
    return gsl_ran_gamma(r , alpha, beta);
}
RowVectorXd Utilities::dirrnd(RowVectorXd q0){
    RowVectorXd resp(q0.cols());
    double sum = 0;
    for( int i =0;i<q0.cols();i++){
        resp(i)  = gammarnd(q0(i),1);
        sum += resp(i);
    }for(int i =0;i<q0.cols();i++) resp(i)  /= sum;
    return resp;
}
Matrix3d Utilities::iwishrnd( Matrix3d tau, double nu, int dimensionality , int df){
    Matrix3d normTransform(dimensionality, dimensionality);
    LLT<Eigen::Matrix3d> cholSolver(tau);
    if (cholSolver.info()==Eigen::Success)
        normTransform = cholSolver.matrixL();
    else{
        SelfAdjointEigenSolver<Eigen::Matrix3d> eigenSolver(tau);
        normTransform = eigenSolver.eigenvectors()* eigenSolver.eigenvalues().cwiseSqrt().asDiagonal();
     }
    int sizeOfMat = normTransform.rows();
    std::random_device rd;
    std::mt19937 gen(rd());
    std::normal_distribution<> d(0,1);
    if(df ==0.0){
        const gsl_rng_type * T;
        gsl_rng * r;
        gsl_rng_env_setup();
        T = gsl_rng_default;
        r = gsl_rng_alloc(T);
        RowVector3d chi2rand(dimensionality);
        chi2rand << gsl_ran_chisq(r , nu) , gsl_ran_chisq(r , nu-1) , gsl_ran_chisq(r , nu-2);
        Matrix3d tempMat = chi2rand.cwiseSqrt().asDiagonal();
        for(int i=0;i<tempMat.rows();i++)
            for(int j=0;j<tempMat.cols();j++)
                if(i>j) tempMat(i,j) = d(gen);
        tempMat = normTransform * tempMat;
        Eigen::Matrix3d wishartstuff= tempMat * tempMat.transpose();
        wishartstuff = wishartstuff.array().inverse();
        return wishartstuff;
    }else{
        MatrixXd tempMat(df, dimensionality);
        for(int i=0;i<tempMat.rows();i++)
            for(int j=0;j<tempMat.cols();j++)
                if(i>j) tempMat(i,j) = d(gen);

         Eigen::MatrixXd tempwishartstuff(tempMat * normTransform);
         Eigen::Matrix3d wishartstuff(tempwishartstuff.transpose()* tempwishartstuff) ;
         return wishartstuff.inverse();
    }
}

Eigen::MatrixXd Utilities::sampleMultivariateNormal(VectorXd mean, Eigen::MatrixXd covar, int samples, int dimensionality){
    Eigen::internal::scalar_normal_dist_op<double> randN; // Gaussian functor
    Eigen::internal::scalar_normal_dist_op<double>::rng.seed(1); // Seed the rng
    Eigen::MatrixXd normTransform(dimensionality,dimensionality);
    Eigen::LLT<Eigen::MatrixXd> cholSolver(covar);
    normTransform = cholSolver.matrixL();
    MatrixXd samplez = (normTransform * MatrixXd::NullaryExpr(dimensionality,samples,randN)).colwise() + mean;
    return samplez;
}
