#include "Utilities.h"
Utilities::Utilities(){}
Utilities::~Utilities(){}

vector < vector < vector<double> > >  Utilities::readFile(string CloudSeperator, string filepath){
    vector< vector< vector<double> > > clouds;
    std::string delimiter = ",";
    std::string line;
    int currentCloud   = -1, currPointIndex = -1;
    ifstream myfile (filepath);
    int lineCount = 0;
    if (myfile.is_open()){
        while ( getline (myfile,line) ){
            if(lineCount ==0){
                lineCount++;
                continue;
            }if(line == CloudSeperator){
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
    }else
        cout << " could not read data file";
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
    for(unsigned int i =0;i<vec->size();i++) vec->at(i) += i>0?vec->at(i-1):0;
    double r = ((double) rand() / (RAND_MAX));
    for(unsigned int i=0;i< vec->size();i++) if( vec->at(i) >= r) return i>=vec->size()?i-1:i;
    return 0;
}
Eigen::MatrixXd Utilities::sampleMultinomial(vector<double> probabilities, int samples){
    // Return samples number of items from a multinomial with probabilities vector defining the probability of samplign each item
    MatrixXd mat(samples, probabilities.size());
    vector<double> cumsum =probabilities;
    for(unsigned int j=1; j<cumsum.size(); ++j)  cumsum[j] += cumsum[j-1];
    for(int i=0;i< samples; i++){
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
    vector<double> norms(probabilities.size());
    for_each(probabilities.begin(), probabilities.end(), [&sum] (double y) mutable { sum +=y; });
    for(unsigned int i =0;i<probabilities.size();i++) norms[i] = probabilities[i]/sum;
    return norms[index];
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
    LLT<Eigen::Matrix3d> cholSolver(tau.inverse());
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
        tempMat = normTransform.transpose() * tempMat;
        Eigen::Matrix3d wishartstuff= tempMat * tempMat.transpose();
        wishartstuff = wishartstuff.array().inverse();
        return wishartstuff;
    }else{
        MatrixXd tempMat(df, dimensionality);
        for(int i=0;i<tempMat.rows();i++)
            for(int j=0;j<tempMat.cols();j++)
                tempMat(i,j) = d(gen);
         MatrixXd tempMat2 = tempMat * normTransform;
         Eigen::MatrixXd tempwishartstuff(tempMat2.transpose() * tempMat2);
         Eigen::Matrix3d wishartstuff(tempwishartstuff.transpose()* tempwishartstuff);
         return wishartstuff.inverse();
    }
}
Eigen::MatrixXd Utilities::sampleMultivariateNormal(VectorXd mean, Eigen::MatrixXd covar, int samples, int dimensionality){
    Eigen::internal::scalar_normal_dist_op<double> randN; // Gaussian functor
    Eigen::internal::scalar_normal_dist_op<double>::rng.seed(1); // Seed the rng
    Eigen::MatrixXd normTransform(dimensionality,dimensionality);
    Eigen::LLT<Eigen::MatrixXd> cholSolver(covar);
    if (cholSolver.info()==Eigen::Success)
            normTransform = cholSolver.matrixL();
    else {
        Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigenSolver(covar);
        normTransform = eigenSolver.eigenvectors()
                   * eigenSolver.eigenvalues().cwiseSqrt().asDiagonal();
    }
    for(int i =0;i< normTransform.rows();i++)
        for(int j =0;j<normTransform.cols();j++)
            normTransform(i,j) = normTransform(i,j)==normTransform(i,j)?normTransform(i,j):0;
    MatrixXd samplez = (normTransform * MatrixXd::NullaryExpr(dimensionality,samples,randN)).colwise() + mean;
    return samplez;
}
//Exponential distribution distances
// Squared hellinger distance as given in wiki
double Utilities::Expsquaredhellinger(double lambda1, double lambda2){
    if(lambda1 ==0 && lambda2 == 0) return 10000000;
    if(lambda1+lambda2<0)          return 10000000;
    return 1-(2*sqrt(lambda1+lambda2))/(lambda1+lambda2);
}
// KL divergence as given in wiki
double Utilities::ExpKLDivergence(double lambda1, double lambda2){
    return log(lambda1) - log(lambda2) + lambda1/lambda2  -1;
}
//KL divergence for gaussians as given in https://tgmstat.wordpress.com/2013/07/10/kullback-leibler-divergence/#ref5
double Utilities::GaussKLDivergence(std::vector<double> tempmean1, Matrix3d covar1, std::vector<double> tempmean2, Matrix3d covar2 ){
    double KL = 0.5;
    Vector3d mean1(3);
    mean1(0) = tempmean1[0];
    mean1(1) = tempmean1[1];
    mean1(2) = tempmean1[2];
    Vector3d mean2(3);
    mean2(0) = tempmean2[0];
    mean2(1) = tempmean2[1];
    mean2(2) = tempmean2[2];
    double firstPart    = log(covar1.determinant()/covar2.determinant());
    MatrixXd second     = (covar1.transpose().array()*covar2.array());
    double secondPart   = second.trace();
    double thirdPart    = (mean2-mean1).transpose()*covar2.inverse()*(mean2-mean1);
    return KL*(firstPart + secondPart+ thirdPart);
}
double Utilities::Wasserstein(std::vector<double> tempmean1, Matrix3d covar1, std::vector<double> tempmean2, Matrix3d covar2 ){
    Vector3d mean1(3);
    mean1(0) = tempmean1[0];
    mean1(1) = tempmean1[1];
    mean1(2) = tempmean1[2];
    Vector3d mean2(3);
    mean2(0) = tempmean2[0];
    mean2(1) = tempmean2[1];
    mean2(2) = tempmean2[2];
    Vector3d diff = mean1 - mean2;
    MatrixXd first  = covar1 + covar2;
    MatrixXd second = covar2.array().sqrt() * covar1.array() * covar2.array().sqrt();
    second          = second.array().sqrt();
    second          = second.array() *2;
    first           = first - second;
    return diff.array().square().sum() + first.trace();
}
std::vector<float>  Utilities::categoricalhistogramCompare( float histA[], float histB[], int N){
     float maxA = *std::max_element(histA, histA+N);
     float maxB = *std::max_element(histB, histB+N);
     cv::Mat M1 = cv::Mat(1,N, cv::DataType<float>::type , histA);
     cv::Mat M2 = cv::Mat(1,N, cv::DataType<float>::type , histB);
     int histSize = N;
     float rangeA[] = {0, maxA+1};//ranges are exclusive hence + 1
     float rangeB[] = {0, maxB+1};// see above
     const float* histRangeA = {rangeA};
     const float* histRangeB = {rangeB};
     bool uniform = true;
     bool accumulate = false;
     cv::Mat a1_hist, a2_hist;
     // normalization means SQRT( sum(component*component)) = 1A
     cv::calcHist(&M1, 1, 0, cv::Mat(), a1_hist, 1, &histSize, &histRangeA, uniform, accumulate );
     cv::calcHist(&M2, 1, 0, cv::Mat(), a2_hist, 1, &histSize, &histRangeB, uniform, accumulate );
     normalize(a1_hist, a1_hist,  0, 1, CV_MINMAX);
     normalize(a2_hist, a2_hist,  0, 1, CV_MINMAX);
     cv::Mat sig1(N, 2, cv::DataType<float>::type);
     cv::Mat sig2(N, 2, cv::DataType<float>::type);
     for(int i=0;i<histSize;i++){
         float binval = a1_hist.at<float>(i);
         sig1.at< float >(i, 0) = binval;
         sig1.at< float >(i, 1) = i;
         binval = a2_hist.at< float>(i);
         sig2.at< float >(i, 0) = binval;
         sig2.at< float >(i, 1) = i;
     }
     float emd           = cv::EMD(sig1, sig2, CV_DIST_L2);
     float compar_hell   = (float)cv::compareHist(a1_hist, a2_hist, CV_COMP_HELLINGER );
     // WARNING!!! KLDivergence changegs the values of the histograms so it should be called last.
     float kld =  categoricalKLDivergence( &a1_hist, &a2_hist);
     vector<float> distances;
     distances.push_back(kld);
     distances.push_back(emd);
     distances.push_back(compar_hell);
     return distances;
}
float Utilities::categoricalKLDivergence( cv::Mat * mat1, cv::Mat * mat2){
     float sum1 = 0,sum2 = 0;
     for(int i=0;i<mat1->rows;i++){
        sum1 += mat1->at<float>(i,0);
        sum2 += mat2->at<float>(i,0);
     }
     for(int i=0;i<mat1->rows;i++){
         mat1->at<float>(i,0) /= sum1;
         mat2->at<float>(i,0) /= sum2;
     }
     float result = 0.;
     for(int i=0;i< mat1->rows;i++)
         if(  mat1->at<float>(0,i) !=0 ){
            float ratio = mat1->at<float>(i,0)/ mat2->at<float>(0,i);
             if(ratio>0 && ratio  != std::numeric_limits<float>::infinity() )
                 result += mat1->at<float>(i,0) * log(ratio);
         }
     return result;
 }
