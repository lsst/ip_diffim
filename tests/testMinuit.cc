#include <lsst/fw/MinimizerFunctionBase.h>
#include <lsst/fw/FunctionLibrary.h>
#include <Minuit/FunctionMinimum.h>
#include <Minuit/MnMigrad.h>
#include <Minuit/MnMinos.h>
#include <Minuit/MnPrint.h>
#include <cstdlib>

using namespace std;

int main( int argc, char** argv )
{
    typedef double funcType;
    const unsigned int order = 3;
    const unsigned int npts = 10;
    vector<double> params(order+1);

    boost::shared_ptr<lsst::fw::function::Chebyshev1Function1<funcType> > chebyFuncPtr(
        new lsst::fw::function::Chebyshev1Function1<funcType>(order)
        );

    const unsigned int polyorder = 0;
    boost::shared_ptr<lsst::fw::function::PolynomialFunction1<funcType> > polyFuncPtr(
        new lsst::fw::function::PolynomialFunction1<funcType>(polyorder)
        );

    params[0] = 0;
    params[1] = 0.1;
    params[2] = 0.2;
    params[3] = 0.3;

    chebyFuncPtr->setParameters(params);

    cout << "Input : Chebychev polynomial of the first kind with parameters: ";
    for (unsigned int ii = 0; ii < params.size(); ++ii) {
        cout << params[ii] << " ";
    }
    cout << endl;

    vector<double> measurements(npts);
    vector<double> variances(npts);
    vector<double> positions(npts);

    double x = -1.;
    for (unsigned int i = 0; i < npts; i++, x += 0.2) {
        measurements[i] = (*chebyFuncPtr)(x);
        variances[i] = 0.1;
        positions[i] = x;
    }


    double def = 1.0;
    lsst::fw::function::MinimizerFunctionBase1<funcType> myFcn(measurements, variances, positions, def, chebyFuncPtr);

    // Initialize paramters
    // Name; value; uncertainty
    MnUserParameters upar;
    upar.add("p0", params[0]+double(rand())/RAND_MAX, 0.1);
    upar.add("p1", params[1]+double(rand())/RAND_MAX, 0.1);
    upar.add("p2", params[2]+double(rand())/RAND_MAX, 0.1);
    upar.add("p3", params[3]+double(rand())/RAND_MAX, 0.1);

    std::cout<<"initial parameters: "<<upar<<std::endl;

    MnMigrad migrad(myFcn, upar);
    FunctionMinimum min = migrad();

    MnMinos minos(myFcn, min);

    std::pair<double,double> e0 = minos(0);
    std::pair<double,double> e1 = minos(1);
    std::pair<double,double> e2 = minos(2);
    std::pair<double,double> e3 = minos(3);

    cout << "Best fit:" << endl;
    std::cout<<"par0: "<<min.userState().value("p0")<<" "<<e0.first<<" "<<e0.second<<std::endl;
    std::cout<<"par1: "<<min.userState().value("p1")<<" "<<e1.first<<" "<<e1.second<<std::endl;
    std::cout<<"par2: "<<min.userState().value("p2")<<" "<<e2.first<<" "<<e2.second<<std::endl;
    std::cout<<"par3: "<<min.userState().value("p3")<<" "<<e3.first<<" "<<e3.second<<std::endl;

    // Try fitting for a number with no variation
    // Try fitting for a number with no variation
    // Try fitting for a number with no variation
    // Try fitting for a number with no variation
    MnUserParameters upar2;
    upar2.add("p0", 1, 0.1);
    vector<double> measurements2(npts);
    vector<double> variances2(npts);
    vector<double> positions2(npts);
    x = -1.;
    for (unsigned int i = 0; i < npts; i++, x += 0.2) {
        measurements[i] = 1.;
        variances[i] = 0.1;
        positions[i] = x;
    }
    def = 1.0;
    lsst::fw::function::MinimizerFunctionBase1<funcType> myFcn2(measurements, variances, positions, def, polyFuncPtr);

    MnMigrad migrad2(myFcn2, upar2);
    FunctionMinimum min2 = migrad2();
    MnMinos minos2(myFcn2, min2);

    std::pair<double,double> e02 = minos2(0);

    cout << "Best fit:" << endl;
    std::cout<<"par0: "<<min2.userState().value("p0")<<" "<<e02.first<<" "<<e02.second<<std::endl;
    

    return 0;
}
