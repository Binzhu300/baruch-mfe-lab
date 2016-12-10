#include <fd.h>
#include <evaluator.h>
#include <option_value.h>
#include <black_scholes.h>
#include <heat_pde.h>
#include <tuple>
#include <vector>
#include <cassert>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <cassert>

FiniteDifference::FiniteDifference( double expiry,
                                    double spot,
                                    double strike,
                                    double rate,
                                    double div,
                                    double vol,
                                    Evaluator & terminalCondition,
                                    Evaluator & leftBoundaryCondition,
                                    Evaluator & rightBoundaryCondition,
                                    const Evaluator & earlyExercisePremium )
                                  : d_expiry(expiry),
                                    d_spot(spot),
                                    d_strike(strike),
                                    d_rate(rate),
                                    d_div(div),
                                    d_vol(vol),
                                    d_f(&terminalCondition),
                                    d_gl(&leftBoundaryCondition),
                                    d_gr(&rightBoundaryCondition),
                                    d_prem(&earlyExercisePremium)
{
    // Calculate the Black-Scholes value
    // in case of a vanilla call or put.
    d_isPut = ((*d_f)(-1)>0);
    std::tuple<OptionValue,OptionValue> res
        = BlackScholes(d_expiry,d_strike,d_spot,d_rate,d_div,d_vol);
    if (d_isPut) { d_BlackScholes = std::get<1>(res);}
    else {d_BlackScholes = std::get<0>(res);}

    std::cout << "C(S,K) = " << d_BlackScholes.price << std::endl;

    std::cout << "S = " << d_spot << std::endl;
    std::cout << "K = " << d_strike << std::endl;

    d_terminalConditionOverriden = false;
}

void FiniteDifference::setToBarrierOption(double B)
{
    double a = 2*(d_rate-d_div)/(d_vol*d_vol)-1;
    double fac = std::pow(B/d_spot, a);

    OptionValue BlackScholes2;
    std::tuple<OptionValue,OptionValue> res2
        = BlackScholes(d_expiry,d_strike, B*B/d_spot, d_rate, d_div, d_vol);
    if (d_isPut) { BlackScholes2 = std::get<1>(res2);}
    else {BlackScholes2 = std::get<0>(res2);}

    d_BlackScholes.price -= fac*BlackScholes2.price;
    std::cout << "Exact barrier solution = " << d_BlackScholes.price << std::endl;
}

void FiniteDifference::setDefaultDomain() 
{
    // terminal time
    d_tf = 0.5*d_expiry*d_vol*d_vol;
    d_ti = 0.0; 

    // x-boundaries
    double xmid 
        = std::log(d_spot/d_strike) + (d_rate-d_div-0.5*d_vol*d_vol)*d_expiry;
    d_xl = xmid - 3*d_vol*std::sqrt(d_expiry);
    d_xr = xmid + 3*d_vol*std::sqrt(d_expiry);

    // set boundaries on boundary conditions
    d_f->setBoundaries(d_ti, d_tf, d_xl, d_xr);
    d_gl->setBoundaries(d_ti, d_tf, d_xl, d_xr);
    d_gr->setBoundaries(d_ti, d_tf, d_xl, d_xr);
}

void FiniteDifference::setExpandedDomain(int M, double alpha, double T, double qDiv)
{
    // terminal time
    d_tf = 0.5*d_expiry*d_vol*d_vol;
    d_ti = 0.0;

    // x-boundaries
    double x0 = std::log(d_spot/d_strike);
    double xCompute = x0 + std::log(1-qDiv);
    double xmid = x0 + (d_rate-d_div-0.5*d_vol*d_vol)*T;
    double xlTemp = xmid - 3*d_vol*std::sqrt(T);
    double xrTemp = xmid + 3*d_vol*std::sqrt(T);

    double dt = (d_tf-d_ti)/M;
    double dx = std::sqrt(dt/alpha);
    int Nl = static_cast<int>(std::ceil((xCompute - xlTemp)/dx));
    int Nr = static_cast<int>(std::ceil((xrTemp - xCompute)/dx));

    std::cout << "Nl = " << Nl << std::endl;

    d_xl = xCompute - Nl*dx;
    d_xr = xCompute + Nr*dx;

    // set boundaries on boundary conditions
    d_f->setBoundaries(d_ti, d_tf, d_xl, d_xr);
    d_gl->setBoundaries(d_ti, d_tf, d_xl, d_xr);
    d_gr->setBoundaries(d_ti, d_tf, d_xl, d_xr);
}

void FiniteDifference::setDomain(double xl, double xr, double tf, double ti)
{
    d_tf = tf;
    d_ti = ti;
    // d_tf-d_ti has to be equal to d_expiry...
    d_xl = xl;
    d_xr = xr;
    
    d_f->setBoundaries(d_ti, d_tf, d_xl, d_xr);
    d_gl->setBoundaries(d_ti, d_tf, d_xl, d_xr);
    d_gr->setBoundaries(d_ti, d_tf, d_xl, d_xr);
}

void FiniteDifference::discretizeDomainByTimeStepsAndAlphaTemp( int M, 
                                                    double alphaTemp )
{
    d_M = M;
    double dt = (d_tf-d_ti)/d_M;
    double dxTemp = std::sqrt(dt/alphaTemp);
    d_N = static_cast<int>(std::floor((d_xr - d_xl)/dxTemp));
    double dx = (d_xr - d_xl)/d_N;
    d_alpha = dt/(dx*dx);
}

void FiniteDifference::discretizeDomainByIntervalsAndAlphaTemp( int N,
                                                   double alphaTemp )
{
    d_N = N;
    double dx = (d_xr - d_xl)/d_N;
    double dtTemp = alphaTemp*dx*dx;
    d_M = static_cast<int>(std::ceil((d_tf-d_ti)/dtTemp));
    double dt = (d_tf-d_ti)/d_M;
    d_alpha = dt/(dx*dx);
}

void FiniteDifference::discretizeDomainByTimeStepsAndAlphaFixed(int M, double alpha)
{
    d_M = M;
    double dt = (d_tf-d_ti)/d_M;
    double dx = std::sqrt(dt/alpha);
    d_N = static_cast<int>((d_xr - d_xl)/dx);
    d_alpha = alpha;
}

void FiniteDifference::discretizeDomain(int M, int N)
{
    d_M = M;
    d_N = N;
    double dx = (d_xr-d_xl)/d_N;
    double dt = (d_tf-d_ti)/d_M;
    d_alpha = dt/(dx*dx);
}

void FiniteDifference::overrideTerminalCondition(const std::vector<double> & u)
{
    d_terminalConditionOverriden = true;
    d_u0 = u;
}

OptionValue FiniteDifference::BlackScholesValue() const
{
    return d_BlackScholes;
}

OptionValue FiniteDifference::getOptionValue( const std::vector<double> & u, 
                                              const std::vector<double> & uOld,
                                              double dt )
{
    // Use the numerical solution to find option values and greeks
    int N = u.size()-1;
    double c = (d_rate-d_div)/(d_vol*d_vol);
    double a = c-0.5;
    double b = (c+0.5)*(c+0.5) + 2*d_div/(d_vol*d_vol);
    double dx = (d_xr-d_xl)/N;
    
    double xCompute = std::log(d_spot/d_strike);
    int spotIndex = (xCompute-d_xl)/dx;

    std::cout << "spot index = " << spotIndex << std::endl;

    double xi = d_xl + spotIndex*dx;
    double xj = d_xl + (spotIndex+1)*dx;
    double si = d_strike*std::exp(xi);
    double sj = d_strike*std::exp(xj);
    double ui = u[spotIndex];
    double uj = u[spotIndex+1];
    double vi = std::exp(-a*xi-b*d_tf)*ui;
    double vj = std::exp(-a*xj-b*d_tf)*uj;

    // Approximation I: interpolation of v
    double vApproximate1 = ((sj-d_spot)*vi + (d_spot-si)*vj)/(sj-si);

    // Approximation II: interpolation of u
    double uApproximate = ((xj-xCompute)*ui + (xCompute-xi)*uj)/(xj-xi);
    double vApproximate2 = std::exp(-a*xCompute-b*d_tf)*uApproximate;

    // Compute Greeks I: Delta and Gamma
    double xii = d_xl + (spotIndex-1)*dx;
    double xjj = d_xl + (spotIndex+2)*dx;
    double sii = d_strike*std::exp(xii);
    double sjj = d_strike*std::exp(xjj);
    double uii = u[spotIndex-1];
    double ujj = u[spotIndex+2];
    double vii = std::exp(-a*xii-b*d_tf)*uii;
    double vjj = std::exp(-a*xjj-b*d_tf)*ujj;

    double delta = (vj-vi)/(sj-si);
    double deltaLeft = (vi-vii)/(si-sii);
    double deltaRight = (vjj-vj)/(sjj-sj);
    double gamma = 2*(deltaRight-deltaLeft)/(sjj+sj-si-sii);

    double deltaCentral = (vj-vii)/(sj-sii);
    double gammaCentral = ((si-sii)*vj-(sj-sii)*vi+(sj-si)*vii)/(0.5*(si-sii)*(sj-si)*(sj-sii));

    // Compute Greeks II: Theta
    double dT = 2*dt/(d_vol*d_vol);
    double uOldi = uOld[spotIndex];
    double uOldj = uOld[spotIndex+1];
    double vOldi = std::exp(-a*xi-b*(d_tf-dt))*uOldi;
    double vOldj = std::exp(-a*xj-b*(d_tf-dt))*uOldj;

    double vOldApproximate1 = ((sj-d_spot)*vOldi + (d_spot-si)*vOldj)/(sj-si);
    // minus sign: backward B-S equation converted to forward heat equation
    double theta = -(vApproximate1-vOldApproximate1)/dT;

    double thetaForward = -(vi-vOldi)/dT;
    std::cout << "Central," 
              << std::fixed
              << std::setprecision(9)
              << ui << ","
              << vi << ","
              << deltaCentral << ","
              << gammaCentral << ","
              << thetaForward
              << std::endl;

    OptionValue optionValue;
    optionValue.price = vApproximate1;
    optionValue.price2 = vApproximate2;
    optionValue.delta = delta;
    optionValue.gamma = gamma;
    optionValue.theta = theta;

    return optionValue;
}

std::vector<double> FiniteDifference::evaluate( FiniteDifferenceMethod fdm, 
                                                double omega, 
                                                bool varReduction, 
                                                double vExact )
{
    if (d_terminalConditionOverriden) {
        d_h = HeatPDE(d_xl, d_xr, d_tf, d_ti, d_u0, (*d_gl), (*d_gr), (*d_prem));
    }
    else {
        d_h = HeatPDE(d_xl, d_xr, d_tf, d_ti, (*d_f), (*d_gl), (*d_gr), (*d_prem));
    }

    // Discretize the computational domain
    double dt = (d_tf-d_ti)/d_M;
    
    // Discretize and solve the PDE
    std::vector<double> u(d_N+1,0);
    if (d_terminalConditionOverriden) { u = d_u0; }

    int dM=0;
    int dN=0;
    bool isAmerican = false;
    bool xReversed = false;
    switch (fdm)
    {
        case EulerForward: 
            isAmerican = false;
            xReversed = false;
            d_h.fdSolveForwardEuler(d_M, d_N, &u, dM, dN);
            break;
        case EulerBackwardByLU:
            isAmerican = false;
            xReversed = false;
            d_h.fdSolveBackwardEulerByLU(d_M, d_N, &u, dM, dN);
            break;
        case EulerBackwardBySOR:
            isAmerican = false;
            xReversed = false;
            d_h.fdSolveBackwardEulerBySOR(d_M, d_N, omega, &u, dM, dN);
            break;
        case CrankNicolsonByLU:
            isAmerican = false;
            xReversed = false;
            d_h.fdSolveCrankNicolsonByLU(d_M, d_N, &u, dM, dN);
            break;
        case CrankNicolsonBySOR:
            isAmerican = false;
            xReversed = false;
            d_h.fdSolveCrankNicolsonBySOR(d_M, d_N, omega, &u, dM, dN);
            break;
        case AmericanEulerForward: 
            isAmerican = true;
            xReversed = false;
            d_h.fdSolveAmericanForwardEuler(d_M, d_N, &u, dM, dN);
            break;
        case AmericanEulerBackwardByLU: 
            isAmerican = true;
            xReversed = true;
            d_h.fdSolveAmericanBackwardEulerByLU(d_M, d_N, &u, dM, dN);
            break;
        case AmericanEulerBackwardBySOR: 
            isAmerican = true;
            xReversed = false;
            d_h.fdSolveAmericanBackwardEulerBySOR(d_M, d_N, omega, &u, dM, dN);
            break;
        case AmericanCrankNicolsonByLU:
            isAmerican = true;
            xReversed = true;
            d_h.fdSolveAmericanCrankNicolsonByLU(d_M, d_N, &u, dM, dN);
            break;
        case AmericanCrankNicolsonBySOR:
            isAmerican = true;
            xReversed = false;
            d_h.fdSolveAmericanCrankNicolsonBySOR(d_M, d_N, omega, &u, dM, dN);
            break;
    }
    
    // get early exercise boundary
    if (isAmerican && dM>0) {
        std::cout << " --- Print Early Exercise Boundary --- " << std::endl;
        double dx = (d_xr-d_xl)/d_N;
        std::vector<double> eeb = d_h.getEarlyExerciseBoundary();
        for (int m=0; m<d_M+1; m++) { 
            if ( m%dM == 0 ) {
                double s1 = d_strike*std::exp(eeb[m]);
                double s2 = d_strike*std::exp(eeb[m]+dx);
                double sOpt = 0.5*(s1+s2);
                double t = d_expiry - 2*m*dt/(d_vol*d_vol);
                if (m==0) { sOpt = d_strike; }
                std::cout << std::fixed << std::setprecision(9)
                          << t << "," << sOpt << std::endl;
            }
        }
        std::cout << " ------------------------------------- " << std::endl;
    }
   
    // get solution at the previous time step 
    std::vector<double> uOld = d_h.getPrevSolution(xReversed);
    
    // Option price and greeks
    OptionValue optionValue = getOptionValue(u, uOld, dt);
    double vApproximate1 = optionValue.price;
    double vApproximate2 = optionValue.price2;
    double delta = optionValue.delta;
    double gamma = optionValue.gamma;
    double theta = optionValue.theta;
    
    // variance reduction for American options
    double vApproximate3 = 0.0;
    if ( isAmerican && varReduction ) {
        std::vector<double> uEuro(d_N+1,0);
        switch (fdm)
        {
            case EulerForward: 
            case EulerBackwardByLU:
            case EulerBackwardBySOR:
            case CrankNicolsonByLU:
            case CrankNicolsonBySOR:
                break;
            case AmericanEulerForward: 
                d_h.fdSolveForwardEuler(d_M, d_N, &uEuro, dM, dN);
                break;
            case AmericanEulerBackwardByLU: 
                d_h.fdSolveBackwardEulerByLU(d_M, d_N, &uEuro, dM, dN);
                break;
            case AmericanEulerBackwardBySOR: 
                d_h.fdSolveBackwardEulerBySOR(d_M, d_N, omega, &uEuro, dM, dN);
                break;
            case AmericanCrankNicolsonByLU:
                d_h.fdSolveCrankNicolsonByLU(d_M, d_N, &uEuro, dM, dN);
                break;
            case AmericanCrankNicolsonBySOR:
                d_h.fdSolveCrankNicolsonBySOR(d_M, d_N, omega, &uEuro, dM, dN);
                break;
        }
    
        std::vector<double> uEuroOld = d_h.getPrevSolution(xReversed);
    
        OptionValue EuroValue = getOptionValue( uEuro, uEuroOld, dt );
        double vEuroApproximate1 = EuroValue.price;
        double vEuroExact = d_BlackScholes.price;
        vApproximate3 = vApproximate1 + (vEuroExact-vEuroApproximate1);
    }
    
    // Error analysis I: Point-wise Error
    double exactPrice = vExact > 0 ? vExact : d_BlackScholes.price;
    double error1 = std::fabs(vApproximate1-exactPrice);
    double error2 = std::fabs(vApproximate2-exactPrice);
    
    // Print the errors and greeks
    std::cout << std::fixed
              << std::setprecision(9)
              << error1 << ",,"
              << error2 << ",,";

    // Error analysis II: Root-Mean-Squared (RMS) Error
    // Note: only for European options.
    double error3 = 0;
    if (!isAmerican) {
        double c = (d_rate-d_div)/(d_vol*d_vol);
        double a = c-0.5;
        double b = (c+0.5)*(c+0.5) + 2*d_div/(d_vol*d_vol);
        double dx = (d_xr-d_xl)/d_N;

        int count = 0;
        for (int n=0; n<d_N+1; n++) {
            double x = d_xl + n*dx;
            double s = d_strike*std::exp(x);
            // Black-Scholes value
            OptionValue vBlackScholes;
            std::tuple<OptionValue,OptionValue> res
                = BlackScholes(d_expiry,d_strike,s,d_rate,d_div,d_vol);
            if (d_isPut) { vBlackScholes = std::get<1>(res);}
            else { vBlackScholes = std::get<0>(res);}

            double bsPrice = vBlackScholes.price;
            if (bsPrice > 0.00001*s) {
                double vApproximate = std::exp(-a*x-b*d_tf)*u[n];
                double diff = (vApproximate - bsPrice)/bsPrice;
                double diff2 = diff*diff;
                error3 += diff2;
                count++;
            }
        }

        error3 /= count;
        error3 = std::sqrt(error3);
        std::cout << error3 << ",,";
    }
    
    std::cout << std::fixed
              << std::setprecision(9)
              << delta << ","
              << gamma << ","
              << theta << ",";
    
    if ( isAmerican && varReduction ) {
        double error4 = std::fabs(vApproximate3-exactPrice);
        std::cout << vApproximate3 << "," << error4 << ",";
    }

    std::cout << std::endl;

    return u;
}

