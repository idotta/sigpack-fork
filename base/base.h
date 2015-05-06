// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.
#ifndef SP_BASE_H
#define SP_BASE_H
namespace sp
{
    //////////////////////////////////////////////////////////////////
    // Math functions
    //////////////////////////////////////////////////////////////////
    const double PI   = 3.14159265358979323846;     // ... or use arma::datum::pi
    const double PI_2 = 6.28318530717958647692;

    ///////////////////////////////////
    // y = sinc(x)
    //      Sinc function
    double sinc( double x )
    {
        if(x==0.0)
            return 1.0;
        else
            return std::sin(PI*x)/(PI*x);
    }

    ///////////////////////////////////
    // y = besseli0()
    //      Modified first kind bessel function order zero
    double besseli0( double x )
    {
        double y=1.0,s=1.0,x2=x*x;
        int n = 1;
        while (s > y*1.0e-9)
        {
            s *= x2/4.0/(n*n);
            y += s;
            n++;
        }
        return y;
    }

    ///////////////////////////////////
    // p = angle(x)
    //      Calculates angle in radians for complex input
    template <typename T>
    double angle( std::complex<T> &x )
    {
        return std::arg(x);
    }
    arma::mat angle( arma::cx_mat &x )
    {
        arma::mat P;
        P.copy_size(x);
        for(unsigned int r=0;r<x.n_rows;r++)
            for(unsigned int c=0;c<x.n_cols;c++)
                P(r,c) = std::arg(x(r,c));
        return P;
    }
    arma::vec angle( arma::cx_vec &x )
    {
        arma::vec P;
        P.copy_size(x);
        for(unsigned int r=0;r<x.n_rows;r++)
            P(r) = std::arg(x(r));
        return P;
    }

    ///////////////////////////////////
    // err_handler("Error string")
    //      Prints an error message, waits for input and
    //      then exits with error
#define err_handler(msg) \
    { \
        std::cout << "SigPack Error [" << __FILE__  << "@" << __LINE__ << "]: " << msg << std::endl; \
        std::cin.get(); \
        exit(EXIT_FAILURE); \
    }

    ///////////////////////////////////
    // wrn_handler("Warning string")
    //      Prints an warning message
#define wrn_handler(msg)  std::cout << "SigPack warning [" << __FILE__ << "@" << __LINE__ << "]: " << msg << std::endl;

} // end namespace
#endif

