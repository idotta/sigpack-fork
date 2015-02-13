// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.
#ifndef SP_RESAMPLING_H
#define SP_RESAMPLING_H
namespace sp
{
    //////////////////////////////////////////////////////////////////
    // Resampling functions
    //////////////////////////////////////////////////////////////////

    ///////////////////////////////////
    // xN = upsample(x,p)
    //      Upsamples x p times
    template <class T1>
    arma::Col<T1> upsample(arma::Col<T1> &x, const int p )
    {
        long int N = x.size();
        arma::Col<T1> y;
        y.set_size(p*N);
        y.zeros();
        for(long int n=0;n<N;n++)
            y[p*n] = x[n];
        return y;
    }

    ///////////////////////////////////
    // xN = upfir(x,p,K)
    //      Upsamples x p times using an K tap antialias filter
    template <class T1>
    arma::Col<T1> upfir(arma::Col<T1> &x, const int p, const int K )
    {
        FIR_filt<T1,double,T1> AA;
        AA.set_coeffs(fir1(K,1/float(p)));
        return AA.filter(upsample(x,p));
    }

    ///////////////////////////////////
    // xN = downsample(x,q)
    //      Downsamples x q times
    template <class T1>
    arma::Col<T1> downsample(arma::Col<T1> &x, const int q )
    {
        arma::Col<T1> y;
        int N = int(floor(1.0*x.size()/q));
        y.set_size(N);
        for(long int n=0;n<N;n++)
            y[n] = x[n*q];
        return y;
    }

    ///////////////////////////////////
    // xN = downfir(x,q,K)
    //      Downsamples x q times using an K tap antialias filter
    template <class T1>
    arma::Col<T1> downfir(arma::Col<T1> &x, const int q, const int K )
    {
        FIR_filt<T1,double,T1> AA;
        AA.set_coeffs(fir1(K,1/float(q)));
        return downsample(AA.filter(x),q);
    }
} // end namespace
#endif