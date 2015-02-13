// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.
#ifndef SP_WINDOW_H
#define SP_WINDOW_H
namespace sp
{   
    //////////////////////////////////////////////////////////////////
    // Windows
    //////////////////////////////////////////////////////////////////

    ///////////////////////////////////
    // b = cos_win(N,a)
    //      Generic fifth order symmetric cos window
    arma::vec cos_win( int N, arma::vec &a )
    {
        arma::vec h(N);
        for(int i=0;i<N;i++)
        {
            h[i] = a[0] - a[1]*std::cos(1.0*PI_2*i/(N-1)) + a[2]*std::cos(2.0*PI_2*i/(N-1)) \
                        - a[3]*std::cos(3.0*PI_2*i/(N-1)) + a[4]*std::cos(4.0*PI_2*i/(N-1));
        }
        return h;
    }

    ///////////////////////////////////
    // b = hamming(N)
    arma::vec hamming( int N )
    {
        arma::vec a=arma::zeros<arma::vec>(5);
        a[0] = 0.54;
        a[1] = 0.46;
        return cos_win(N,a);
    }

    ///////////////////////////////////
    // b = hann(N)
    arma::vec hann( int N )
    {
        arma::vec a=arma::zeros<arma::vec>(5);
        a[0] = 0.5;
        a[1] = 0.5;
        return cos_win(N,a);
    }

    ///////////////////////////////////
    // b = blackman(N)
    arma::vec blackman( int N )
    {
        arma::vec a=arma::zeros<arma::vec>(5);
        a[0] = 0.42; // 7938/18608.0
        a[1] = 0.5;  // 9240/18608.0
        a[2] = 0.08; // 1430/18608.0
        return cos_win(N,a);
    }

    ///////////////////////////////////
    // b = blackmanharris(N)
    //      Symmetric BH4 window
    arma::vec blackmanharris( int N )
    {
        arma::vec a=arma::zeros<arma::vec>(5);
        a[0] = 0.35875;
        a[1] = 0.48829;
        a[2] = 0.14128;
        a[3] = 0.01168;
        return cos_win(N,a);
    }

    ///////////////////////////////////
    // b = flattopwin(N)
    arma::vec flattopwin( int N )
    {
        arma::vec a=arma::zeros<arma::vec>(5);
        a[0] = 0.21557895;
        a[1] = 0.41663158;
        a[2] = 0.277263158;
        a[3] = 0.083578947;
        a[4] = 0.006947368;
        return cos_win(N,a);
    }

    ///////////////////////////////////
    // b = hanning(N)
    arma::vec hanning( int N )
    {
        arma::vec h(N);
        for(int i=0;i<N;i++)
        {
            h[i] = 0.5-0.5*std::cos(PI_2*(i+1)/(N+1));
        }
        return h;
    }

    ///////////////////////////////////
    // b = kaiser(N,beta)
    arma::vec kaiser( int N, double beta )
    {
        arma::vec h(N);
        double bb = besseli0(beta);
        for(int i=0;i<N;i++)
        {
            h[i] = besseli0(beta*sqrt(4.0*i*(N-1-i))/(N-1))/bb;
        }
        return h;
    }

    ///////////////////////////////////
    // b = triang(N)
    arma::vec triang( int N )
    {
        arma::vec h(N);
        if(N%2)    // Odd
        {
            for(int i=0;i<(N-1)/2;i++)
            {
                h[i]     = 2.0*(i+1)/(N+1);
                h[N-i-1] = h[i];
            }
            h[(N-1)/2] = 1.0;
        }
        else      // Even
        {
            for(int i=0;i<N/2;i++)
            {
                h[i]     = (2.0*i+1)/N;
                h[N-i-1] = h[i];
            }
        }
        return h;
    }
} // end namespace
#endif
