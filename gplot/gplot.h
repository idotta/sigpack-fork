// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.
#ifndef SP_GPLOT_H
#define SP_GPLOT_H
namespace sp
{
    //////////////////////////////////////////////////////////////////
    // gplot Class
    //
    //    Verified with Gnuplot 4.6.5 for Win64 and Linux
    //    Inspiration from https://code.google.com/p/gnuplot-cpp/
    //
    //////////////////////////////////////////////////////////////////
    class gplot
    {
    private:
        FILE           *gnucmd;
        std::string    linestyle;

    public:
        gplot();
        ~gplot();
        void send2gp(const std::string &cmdstr);
        void figure(const int fig);
        void window(const int fig, const std::string &name,const int x=20,const int y=20,const int width=500,const int height=400);
        void window(const std::string &name,const int x=20,const int y=20,const int width=500,const int height=400);
        void set_linestyle(const std::string& style);
        void title(const std::string& label);
        void xlabel(const std::string& label);
        void ylabel(const std::string& label);
        void xlim(const double xmin, const double xmax);
        void ylim(const double ymin, const double ymax);
        void plot_str2(arma::vec &x, arma::vec &y);
        void plot(arma::vec &x, arma::vec &y, const std::string& label);
        void plot(arma::vec &y, const std::string& label);
        void plot(arma::vec &x, arma::vec &y1, arma::vec &y2, const std::string& label1, const std::string& label2);
        void plot(arma::vec &y1, arma::vec &y2, const std::string& label1, const std::string& label2);
        void scatter( arma::vec &x, arma::vec &y, const std::string& label);
        void scatter( arma::vec &x1, arma::vec &y1, const std::string& label1,
            arma::vec &x2, arma::vec &y2, const std::string& label2);
        void image( arma::mat &x);
        void mesh( arma::mat &x);
        void surf( arma::mat &x);
    }; // End Gnuplot Class

    ///////////////////////////////////
    // Constructor
    gplot::gplot()
    {
#if defined(WIN32)
        gnucmd = _popen("gnuplot -persist 2> NUL","w");
#define GP_TERM "win"
#elif defined(unix)
//        gnucmd = popen("gnuplot -persist &> /dev/null","w");
        gnucmd = popen("gnuplot -persist","w");
#define GP_TERM "x11"
        //#elif defined(_APPLE_)
        //            gnucmd = popen("gnuplot -persist &> /dev/null","w");
        //#define GP_TERM "aqua"
#else
#error Only Windows and Linux/Unix is supported so far!
#endif
        if(!gnucmd)
        {
            std::cout << "Could not start gnuplot" << std::endl;
        }
        setvbuf(gnucmd, NULL, _IOLBF, 512);

        // Set global params
        linestyle = "lines";
    }

    ///////////////////////////////////
    // Destructor
    gplot::~gplot()
    {
#if defined(WIN32)
        _pclose(gnucmd);
#elif defined(unix)
        pclose(gnucmd);
#endif
    }

    ///////////////////////////////////
    // Send command to Gnuplot pipe
    void gplot::send2gp(const std::string &cmdstr)
    {
        std::string tmp=cmdstr+"\n";
        std::fputs(tmp.c_str(), gnucmd );
        //std::cout << tmp.c_str() << std::endl;
    }

    ///////////////////////////////////
    // Set active figure
    void gplot::figure(const int fig)
    {
        std::ostringstream tmp_s;
        tmp_s << "set term " << GP_TERM << " " << fig;
        send2gp(tmp_s.str());
        send2gp("reset");
    }

    ///////////////////////////////////
    // Configurate the window
    void gplot::window(const int fig, const std::string &name,const int x,const int y,const int width,const int height)
    {
        std::ostringstream tmp_s;
        tmp_s << "set term " << GP_TERM << " " << fig << " title \"" << name << "\" position " << x << "," << y << " size " << width << "," << height;
        send2gp(tmp_s.str());
        send2gp("reset");
    }
    void gplot::window(const std::string &title,const int x,const int y,const int width,const int height)
    {
        window(0,title,x,y,width,height);
    }

    ///////////////////////////////////
    // Set line style
    //      lines,points,linespoints,dots,steps
    void gplot::set_linestyle(const std::string& style)
    {
        linestyle = style;
    }

    ///////////////////////////////////
    // Set x-axis label
    void gplot::xlabel(const std::string& label)
    {
        std::ostringstream tmp_s;
        tmp_s << "set xlabel \"" << label << "\" ";
        send2gp(tmp_s.str());
    }

    ///////////////////////////////////
    // Set y-axis label
    void gplot::ylabel(const std::string& label)
    {
        std::ostringstream tmp_s;
        tmp_s << "set ylabel \"" << label << "\" ";
        send2gp(tmp_s.str());
    }

    ///////////////////////////////////
    // Set title
    void gplot::title(const std::string& name)
    {
        std::ostringstream tmp_s;
        tmp_s << "set title \"" << name << " \" ";
        send2gp(tmp_s.str());
    }

    ///////////////////////////////////
    // Set x range
    void gplot::xlim(const double xmin, const double xmax)
    {
        std::ostringstream tmp_s;
        tmp_s << "set xrange [" << xmin << ":" << xmax << "]";
        send2gp(tmp_s.str());
    }

    ///////////////////////////////////
    // Set y range
    void gplot::ylim(const double ymin, const double ymax)
    {
        std::ostringstream tmp_s;
        tmp_s << "set yrange [" << ymin << ":" << ymax << "]";
        send2gp(tmp_s.str());
    }

    ///////////////////////////////////
    // Plot stream x and y
    void gplot::plot_str2( arma::vec &x, arma::vec &y)
    {
        std::ostringstream tmp_s;
        for(unsigned int n=0;n<y.size();n++)
        {
            tmp_s << x(n) << " " << y(n);
            send2gp(tmp_s.str());
            tmp_s.str(""); // Clear buffer
        }
        send2gp("e");
    }

    ///////////////////////////////////
    // Plot
    void gplot::plot( arma::vec &x, arma::vec &y, const std::string& label)
    {
        std::ostringstream tmp_s;
        send2gp("set key noautotitle");
        send2gp("set grid");
        tmp_s << "plot '-' title \"" << label << "\" with " << linestyle;
        send2gp(tmp_s.str());
        plot_str2(x,y);
    }

    void gplot::plot(arma::vec &y, const std::string& label)
    {
        arma::vec t(y.size());
        t = arma::linspace(1,y.size(),y.size());
        plot(t,y,label);
    }

    void gplot::plot( arma::vec &x, arma::vec &y1, arma::vec &y2, const std::string& label1, const std::string& label2)
    {
        std::ostringstream tmp_s;
        send2gp("set key noautotitle");
        send2gp("set grid");
        tmp_s << "plot '-' title \"" << label1 << "\" with " << linestyle << ", '-' title \"" << label2 << "\" with " << linestyle;
        send2gp(tmp_s.str());
        plot_str2(x,y1);
        plot_str2(x,y2);
    }

    void gplot::plot(arma::vec &y1, arma::vec &y2, const std::string& label1, const std::string& label2)
    {
        arma::vec t(y1.size());
        t = arma::linspace(1,y1.size(),y1.size());
        plot(t,y1,label1);
        plot(t,y2,label2);
    }

    ///////////////////////////////////
    // Scatter plot
    void gplot::scatter( arma::vec &x, arma::vec &y, const std::string& label)
    {
        set_linestyle("points");
        plot(x,y,label);
    }
    void gplot::scatter( arma::vec &x1, arma::vec &y1, const std::string& label1,
        arma::vec &x2, arma::vec &y2, const std::string& label2)
    {
        set_linestyle("points");
        std::ostringstream tmp_s;
        send2gp("set key noautotitle");
        send2gp("set grid");
        tmp_s << "plot '-' title \"" << label1 << "\" with " << linestyle << ", '-' title \"" << label2 << "\" with " << linestyle;
        send2gp(tmp_s.str());
        plot_str2(x1,y1);
        plot_str2(x2,y2);
    }

    ///////////////////////////////////
    // Image
    void gplot::image( arma::mat &x)
    {
        std::ostringstream tmp_s;
        xlim(-0.5,x.n_cols-0.5);
        ylim(x.n_rows-0.5,-0.5);
        tmp_s.str(""); // Clear buffer
        tmp_s << "plot '-' matrix with image";
        send2gp(tmp_s.str());
        for(unsigned int r=0;r<x.n_rows;r++)
        {
            tmp_s.str("");  // Clear buffer
            for(unsigned int c=0;c<x.n_cols;c++)
            {
                tmp_s << x(r,c) << " " ;
            }
            send2gp(tmp_s.str());
        }
        send2gp("e");
        send2gp("e");
    }

    ///////////////////////////////////
    // Mesh
    void gplot::mesh( arma::mat &x)
    {
        std::ostringstream tmp_s;
        send2gp("unset key");
        send2gp("set hidden3d");
        tmp_s << "splot '-' with lines";
        send2gp(tmp_s.str());

        for(unsigned int r=0;r<x.n_rows;r++)
        {
            for(unsigned int c=0;c<x.n_cols;c++)
            {
                tmp_s.str("");  // Clear buffer
                tmp_s << r << " " << c << " "<< x(r,c);
                send2gp(tmp_s.str());
            }
            send2gp("");
        }
        send2gp("e");
    }

    ///////////////////////////////////
    // Surf
    void gplot::surf( arma::mat &x)
    {
        send2gp("set pm3d");
        mesh(x);
    }
} // end namespace
#endif
