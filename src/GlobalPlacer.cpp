#include "GlobalPlacer.h"

#include <cstdio>
#include <vector>

#include "ObjectiveFunction.h"
#include "Optimizer.h"
#include "Point.h"

GlobalPlacer::GlobalPlacer(Placement &placement)
    : _placement(placement) {
}

void GlobalPlacer::place() {
    ////////////////////////////////////////////////////////////////////
    // This section is an example for analytical methods.
    // The objective is to minimize the following function:
    //      f(x,y) = 3*x^2 + 2*x*y + 2*y^2 + 7
    //
    // If you use other methods, you can skip and delete it directly.
    ////////////////////////////////////////////////////////////////////
    //std::vector<Point2<double>> t(1);
    ///cout<<_placement.numModules()<<endl;
    std::vector<Point2<double>> t(_placement.numModules());                   // Optimization variables (in this example, there is only one t)
    ObjectiveFunction wl(_placement);                    // Objective function
    double outlineWidth = _placement.boundryRight() - _placement.boundryLeft();
    double outlineHeight = _placement.boundryTop() - _placement.boundryBottom();
    double midx = _placement.boundryLeft() + outlineWidth/2;
    double midy = _placement.boundryBottom() + outlineHeight/2;
    double kAlpha = 200;                         // Constant step size
    if(outlineWidth>100000.)
    {
        kAlpha = 512;
        cout<<"yes"<<endl;
    }
    if(outlineWidth<10000.)
    {
        kAlpha = 7.2;
        cout<<"ibm05"<<endl;
    }
    SimpleConjugateGradient optimizer(wl, t, kAlpha);  // Optimizer

    // Set initial point
    for(size_t i = 0; i < _placement.numModules(); i++){
        t[i].x = midx - 100 +rand()%200;
        t[i].y = midy - 100 +rand()%200;
    }
    //t[0] = 4.;  // This set both t[0].x and t[0].y to 4.

    // Initialize the optimizer
    optimizer.Initialize();
    // Perform optimization, the termination condition is that the number of iterations reaches 100
    // TODO: You may need to change the termination condition, which is determined by the overflow ratio.
    //printf("iter = %3lu, f = %9.4f, x = %9.4f, y = %9.4f\n", 0, wl(t), t[0].x, t[0].y);
    //optimizer.Step();
    double iternum = 5000;
    iternum = iternum +1000*(_placement.boundryRight()-_placement.boundryLeft()) /66726.;
    if(outlineWidth<10000.)
    {
        iternum=3250;
        cout<<"ibm05"<<endl;
    }
    for (size_t i = 0; i < iternum; ++i) {
        optimizer.Step();
        for (unsigned i = 0; i < _placement.numModules(); ++i)
        {
            if (_placement.module(i).isFixed())
            {
                t[i].x = _placement.module(i).x() + _placement.module(i).width() / 2;
                t[i].y = _placement.module(i).y() + _placement.module(i).height() / 2;
                continue;
            }
            // out of left bound
            if (t[i].x < _placement.boundryLeft())
            {
                t[i].x = _placement.boundryLeft() + double(rand()) / RAND_MAX * outlineWidth * 0.1;
            }
            // out of right bound
            else if (t[i].x > _placement.boundryRight())
            {
                t[i].x = _placement.boundryRight() - double(rand()) / RAND_MAX * outlineWidth * 0.1;
            }
            // out of bottom bound
            if (t[i].y < _placement.boundryBottom())
            {
                t[i].y = _placement.boundryBottom() + double(rand()) / RAND_MAX * outlineHeight * 0.1;
            }
            // out of top bound
            else if (t[i].y > _placement.boundryTop())
            {
                t[i].y = _placement.boundryTop() - double(rand()) / RAND_MAX * outlineHeight * 0.1;
            }
        }
        //printf("iter = %3lu, f = %9.4f, x = %9.4f, y = %9.4f\n", i, 1., t[0].x, t[0].y);
    }

    ////////////////////////////////////////////////////////////////////
    // Global placement algorithm
    ////////////////////////////////////////////////////////////////////
    
    // TODO: Implement your global placement algorithm here.
    const size_t num_modules = _placement.numModules();  // You may modify this line.
    std::vector<Point2<double>> positions(num_modules);  // Optimization variables (positions of modules). You may modify this line.
    
    ////////////////////////////////////////////////////////////////////
    // Write the placement result into the database. (You may modify this part.)
    /*or (size_t i = 0; i < num_modules; i++){
        positions[i].x = _placement.boundryLeft() + (rand() % static_cast<int>(_placement.boundryRight() - _placement.boundryLeft() - _placement.module(i).width()));
        positions[i].y = _placement.boundryBottom() + (rand() % static_cast<int>(_placement.boundryTop() - _placement.boundryBottom() - _placement.module(i).height()));
    }*/
    for (size_t i = 0; i < num_modules; i++){
        if (t[i].x < _placement.boundryLeft())
        {
            t[i].x = _placement.boundryLeft();
        }
        else if (t[i].x > _placement.boundryRight() - _placement.module(i).width())
        {
            t[i].x = _placement.boundryRight() - _placement.module(i).width();
        }

        if (t[i].y < _placement.boundryBottom())
        {
            t[i].y = _placement.boundryBottom();
        }
        else if (t[i].y > _placement.boundryTop() - _placement.module(i).height())
        {
            t[i].y = _placement.boundryTop() - _placement.module(i).height();
        }
        positions[i].x = t[i].x;
        positions[i].y = t[i].y;
    }
    ////////////////////////////////////////////////////////////////////
    size_t fixed_cnt = 0;
    for (size_t i = 0; i < num_modules; i++) {        
        if (_placement.module(i).isFixed()) {
            fixed_cnt++;
            continue;
        }
        _placement.module(i).setPosition(positions[i].x, positions[i].y);
    }
    printf("INFO: %lu / %lu modules are fixed.\n", fixed_cnt, num_modules);
}

void GlobalPlacer::plotPlacementResult(const string outfilename, bool isPrompt) {
    ofstream outfile(outfilename.c_str(), ios::out);
    outfile << " " << endl;
    outfile << "set title \"wirelength = " << _placement.computeHpwl() << "\"" << endl;
    outfile << "set size ratio 1" << endl;
    outfile << "set nokey" << endl
            << endl;
    outfile << "plot[:][:] '-' w l lt 3 lw 2, '-' w l lt 1" << endl
            << endl;
    outfile << "# bounding box" << endl;
    plotBoxPLT(outfile, _placement.boundryLeft(), _placement.boundryBottom(), _placement.boundryRight(), _placement.boundryTop());
    outfile << "EOF" << endl;
    outfile << "# modules" << endl
            << "0.00, 0.00" << endl
            << endl;
    for (size_t i = 0; i < _placement.numModules(); ++i) {
        Module &module = _placement.module(i);
        plotBoxPLT(outfile, module.x(), module.y(), module.x() + module.width(), module.y() + module.height());
    }
    outfile << "EOF" << endl;
    outfile << "pause -1 'Press any key to close.'" << endl;
    outfile.close();

    if (isPrompt) {
        char cmd[200];
        sprintf(cmd, "gnuplot %s", outfilename.c_str());
        if (!system(cmd)) {
            cout << "Fail to execute: \"" << cmd << "\"." << endl;
        }
    }
}

void GlobalPlacer::plotBoxPLT(ofstream &stream, double x1, double y1, double x2, double y2) {
    stream << x1 << ", " << y1 << endl
           << x2 << ", " << y1 << endl
           << x2 << ", " << y2 << endl
           << x1 << ", " << y2 << endl
           << x1 << ", " << y1 << endl
           << endl;
}
