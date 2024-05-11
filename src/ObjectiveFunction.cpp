#include "ObjectiveFunction.h"

#include "cstdio"

ExampleFunction::ExampleFunction(Placement &placement) : BaseFunction(placement.numModules()), placement_(placement)
{
    printf("Fetch the information you need from placement database.\n");
    printf("For example:\n");
    printf("    Placement boundary: (%.f,%.f)-(%.f,%.f)\n", placement_.boundryLeft(), placement_.boundryBottom(),
           placement_.boundryRight(), placement_.boundryTop());
}

const double &ExampleFunction::operator()(const std::vector<Point2<double>> &input)
{
    // Compute the value of the function
    value_ = 3. * input[0].x * input[0].x + 2. * input[0].x * input[0].y +
             2. * input[0].y * input[0].y + 7.;
    cout<<input.size()<<endl;
    input_ = input;
    return value_;
}

const std::vector<Point2<double>> &ExampleFunction::Backward()
{
    // Compute the gradient of the function
    grad_[0].x = 6. * input_[0].x + 2. * input_[0].y;
    grad_[0].y = 2. * input_[0].x + 4. * input_[0].y;
    
    return grad_;
}

Wirelength::Wirelength(Placement &placement) : BaseFunction(placement.numModules()), placement_(placement)
{
    printf("Fetch the information you need from placement database.\n");
    printf("For example:\n");
    printf("    Placement boundary: (%.f,%.f)-(%.f,%.f)\n", placement_.boundryLeft(), placement_.boundryBottom(),
           placement_.boundryRight(), placement_.boundryTop());
}

const double &Wirelength::operator()(const std::vector<Point2<double>> &input)
{
    // Compute the value of the function
    double gamma = 500;
    value_ = 0;
    /*value_ = 3. * input[0].x * input[0].x + 2. * input[0].x * input[0].y +
             2. * input[0].y * input[0].y + 7.;*/
    for(unsigned i=0;i<placement_.numNets();i++){
        Net net = placement_.net(i);
        double xmax = 0;
        double ymax = 0;
        double xmin = 0;
        double ymin = 0;
        for(unsigned j=0;j<net.numPins();j++){
            Pin pin = net.pin(j);
            int id = pin.moduleId();
            xmax += exp(input[id].x/gamma);
            xmin += exp(-input[id].x/gamma);
            ymax += exp(input[id].y/gamma);
            ymin += exp(-input[id].y/gamma);
            
        }
        //cout<<exp(-input[0].x/gamma)<<endl;
        //if(i==0)cout<<"f:xmax "<<xmax<<" xmin "<<xmin<<" ymax "<<ymax<<" ymin "<<ymin<<endl;
        xmax = log(xmax);
        xmin = log(xmin);
        ymax = log(ymax);  
        ymin = log(ymin);

        value_ = value_ + xmax + xmin + ymax + ymin;
    }
    value_ = value_ * gamma;
    input_ = input;
    return value_;
}

const std::vector<Point2<double>> &Wirelength::Backward()
{
    for(unsigned i=0;i<placement_.numModules();i++){
        grad_[i]=0;
    }
    double gamma = 500;
    // Compute the gradient of the function

    //grad_[i].x=0; grad_[i].y=0;
    for(unsigned i=0;i<placement_.numNets();i++){
        Net net = placement_.net(i);
        double xmax = 0;
        double ymax = 0;
        double xmin = 0;
        double ymin = 0;
        for(unsigned j=0;j<net.numPins();j++){
            Pin pin = net.pin(j);
            int id = pin.moduleId();
            xmax += exp(input_[id].x/gamma);
            xmin += exp(-input_[id].x/gamma);
            ymax += exp(input_[id].y/gamma);
            ymin += exp(-input_[id].y/gamma);
            //Module module = _modules[pin.moduleId()];
            //Point2<double> pos = module.position();
            
        }
        //if(i==0)cout<<"xmax "<<xmax<<" xmin "<<xmin<<" ymax "<<ymax<<" ymin "<<ymin<<endl;
        for(unsigned j=0;j<net.numPins();j++){
            Pin pin = net.pin(j);
            int id = pin.moduleId();
            grad_[id].x +=  ((exp(input_[id].x/gamma) ) / xmax ) - ((exp(-input_[id].x/gamma) ) / xmin );
            grad_[id].y +=  ((exp(input_[id].y/gamma) ) / ymax ) - ((exp(-input_[id].y/gamma) ) / ymin );
            
        }

    }
    //cout<<"0_grad "<<grad_[1].x<<" "<<grad_[1].y<<endl;
    return grad_;
}


Density::Density(Placement &placement) : BaseFunction(placement.numModules()), placement_(placement)
{
    cout<<"do "<<endl;
    binSize_=480;
    double outlineWidth = placement_.boundryRight() - placement_.boundryLeft();
    //double outlineHeight = placement_.boundryTop() - placement_.boundryBottom();
    if(outlineWidth<10000){
        cout<<"ibm05binsize"<<endl;
        binSize_=14;
    }
    if(outlineWidth>100000){
        cout<<"ibm07binsize"<<endl;
        binSize_=960;
    }
    Mb = 0.85;
    widthBinNum_ = (placement_.boundryRight()-placement_.boundryLeft())/binSize_;
    heightBinNum_ = (placement_.boundryTop()-placement_.boundryBottom())/binSize_;
    binDensity.resize(widthBinNum_);
    for(int i=0;i<widthBinNum_;i++){
        binDensity[i].resize(heightBinNum_);
    }
}

const double &Density::operator()(const std::vector<Point2<double>> &input)
{
    // Compute the value of the function
    
    value_=0;
    overflowRatio_ = 0;
    for(int i=0;i<widthBinNum_;i++){
        for(int j=0;j<heightBinNum_;j++){
            binDensity[i][j]=0;
        }
    }
    int dx ,dy =0;
    int binleft,binright,bintop,binbottom =0;
    int nearbinleft,nearbinright,nearbintop,nearbinbottom =0;
    for(unsigned k=0;k<placement_.numModules();k++){
        binleft = ((input[k].x-placement_.boundryLeft())/binSize_);
        binright = ((input[k].x+placement_.module(k).width()-placement_.boundryLeft())/binSize_);
        bintop = ((input[k].y+placement_.module(k).height()-placement_.boundryBottom())/binSize_);
        binbottom = ((input[k].y-placement_.boundryBottom())/binSize_);
        if((binleft-3)>0)nearbinleft=binleft-3;
        else nearbinleft=0;
        if((binright+3)<widthBinNum_)nearbinright=binright+3;
        else nearbinright=widthBinNum_;
        if((bintop+3)<heightBinNum_)nearbintop=bintop+3;
        else nearbintop=heightBinNum_;
        if((binbottom-3)>0)nearbinbottom=binbottom-3;
        else nearbinbottom=0;
        double mWidth = placement_.module(k).width(), mHeight = placement_.module(k).height();
        double alphax = 4/((mWidth+2*binSize_)*(mWidth+4*binSize_));
        double betax = 2/(binSize_*(mWidth+4*binSize_));
        double alphay = 4/((mHeight+2*binSize_)*(mHeight+4*binSize_));
        double betay = 2/(binSize_*(mHeight+4*binSize_));
        for(int i=nearbinleft;i<nearbinright;i++){
            for(int j=nearbinbottom;j<nearbintop;j++){
                
                dx=abs((input[k].x+mWidth/2)-(placement_.boundryLeft()+i*binSize_+0.5 * binSize_));
                dy=abs((input[k].y+mHeight/2)-(placement_.boundryBottom()+j*binSize_+0.5 * binSize_));
                double xl,yl=0;
                if( (dx>=0) && dx<(mWidth/2+binSize_)){
                    xl=1-alphax*dx*dx;
                }
                else if( (dx>(mWidth/2+binSize_)) && dx<(mWidth/2+2*binSize_)){
                    xl=betax*(dx-mWidth/2-2*binSize_)*(dx-mWidth/2-2*binSize_);
                }
                else{
                    xl=0;
                }
                if(dy>=0 && dy<(mHeight/2+binSize_)){
                    yl=1-alphay*dy*dy;
                }
                else if(dy>(mHeight/2+binSize_) && dy<(mHeight/2+2*binSize_)){
                    yl=betay*(dy-mHeight/2-2*binSize_)*(dy-mHeight/2-2*binSize_);
                }
                else{
                    yl=0;
                }
                binDensity[i][j] += xl*yl;
            }
        }
    }


    /*for(int i=0;i<widthBinNum_;i++){
        for(int j=0;j<heightBinNum_;j++){
            for(unsigned k=0;k<placement_.numModules();k++){
                
                double mWidth = placement_.module(k).width(), mHeight = placement_.module(k).height();
                double alpha = 4/((mWidth+2*binSize_)*(mWidth+4*binSize_));
                double beta = 2/(binSize_*(mWidth+4*binSize_));
                dx=abs((input[k].x+mWidth/2)-(placement_.boundryLeft()+i*binSize_+0.5 * binSize_));
                dy=abs((input[k].y+mHeight/2)-(placement_.boundryBottom()+j*binSize_+0.5 * binSize_));
                double xl,yl=0;
                if( (dx>=0) && dx<(mWidth/2+binSize_)){
                    xl=1-alpha*dx*dx;
                }
                else if( (dx>(mWidth/2+binSize_)) && dx<(mWidth/2+2*binSize_)){
                    xl=beta*(dx-mWidth/2-2*binSize_)*(dx-mWidth/2-2*binSize_);
                }
                else{
                    xl=0;
                }
                if(dy>=0 && dy<(mHeight/2+binSize_)){
                    yl=1-alpha*dy*dy;
                }
                else if(dy>(mHeight/2+binSize_) && dy<(mHeight/2+2*binSize_)){
                    yl=beta*(dy-mHeight/2-2*binSize_)*(dy-mHeight/2-2*binSize_);
                }
                else{
                    yl=0;
                }
                binDensity[i][j] += xl*yl;
            }
        }
    }*/
    for(int i=0;i<widthBinNum_;i++){
        for(int j=0;j<heightBinNum_;j++){
            value_ += (binDensity[i][j]-Mb)*(binDensity[i][j]-Mb);
            overflowRatio_ += max(0.0, binDensity[i][j] - Mb);
        }
    }
    overflowRatio_ /= widthBinNum_ * heightBinNum_;
    cout<<"overflowRatio: "<<overflowRatio_<<endl;
    input_ = input;
    return value_;
}

const std::vector<Point2<double>> &Density::Backward()
{

    
    for(unsigned i=0;i<placement_.numModules();i++){
        grad_[i].x = 0;
        grad_[i].y = 0;
    }
    int dx ,dy =0;
    int signdx ,signdy =0;
    int binleft,binright,bintop,binbottom =0;
    int nearbinleft,nearbinright,nearbintop,nearbinbottom =0;
    for(unsigned k=0;k<placement_.numModules();k++){
        binleft = ((input_[k].x-placement_.boundryLeft())/binSize_);
        binright = ((input_[k].x+placement_.module(k).width()-placement_.boundryLeft())/binSize_);
        bintop = ((input_[k].y+placement_.module(k).height()-placement_.boundryBottom())/binSize_);
        binbottom = ((input_[k].y-placement_.boundryBottom())/binSize_);
        if((binleft-3)>0)nearbinleft=binleft-3;
        else nearbinleft=0;
        if((binright+3)<widthBinNum_)nearbinright=binright+3;
        else nearbinright=widthBinNum_;
        if((bintop+3)<heightBinNum_)nearbintop=bintop+3;
        else nearbintop=heightBinNum_;
        if((binbottom-3)>0)nearbinbottom=binbottom-3;
        else nearbinbottom=0;
        
        /*if((binleft)>0)nearbinleft=binleft;
        else nearbinleft=0;
        if((binright)<widthBinNum_)nearbinright=binright;
        else nearbinright=widthBinNum_-1;
        if((bintop)<heightBinNum_)nearbintop=bintop;
        else nearbintop=heightBinNum_-1;
        if((binbottom)>0)nearbinbottom=binbottom;
        else nearbinbottom=0;*/
        double mWidth = placement_.module(k).width(), mHeight = placement_.module(k).height();
        double alphax = 4/((mWidth+2*binSize_)*(mWidth+4*binSize_));
        double betax = 2/(binSize_*(mWidth+4*binSize_));
        double alphay = 4/((mHeight+2*binSize_)*(mHeight+4*binSize_));
        double betay = 2/(binSize_*(mHeight+4*binSize_));
        for(int i=nearbinleft;i<nearbinright;i++){
            for(int j=nearbinbottom;j<nearbintop;j++){
                
                signdx=(input_[k].x+mWidth/2)-(placement_.boundryLeft()+i*binSize_+0.5 * binSize_);
                signdy=(input_[k].y+mHeight/2)-(placement_.boundryBottom()+j*binSize_+0.5 * binSize_);
                dx=abs(signdx);
                dy=abs(signdy);
                double xl,yl,partialxl,partialyl =0;
                if( (dx>=0) && dx<(mWidth/2+binSize_)){
                    xl=1-alphax*dx*dx;
                }
                else if( (dx>(mWidth/2+binSize_)) && (dx<(mWidth/2+2*binSize_))){
                    xl=betax*(dx-mWidth/2-2*binSize_)*(dx-mWidth/2-2*binSize_);
                }
                else{
                    xl=0;
                }
                if( (dy>=0) && dy<(mHeight/2+binSize_)){
                    yl=1-alphay*dy*dy;
                }
                else if( (dy>(mHeight/2+binSize_)) && (dy<(mHeight/2+2*binSize_))){
                    yl=betay*(dy-mHeight/2-2*binSize_)*(dy-mHeight/2-2*binSize_);
                }
                else{
                    yl=0;
                }
                if( (dx>=0) && (dx<(mWidth/2+binSize_))){
                    if(signdx>0)partialxl=-2*alphax*dx;
                    else partialxl=2*alphax*dx;
                    //partialxl=-2*alpha*dx;
                }
                else if( ((dx>(mWidth/2+binSize_))) && (dx<(mWidth/2+2*binSize_))){
                    if(signdx>0)partialxl=2*betax*(dx-mWidth/2-2*binSize_);
                    else partialxl=-2*betax*(dx-mWidth/2-2*binSize_);
                    //partialxl=2*beta*(dx-mWidth/2-2*binSize_);
                }
                else{
                    partialxl=0;
                }
                if(dy>=0 && (dy<(mHeight/2+binSize_))){
                    if(signdy>0)partialyl=-2*alphay*dy;
                    else partialyl=2*alphay*dy;
                    //partialyl=-2*alpha*dy;
                }
                else if((dy>(mHeight/2+binSize_)) && (dy<(mHeight/2+2*binSize_))){
                    if(signdy>0)partialyl=2*betay*(dy-mHeight/2-2*binSize_);
                    else partialyl=-2*betay*(dy-mHeight/2-2*binSize_);
                    //partialyl=2*beta*(dy-mHeight/2-2*binSize_);
                }
                else{
                    partialyl=0;
                }
                //grad_[k].x += 2*(binDensity[i][j]-Mb)*partialxl*yl*mHeight/mWidth;
                grad_[k].x += 2*(binDensity[i][j]-Mb)*partialxl*yl;
                grad_[k].y += 2*(binDensity[i][j]-Mb)*xl*partialyl;
            }
        }
    }
    /*for(int i=0;i<widthBinNum_;i++){
        for(int j=0;j<heightBinNum_;j++){
            for(unsigned k=0;k<placement_.numModules();k++){
                
                double mWidth = placement_.module(k).width(), mHeight = placement_.module(k).height();
                double alpha = 4/((mWidth+2*binSize_)*(mWidth+4*binSize_));
                double beta = 2/(binSize_*(mWidth+4*binSize_));
                dx=abs((input_[k].x+mWidth/2)-(placement_.boundryLeft()+i*binSize_+0.5 * binSize_));
                dy=abs((input_[k].y+mHeight/2)-(placement_.boundryBottom()+j*binSize_+0.5 * binSize_));
                double xl,yl,partialxl,partialyl =0;
                if( (dx>=0) && dx<(mWidth/2+binSize_)){
                    xl=1-alpha*dx*dx;
                }
                else if( (dx>(mWidth/2+binSize_)) && dx<(mWidth/2+2*binSize_)){
                    xl=beta*(dx-mWidth/2-2*binSize_)*(dx-mWidth/2-2*binSize_);
                }
                else{
                    xl=0;
                }
                if(dy>=0 && dy<(mHeight/2+binSize_)){
                    yl=1-alpha*dy*dy;
                }
                else if(dy>(mHeight/2+binSize_) && dy<(mHeight/2+2*binSize_)){
                    yl=beta*(dy-mHeight/2-2*binSize_)*(dy-mHeight/2-2*binSize_);
                }
                else{
                    yl=0;
                }
                if( (dx>=0) && dx<(mWidth/2+binSize_)){
                    partialxl=-2*alpha*dx;
                }
                else if( (dx>(mWidth/2+binSize_)) && dx<(mWidth/2+2*binSize_)){
                    partialxl=2*beta*(dx-mWidth/2-2*binSize_);
                }
                else{
                    partialxl=0;
                }
                if(dy>=0 && dy<(mHeight/2+binSize_)){
                    partialyl=-2*alpha*dy;
                }
                else if(dy>(mHeight/2+binSize_) && dy<(mHeight/2+2*binSize_)){
                    partialyl=2*beta*(dy-mHeight/2-2*binSize_);
                }
                else{
                    partialyl=0;
                }
                grad_[k].x += 2*(binDensity[i][j]-Mb)*partialxl*yl;
                grad_[k].y += 2*(binDensity[i][j]-Mb)*xl*partialyl;
            }
        }
    }*/
    return grad_;
}



const double &ObjectiveFunction::operator()(const std::vector<Point2<double>> &input)
{
    // Compute the value of the function
    //lambda = 1;
    lastvalue = currentvalue;
    double wlcost =wirelength_(input);
    double densitycost = density_(input);
    //cout<<"wirelength: "<<wlcost<<endl;
    //cout<<"density: "<<densitycost<<endl;
   value_=wlcost+lambda*densitycost;
    currentvalue = value_;
    //cout<<"value: "<<currentvalue<<" lastvalue: "<<lastvalue<<" difference: "<<currentvalue-lastvalue<<endl;
    
    input_ = input;
    return value_;
}

const std::vector<Point2<double>> &ObjectiveFunction::Backward()
{
    cout<<"lambda: "<<lambda<<" iter "<<iterNum_<<endl;
    vector<Point2<double>> grad1 = wirelength_.Backward();
    vector<Point2<double>> grad2 = density_.Backward();
    if (this->getOverflowRatio() < 0.2)
    {
        spreadEnough_ = true;
    }
    if (iterNum_ == 30)
    {
        double wirelengthGradNorm = 0, densityGradNorm = 0;
        for (unsigned i = 0; i <placement_.numModules(); ++i)
        {
            wirelengthGradNorm += sqrt(grad1[i].x * grad1[i].x + grad1[i].y * grad1[i].y);
            densityGradNorm += sqrt(grad2[i].x * grad2[i].x + grad2[i].y * grad2[i].y);
        }
        lambda = wirelengthGradNorm / densityGradNorm * 0.8;
    }
    bool ifproduct=1;
    if(lambda>100000000000000000000000000000000000000000000000000.)ifproduct=0;
    if((currentvalue-lastvalue)>=0 && ifproduct)
    {

        lambda *= 1.1;
        //cout<<"di>=0: "<<(currentvalue-lastvalue)<<endl;
    }
    for(unsigned i=0;i<placement_.numModules();i++){
        grad2[i].x = grad2[i].x*lambda;
        grad2[i].y = grad2[i].y*lambda;
    }
    for(unsigned i=0;i<placement_.numModules();i++){
        grad_[i]=grad1[i]+grad2[i];
    }
    //cout<<"gradwl: "<<grad1[0].x<<" "<<grad1[0].y<<" graddt: "<<grad2[0].x<<" "<<grad2[0].y<<endl;
    iterNum_++;
    return grad_;
    
}

