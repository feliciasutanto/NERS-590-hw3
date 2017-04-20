#include <cmath>
#include <iostream>

#include "Estimator.h"
#include "Material.h"
#include "Particle.h"

//class to bin my tallies
std::vector < double> pathLengthBin(particle* p, double s, double macroXs, std::vector<double> myBin, double lowBin, double dBin, int numBin, double timeOld){
    
    //then find the index that correspond to timeOld and timeNew and their delta time
    int indexOld = 0;
    int indexNew = 0;
    double dtOld = 0.0;
    double dtNew = 0.0;
    
    for(int f = 0; f<numBin+1 ; f++){
        if(timeOld < f*dBin + lowBin ){
            indexOld = f-1;
            dtOld = f*dBin + lowBin - timeOld ;
            break;
        }
    }
    for(int g = 0; g<numBin+1 ; g++){
        if(p->time() < g*dBin + lowBin ){
            indexNew = g-1;
            dtNew = p->time() - ((g-1)*dBin + lowBin);
            break;
        }
    }
    
    if(indexNew != indexOld){
        
        //fill up the tally_hist_scat vector
        for(int h = indexOld; h<( indexNew+1 ) ; h++){
            if(h == indexOld){
                myBin[h] += p->speed() * dtOld * macroXs * p->wgt();
            }
            if(h > indexOld && h < indexNew ){
                myBin[h] += p->speed() * dBin  * macroXs * p->wgt();
            }
            if(h == indexNew){
                myBin[h] += p->speed() * dtNew * macroXs * p->wgt();
            }
        }
    }
    else{
        myBin[indexOld] += p->speed() * (p->time() - timeOld) * macroXs * p->wgt();
    }
    
    return myBin;
}


//counting_estimator--------------------------------------------------------------
void counting_estimator::score( particle* p, double null, double macroCapXs, double macroScatXs, double timeOld) {
    count_hist++;
}

void counting_estimator::endHistory() {
    if ( tally.size() < count_hist + 1 ) { tally.resize( count_hist + 1, 0.0 ); }
    tally[ count_hist ] += 1.0;
    nhist++;
    count_hist = 0;
}

void counting_estimator::report(int) {
    std::cout <<"===================="<<theTallyType<<"===================="<<std::endl;
    //std::cout << name() << std::endl;
    double s1 = 0.0, s2 = 0.0;
    for ( int i = 0 ; i < tally.size() ; i++ ) {
        double p = tally[i] / nhist;
        std::cout << i << " " << p << "   " << std::sqrt( p * ( 1.0 - p ) / nhist ) / p <<  std::endl;
        s1 += p * i;
        s2 += p * i * i;
    }
    std::cout << "   mean = " << s1 << std::endl;
    std::cout << "   var  = " << s2 - s1*s1 << std::endl;
}

//surface_current_estimator------------------------------------------------------
void surface_current_estimator::score( particle* p, double null, double macroCapXs, double macroScatXs, double timeOld) {
    
    //we are just need to count at the surface
    p->setSpeed(1.0);
    p->setTime(1.0);
    double macroXs = 1.0;
    timeOld        = 0.0;
    tally_hist_vec = pathLengthBin( p, 0.0 , macroXs, tally_hist_vec, lowBin, dBin, numBin, timeOld);
}

//cell_pathLengthFlux_estimator---------------------------------------------------
void cell_pathLengthFlux_estimator::score( particle* p , double path_length, double macroCapXs, double macroScatXs, double timeOld){
    
    //make a choice, scatter or absorption?
    if(type == "scatter"){
        macroXs = macroScatXs;
    }
    else if(type == "absorption"){
        macroXs = macroCapXs;
    }
    else{
        std::cout << "Such reaction has no macroscopic cross section"<< std::endl;
        std::cin.ignore(); //stop the code
    }
    
    tally_hist_vec = pathLengthBin( p, path_length, macroXs, tally_hist_vec, lowBin, dBin, numBin, timeOld);

}


