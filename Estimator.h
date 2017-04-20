#ifndef _ESTIMATOR_HEADER_
#define _ESTIMATOR_HEADER_

#include <cmath>
#include <iostream>
#include <string>
#include <vector>
#include <cassert>
#include <typeinfo>

#include "Particle.h"
#include "Material.h"
#include "Reaction.h"

//binning class---------------------------------------------------------------------
class bin{
    
public:
    bin() {};
    ~bin() {};
    
    std::vector <double> pathLengthBin(particle*, double, double, std::vector<double>, double , double, int, double) ;
    
};

// base estimator class------------------------------------------------------------
class estimator {
private:
    std::string estimator_name;
protected:
    unsigned long long nhist;
public:
    estimator( std::string label, std::string name ) : estimator_name(label)  {};
    ~estimator() {};
    
    //class to say name of the estimator
    virtual std::string name() final { return estimator_name; };
    
    //class to score
    virtual void score( particle*, double, double, double, double ) = 0;
    template< typename T >
    // particle, distance travelled, macro cap xs, cell name
    void score( particle*, double,  double, double, double,  T ) { assert(false); };
    void score( particle*, double,  double, double, double, std::shared_ptr< material > )  {};
    
    //class to end history
    virtual void endHistory()       = 0;
    
    //class to report the result of the estimator
    virtual void report( int )      = 0;
    
};

//counting_estimator---------------------------------------------------------------
class counting_estimator : public estimator {
private:
    int count_hist;
    std::vector< double > tally;
    std::string theSurface;
    std::string theTallyType;
public:
    counting_estimator( std::string label, std::string nameS ) : estimator(label, nameS ), theSurface(nameS), theTallyType(label) { count_hist = 0; };
    ~counting_estimator() {};
    
    void score( particle*, double, double, double, double);
    void endHistory();
    void report(int);
};

// derived class for simple estimators like current or scalar flux--------------------
// future estimators could get spectra or flux on a mesh
class single_valued_estimator : public estimator {
private:
    
protected:
    std::vector <double> tally_hist_vec;
    std::vector <double> tally_sum_vec;
    std::vector <double> tally_squared_vec;
    int numBin = 0;
    
public:
    using estimator::score;
    
    single_valued_estimator(std::string label, int number ) : estimator(label, "dummy" ), numBin(number) {
        
        //initialization
        nhist              = 0;
        for(int a = 0; a<numBin ; a++){
            tally_hist_vec.push_back(0.0);
            tally_sum_vec.push_back(0.0);
            tally_squared_vec.push_back(0.0);
        }
        
    };
    ~single_valued_estimator() {};
    
    virtual void endHistory()    final {
        nhist++;
        for(int b = 0 ; b<numBin ; b++){
            tally_sum_vec[b]     += tally_hist_vec[b];
            tally_squared_vec[b] += tally_hist_vec[b] * tally_hist_vec[b];
            tally_hist_vec[b]     = 0.0;
        }
        
    }
    
    virtual void score( particle*, double, double, double, double) = 0;
    virtual void report(int) = 0;
    
};

// current crossing a surface-----------------------------------------------------
class surface_current_estimator : public single_valued_estimator {
private:
    double lowBin;
    double dBin;
    int    numBin;
    std::vector <double> myTimeVector;
    std::vector <double> mean_vec;
    std::vector <double> var_vec;
    std::string theSurface;
    std::string theTallyType;
    
public:
    surface_current_estimator( std::string label, int number, double lBin, double Bin, int nBin, std::string nameS ) : single_valued_estimator(label , number ), lowBin(lBin), dBin(Bin), numBin(nBin), theSurface(nameS), theTallyType(label) {};
    ~surface_current_estimator() {};
    
    void score( particle*, double, double, double, double );
    
    void report(int) {
        
        for(int c = 0 ; c<numBin ; c++){
            double mean = tally_sum_vec[c] / nhist;
            mean_vec.push_back(  mean  );
            var_vec.push_back( (tally_squared_vec[c] / nhist - mean*mean )/ nhist);
        }
        
        //create the time vector
        for(int d = 0 ; d<numBin ; d++){
            myTimeVector.push_back(d * dBin + (dBin/2.0) + lowBin );
        }
        
        //output the results
        std::cout <<"===================="<<theTallyType<<"===================="<<std::endl;
        std::cout << "Time [ns] | Prob |  RelativeErr |  FOM" <<std::endl;
        for(int e = 0 ;e<numBin ; e++){
            std::cout <<myTimeVector[e]<<"     "<<mean_vec[e]<<"     "<< std::sqrt(var_vec[e]) / mean_vec[e] << std::endl;
        }
        
    };
    
};

// volume-averaged scalar flux in a cell--------------------------------------------
class cell_pathLengthFlux_estimator : public single_valued_estimator {
private:
    double volume; //we don't use the volume
    double lowBin;
    double dBin;
    int    numBin;
    std::string type;
    std::vector <double> myTimeVector;
    std::vector <double> mean_vec;
    std::vector <double> var_vec;
    double  macroXs;
    std::string theCell;
    std::string theTallyType;
    
public:
    cell_pathLengthFlux_estimator( std::string label, int nBin , double vol, double lBin, double deltaBin, std::string t, std::string nameC ):
    single_valued_estimator(label, nBin ) , volume(vol), lowBin(lBin), dBin(deltaBin), type(t), numBin(nBin), theCell(nameC), theTallyType(label) {};
    ~cell_pathLengthFlux_estimator() {};
    
    void score( particle* , double, double, double, double );
    
    void report(int numMove) {
        
        for(int c = 0 ; c<numBin ; c++){
            double mean = tally_sum_vec[c] / nhist;
            mean_vec.push_back(  mean  );
            var_vec.push_back( (tally_squared_vec[c] / nhist - mean*mean )/ nhist);
        }
        
        //create the time vector
        for(int d = 0 ; d<numBin ; d++){
            myTimeVector.push_back(d * dBin + (dBin/2.0) + lowBin );
        }
        
        //output the results
        std::cout <<"===================="<<theTallyType<<"===================="<<std::endl;
        std::cout << "Time[ns] | Prob-"<<type<<" | RelativeErr | FOM" <<std::endl;
        for(int e = 0 ;e<numBin ; e++){
            std::cout <<myTimeVector[e]<<"     "<<mean_vec[e]<<"     "<< std::sqrt(var_vec[e]) / mean_vec[e] <<"       "<< 1.0/(numMove * pow( (std::sqrt( var_vec[e] ) / mean_vec[e] ) ,2) ) << std::endl;
        }
        
    };
};


#endif
