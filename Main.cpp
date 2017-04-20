#include <string>
#include <memory>
#include <vector>
#include <iostream>
#include <cassert>
#include <stack>

#include "pugixml.hpp"
#include "Distribution.h"
#include "Reaction.h"
#include "Nuclide.h"
#include "Material.h"
#include "Surface.h"
#include "Cell.h"
#include "Source.h"
#include "readXML.h"


int main(){
    
    //read xml file---------------------------------------------------------------------------
    unsigned long long                                      nsamples;
    double                                                  tCut;
    std::vector < std::shared_ptr<surface>     >            surfaces;
    std::vector < std::shared_ptr<cell>        >            cells;
    std::vector < std::shared_ptr<nuclide>     >            nuclides;
    std::vector < std::shared_ptr<material>    >            materials;
    std::vector < std::shared_ptr<estimator>   >            estimators;
    std::shared_ptr< source >                               src;
    std::vector < std::shared_ptr<distribution<double>> >   double_distributions;
    std::vector < std::shared_ptr<distribution<int>>    >   int_distributions;
    std::vector < std::shared_ptr<distribution<point>>  >   point_distributions;

    readXML
    ( nsamples, tCut, surfaces, cells, nuclides, materials, estimators, src, double_distributions, int_distributions, point_distributions );
    std::cout.flush();

    //transport start--------------------------------------------------------------------------
    int cellNum;
    unsigned long long numMove;
    
    for ( unsigned long long ientry = 0 ; ientry < nsamples ; ientry++ ) {
        
        std::stack<particle>pbank = src->sample();
        
        while( ! pbank.empty() ){
            particle p = pbank.top() ; // copy p to working particle
            pbank.pop();               // pop particle off bank
            
            while ( p.alive() ) {
                
                //where are we? which cell?
                for(int j = 0 ; j<cells.size() ; j++){
                    bool whichCell = cells[j]->testPoint( p.pos() );
                    if(whichCell ==  true){cellNum = j;}
                }
                double impNow = cells[cellNum]->getImportance();
                
                if(impNow > 0.0){
                    
                    //get distance to collision. need to get macrocopic xs of that cell
                    double sigma = cells[cellNum]->macro_xs( &p );
                    double dist_collision = -std::log( Urand() ) / sigma;
                    
                    //get distance to the closest surface
                    std::pair< std::shared_ptr< surface >, double >
                    S = cells[cellNum]->surfaceIntersect( p.getRay() );
                    double dist_surface = S.second;
                    
                    // advance the particle
                    // note: moveParticle calls the cell estimator
                    double distance = std::fmin( dist_collision, dist_surface );
                    cells[cellNum]->moveParticle( &p, distance, tCut  );
                    numMove = numMove + 1;
                    
                    // process event
                    if ( distance == dist_surface ) {
                        
                        // cross surface, note that this calls the surface estimator
                        S.first->crossSurface( &p );
                        
                        //check where is it now and what's the importance
                        for(int j = 0 ; j<cells.size() ; j++){
                            bool whichCell = cells[j]->testPoint( p.pos() );
                            if(whichCell ==  true){cellNum = j;}
                        }
                        double impNew = cells[cellNum]->getImportance();
                        
                        //apply variance reduction and get the extra particles if they exist
                        std::stack < particle> tempBank = S.first->varianceReduction( &p, impNow, impNew );
                        while( ! tempBank.empty() ){
                            particle par = tempBank.top() ;
                            tempBank.pop();
                            pbank.push ( par );
                        }
                    }
                    else {
                        //the particle undergoes reaction
                        cells[cellNum]->sampleCollision( &p, &pbank);
                    }
                }
                //the impNow of the current cell is 0.0, kill it
                else{
                    p.kill();
                }
            } //particle is alive?
        }//bank is not empty?
        
        // closeout all estimators
        for ( auto e : estimators ) { e->endHistory(); }
        
    }// history loop
    
    std::cout << "Results: " << std::endl;
    for ( auto e : estimators ) { e->report( numMove ); }
    
    return 0;
}

