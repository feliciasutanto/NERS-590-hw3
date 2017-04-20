#include <iostream>
#include <vector>
#include <utility>
#include <memory>
#include <cassert>
#include <cmath>
#include <limits>

#include "Cell.h"
#include "Particle.h"

// take in a pair of surface pointer and integer describing sense (must not be zero!)
// and append to vector of surfaces
void cell::addSurface( std::shared_ptr< surface > S, int sense ) {
    
    assert( sense != 0 );
    int sgn = std::copysign( 1, sense );
    
    surfaces.push_back( std::make_pair( S, sgn ) );
}

// test if point p inside the current cell
bool cell::testPoint( point p ) {
    
    // loop over surfaces in cell, if not on correct side return false
    // if on correct side of all surfaces, particle is in the cell and return true
    for ( auto s : surfaces ) {
        // first = surface pointer, second = +/- 1 indicating sense
        if ( s.first->eval( p ) * s.second < 0 ) { return false; }
    }
    return true;
}

// find first intersecting surface of ray r and distance to intersection
std::pair< std::shared_ptr< surface >, double > cell::surfaceIntersect( ray r ) {
    
    double dist = std::numeric_limits<double>::max();
    std::shared_ptr< surface > S = nullptr;
    for ( auto s : surfaces ) {
        // first is a surface pointer; distance is always positive or huge if invalid
        double d = s.first->distance( r );
        if ( d < dist ) {
            // current surface intersection is closer
            dist = d;
            S    = s.first;
        }
    }
    return std::make_pair( S, dist );
}

void cell::moveParticle( particle* p, double s, double tCut ){
    double timeOld = p->time();
    double energyOld = p->energy();
    p->move( s );
    
    //update the time of the particle
    double dt = s / p->speed(); //speed in cm/ns
    p->updateTime(dt);          // unit is ns
    
    //pass the value of the macroCap_xs to the estimators along with the path length travelled
    for ( auto e : cell_estimators ){
        
        //get the material of this cell and then get the macro Capture xs and macro Scat xs
        double macroCapXs  = cell_material->macroCap_xs( p );
        double macroScatXs = cell_material->macroScat_xs( p );
        
        e->score( p, s , macroCapXs, macroScatXs, timeOld );

    }
    
    //check the lifetime of the particle, if it is already more than tCut [ns], just kill it
    if(p->time() > tCut){
        p->kill();
    }
}

void cell::sampleCollision( particle* p, std::stack<particle>* bank ) {
    cell_material->sample_collision( p, bank );
}

