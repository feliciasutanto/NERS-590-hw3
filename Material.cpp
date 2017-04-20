#include <vector>
#include <utility>
#include <memory>
#include <cassert>
#include <iostream>

#include "Random.h"
#include "Particle.h"
#include "Reaction.h"
#include "Nuclide.h"
#include "Material.h"

// add a nuclide and atom fraction to the material
void material::addNuclide( std::shared_ptr< nuclide > N, double frac ) {
    nuclides.push_back( std::make_pair( N, frac ) );
}

// private utility function that returns sum of atomic fraction * microscopic total xs
// multiply this by atomic density to get macroscopic cross section
// (this function is useful to accelerate sampling the nuclide)
double material::micro_xs( particle* p ) {
    double xs = 0.0;
    for ( auto n : nuclides ) {
        // first is pointer to nuclide, second is atomic fraction
        xs += n.first->total_xs( p ) * n.second;
    }
    return xs;
}

// return the macroscopic cross section
double material::macro_xs( particle* p ) {
    return atom_density() * micro_xs( p );
}

//----------------------------------------------------------------
//return the material's total capture microscopic cross section
double material::microCap_xs( particle* p ) {
    double xs = 0.0;
    for ( auto n : nuclides ) {
        // first is pointer to nuclide, second is atomic fraction
        xs += n.first->cap_xs( p ) * n.second;
    }
    return xs;
}

//return the macroscopic capture cross section
double material::macroCap_xs( particle* p ) {
    return atom_density() * microCap_xs( p );
}

//----------------------------------------------------------------
//return the material's total scatter microscopic cross section
double material::microScat_xs( particle* p ) {
    double xs = 0.0;
    for ( auto n : nuclides ) {
        // first is pointer to nuclide, second is atomic fraction
        xs += n.first->scat_xs( p ) * n.second;
    }
    return xs;
}

//return the macroscopic scatter cross section
double material::macroScat_xs( particle* p ) {
    return atom_density() * microScat_xs( p );
}

//---------------------------------------------------------------- 

// randomly sample a nuclide based on total cross sections and atomic fractions
std::shared_ptr< nuclide > material::sample_nuclide( particle* p ) {
    double u = micro_xs( p ) * Urand();
    double s = 0.0;
    for ( auto n : nuclides ) {
        // first is pointer to nuclide, second is atomic fraction
        s += n.first->total_xs( p ) * n.second;
        if ( s > u ) { return n.first; }
    }
    assert( false ); // should never reach here
    return nullptr;
}

// function that samples an entire collision: sample nuclide, then its reaction,
// and finally process that reaction with input pointers to the working particle p
// and the particle bank

void material::sample_collision( particle* p, std::stack<particle>* bank) {
    
    // first sample nuclide
    std::shared_ptr< nuclide >  N = sample_nuclide( p );
    
    // now get the reaction
    std::shared_ptr< reaction > R = N->sample_reaction( p );
    
    // finally process the reaction
    R->sample( p, bank, N->getA() );
}
