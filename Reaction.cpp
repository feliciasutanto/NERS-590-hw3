#include "Reaction.h"
#include "Particle.h"
#include "Distribution.h"
#include <iostream>

void  capture_reaction::sample( particle* p, std::stack< particle >* bank, double A ) {
    p->kill();
}

void  scatter_reaction::sample( particle* p, std::stack< particle >* bank, double A ) {
    double mu0 = scatter_dist->sample();
    p->scatter( mu0, A );
}

void  fission_reaction::sample( particle* p, std::stack< particle >* bank, double A ) {
    int n = multiplicity_dist->sample();
    if ( n <= 0 ) {
        p->kill();
    }
    else {
        // bank all but last particle (skips if n = 1)
        for ( int i = 0 ; i < (n - 1) ; i++ ) {
            particle q( p->pos(), isotropic->sample(), p->energy() ); //energy dist for fission? watt?
            q.recordCell( p->cellPointer() );
            q.adjustWeight(p->wgt()); //this new particle should have the same weight as its parent
            bank->push( q );
        }
        // set working particle to last one
        particle q( p->pos(), isotropic->sample(), p->energy() ); //energy dist for fission? watt?
        q.recordCell( p->cellPointer() );
        q.adjustWeight(p->wgt());
        *p = q;
    }
}

double capture_reaction::xs( particle *p ){
    double rxn_xs = xsa + xsb / std::sqrt( p->energy() ); // E in MeV
    return rxn_xs;
}

double scatter_reaction::xs( particle *p ){
    double rxn_xs = xsa + xsb / std::sqrt( p->energy() ); // E in MeV
    return rxn_xs;
}

double fission_reaction::xs( particle *p ){
    double rxn_xs = xsa + xsb / std::sqrt( p->energy() ); // E in MeV
    return rxn_xs;
}
