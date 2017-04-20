#ifndef _REACTION_HEADER_
#define _REACTION_HEADER_

#include <vector>
#include <iostream>
#include <string>
#include <memory>
#include <stack>
#include <utility>

#include "Particle.h"
#include "Distribution.h"

class reaction {
private:
protected:
    std::string rxn_name;
    double rxn_xs;
    double xsa;
    double xsb;
public:
    reaction( double xa, double xb ) : xsa(xa), xsb(xb) {};
    ~reaction() {};
    
    virtual std::string name() final { return rxn_name; };
    virtual double xs( particle *p ) = 0; //final { return rxn_xs; };
    virtual void sample( particle* p, std::stack<particle>* bank, double A ) = 0;
};

class capture_reaction : public reaction {
private:
public:
    capture_reaction( double xa, double xb ) : reaction(xa, xb) { rxn_name = "capture"; } ;
    ~capture_reaction() {};
    
    double xs( particle *p );
    void sample( particle* p, std::stack<particle>* bank, double A );
};

class scatter_reaction : public reaction {
private:
    std::shared_ptr< distribution<double> > scatter_dist;
public:
    scatter_reaction( double xa, double xb, std::shared_ptr< distribution<double> > D ) :
    reaction(xa, xb) , scatter_dist(D) { rxn_name = "scatter"; };
    ~scatter_reaction() {};
    
    double xs( particle *p );
    void sample( particle* p, std::stack<particle>* bank, double A );
};

class fission_reaction : public reaction {
private:
    std::shared_ptr< distribution<int> >   multiplicity_dist;
    std::shared_ptr< distribution<point> > isotropic;
public:
    fission_reaction( double xa, double xb, std::shared_ptr< distribution<int> > D ) :
    reaction(xa, xb) , multiplicity_dist(D) {
        rxn_name = "fission";
        isotropic = std::make_shared< isotropicDirection_distribution > ( "isotropic" );
    };
    ~fission_reaction() {};
    
    double xs( particle *p );
    void sample( particle* p, std::stack<particle>* bank, double A );
};

#endif
