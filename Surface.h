#ifndef _SURFACE_HEADER_
#define _SURFACE_HEADER_

#include <string>
#include <vector>
#include <limits>
#include <stack>

#include "Point.h"
#include "Particle.h"
#include "Estimator.h"
//#include "Source.h"

class surface {
private:
    std::string surface_name;
    bool        reflect_bc;
    std::vector< std::shared_ptr< estimator > > surface_estimators;
public:
    surface( std::string label ) : surface_name(label) { reflect_bc = false; };
    ~surface() {};
    
    virtual std::string name()    final { return surface_name; };
    virtual void makeReflecting() final { reflect_bc = true; };
    
    virtual void attachEstimator( std::shared_ptr< estimator > E ) final {
        surface_estimators.push_back( E );
    }
    virtual void scoreEstimators( particle* p ) final {
        for ( auto e : surface_estimators ) { e->score( p, 0.0, 0.0, 0.0, 0.0 ); }
    }
    
    virtual void crossSurface( particle* p ) final {
        
        // score estimators
        for ( auto e : surface_estimators ) { e->score( p, 0.0, 0.0, 0.0, 0.0 ); }
        
        // reflect if needed
        if ( reflect_bc ) {
            point d = reflect( p->getRay() );
            p->setDirection( d );
        }
        
        // advance particle off the surface
        p->move( std::numeric_limits<float>::epsilon() );
        
    }
    
    virtual std::stack<particle> varianceReduction( particle *p, double impNow, double impNew ) final{
        
        std::stack<particle> tempBank;
        
        //apply variance reduction
        if(impNow > impNew){
            
            //do roulette and bump the weight
            if(Urand() < impNew/impNow){
                p->adjustWeight( impNow/impNew );
            }
            else{
                p->kill();
            }
        }
        else if(impNow < impNew){
            
            //split the particle and reduce the weight
            int numSplit = floor(impNew/impNow + Urand());
            
            //keep one particle to work and keep the rest in the bank
            //also we need to reduce their weight
            for(int jentry = 0 ; jentry< (numSplit-1) ; jentry++){
                particle q(  p->pos(), p->dir(), p->energy()  ); //create new particles with weight 1
                q.adjustWeight( p->wgt() / numSplit  ); //multiply with p.wgt() bc they are created with weight 1
                tempBank.push( q );
            }
            
            //adjust the weight of this one working particle
            p->adjustWeight( 1.0 / numSplit );
        }
        else{
            //if impNow == impNew, do nothing
        }
        
        return tempBank;
    }
    
    virtual double eval( point p )   = 0;
    virtual double distance( ray r ) = 0;
    virtual point  reflect( ray r )  = 0;
};

class plane : public surface {
private:
    double a, b, c, d;
public:
    plane( std::string label, double p1, double p2, double p3, double p4 ) :
    surface(label), a(p1), b(p2), c(p3), d(p4) {};
    ~plane() {};
    
    double eval( point p );
    double distance( ray r );
    point  reflect( ray r );
};

class sphere : public surface {
private:
    double x0, y0, z0, rad;
public:
    sphere( std::string label, double p1, double p2, double p3, double p4 ) :
    surface(label), x0(p1), y0(p2), z0(p3), rad(p4) {};
    ~sphere() {};
    
    double eval( point p );
    double distance( ray r );
    point  reflect( ray r );
};

class cylinderx : public surface {
private:
    double y0, z0, rad;
public:
    cylinderx( std::string label, double p1, double p2, double p3 ) :
    surface(label), y0(p1), z0(p2), rad(p3) {};
    ~cylinderx() {};
    
    double eval( point p );
    double distance( ray r );
    point  reflect( ray r );
};

class cylinderz : public surface {
private:
    double x0, y0, rad;
public:
    cylinderz( std::string label, double p1, double p2, double p3 ) :
    surface(label), x0(p1), y0(p2), rad(p3) {};
    ~cylinderz() {};
    
    double eval( point p );
    double distance( ray r );
    point  reflect( ray r );
};

class conex : public surface {
private:
    double x0, y0, z0, aDenom, bDenom, cDenom;
    double a2, b2, c2;
public:
    conex( std::string label, double p1, double p2, double p3, double p4, double p5, double p6 ) :
    surface(label), x0(p1), y0(p2), z0(p3), aDenom(p4), bDenom(p5), cDenom(p6) {};
    ~conex() {};
    
    double eval( point p );
    double distance( ray r );
    point  reflect( ray r );
};


#endif
