#ifndef _DISTRIBUTION_HEADER_
#define _DISTRIBUTION_HEADER_

#include <iostream>
#include <utility>
#include <vector>
#include <cmath>
#include <cassert>
#include <string>

#include "Random.h"
#include "Point.h"

template <class T>
class distribution {
private:
    std::string distribution_name;
protected:
    
public:
    distribution( std::string label ) : distribution_name(label) {};          // constructor
    ~distribution() {};   // destructor
    
    virtual std::string name() final { return distribution_name; };
    virtual T sample() = 0;  // dummy function that must be implemented in each case
};

template <class T>
class arbitraryDelta_distribution : public distribution<T> {
private:
    T result;
public:
    arbitraryDelta_distribution( std::string label, T val ) : distribution<T>(label), result(val) {};
    ~arbitraryDelta_distribution() {};
    
    T sample() { return result; }
};

template <class T>
class arbitraryDiscrete_distribution : public distribution<T> {
private:
    std::vector< std::pair< T, double > > cdf;
public:
    arbitraryDiscrete_distribution( std::string label, std::vector< std::pair< T, double > > data );
    ~arbitraryDiscrete_distribution() {};
    T sample();
};

template < class T >
arbitraryDiscrete_distribution<T>::arbitraryDiscrete_distribution( std::string label, std::vector< std::pair< T, double > > data )
: distribution<T>(label) {
    // convert pdf to cdf
    double c = 0.0;
    for ( auto d : data ) {
        // first is pointer to data type T and second is pdf input
        cdf.push_back( std::make_pair( d.first, d.second + c ) );
        c += d.second;
    }
}

template < class T >
T arbitraryDiscrete_distribution<T>::sample() {
    double   r = Urand() * cdf.back().second;
    for ( auto c : cdf ) {
        // first is pointer to data type T and second is cdf
        if ( r < c.second ) { return c.first; };
    }
    assert( false ); // should not get here
    return cdf.back().first;
}

class delta_distribution : public distribution<double> {
private:
    double a;
public:
    delta_distribution( std::string label, double p1 ) : distribution(label), a(p1) {};
    ~delta_distribution() {};
    double sample() { return a; };
};

class uniform_distribution : public distribution<double> {
private:
    double a, b;
public:
    uniform_distribution( std::string label, double p1, double p2 ) : distribution(label), a(p1), b(p2) {};
    ~uniform_distribution() {};
    double sample();
};

class linAnisotropicDirection_distribution : public distribution<double> {
private:
    double muBar;
public:
    linAnisotropicDirection_distribution( std::string label , double p1 ) : distribution(label ), muBar(p1) {};
    ~linAnisotropicDirection_distribution() {};
    double sample();
};

class linear_distribution : public distribution<double> {
private:
    double a, b, fa, fb;
public:
    linear_distribution( std::string label, double x1, double x2, double y1, double y2 )
    : distribution(label), a(x1), b(x2), fa(y1), fb(y2) {};
    ~linear_distribution() {};
    double sample();
};

class exponential_distribution : distribution<double> {
private:
    double lambda;
public:
    exponential_distribution( std::string label, double p1 ) : distribution(label), lambda(p1) {};
    ~exponential_distribution() {};
    double sample();
};

class normal_distribution : public distribution<double> {
private:
    const double twopi = 2.0 * std::acos(-1.0);
    double mu, sigma;
public:
    normal_distribution( std::string label, double p1, double p2 ) : distribution(label), mu(p1), sigma(p2) {};
    ~normal_distribution() {};
    double sample();
};

class HenyeyGreenstein_distribution : public distribution<double> {
private:
    double a;
public:
    HenyeyGreenstein_distribution( std::string label, double p1 ) : distribution(label), a(p1) {};
    ~HenyeyGreenstein_distribution() {};
    double sample();
};

class meanMultiplicity_distribution : public distribution<int> {
private:
    double nu;
public:
    meanMultiplicity_distribution( std::string label, double p1 ) : distribution(label), nu(p1) {};
    ~meanMultiplicity_distribution() {};
    int sample();
};

class TerrellFission_distribution : public distribution<int> {
private:
    double nubar, sigma, b;
    
    std::vector<double> cdf;
public:
    TerrellFission_distribution( std::string label, double p1, double p2, double p3 );
    ~TerrellFission_distribution() {};
    int sample();
};

class isotropicDirection_distribution : public distribution<point> {
private:
    const double twopi = 2.0 * std::acos(-1.0);
public:
    isotropicDirection_distribution( std::string label ) : distribution(label) {};
    ~isotropicDirection_distribution() {};
    point sample();
};

class forwardPeakDirection_distribution : public distribution<point> {
private:
    const double twopi = 2.0 * std::acos(-1.0);
public:
    forwardPeakDirection_distribution( std::string label ) : distribution(label) {};
    ~forwardPeakDirection_distribution() {};
    point sample();
};

class anisotropicDirection_distribution : public distribution<point> {
private:
    const double twopi = 2.0 * std::acos(-1.0);
    double sin_t;
    
    point axis;
    std::shared_ptr< distribution<double> > dist_mu;
public:
    anisotropicDirection_distribution( std::string label, point p, std::shared_ptr< distribution<double> > dmu )
    : distribution(label), axis(p), dist_mu(dmu)
    { axis.normalize(); sin_t = std::sqrt( 1.0 - axis.z * axis.z ); };
    ~anisotropicDirection_distribution() {};
    point sample();
};

class independentXYZ_distribution : public distribution<point> {
private:
    std::shared_ptr< distribution<double> > dist_x, dist_y, dist_z;
public:
    independentXYZ_distribution( std::string label, std::shared_ptr< distribution<double> > dx,
                                std::shared_ptr< distribution<double> > dy, std::shared_ptr< distribution<double> > dz )
    : distribution(label), dist_x(dx), dist_y(dy), dist_z(dz) {};
    ~independentXYZ_distribution() {};
    
    point sample();
};

class sampleThreePoint_distribution : public distribution<point> {
private:
    std::shared_ptr< distribution<double> > dist_1x, dist_2x, dist_3x;
    std::shared_ptr< distribution<double> > dist_1y, dist_2y, dist_3y;
    std::shared_ptr< distribution<double> > dist_1z, dist_2z, dist_3z;
    double prob1, prob2, prob3;
public:
    sampleThreePoint_distribution( std::string label, std::shared_ptr< distribution<double> > d1x, std::shared_ptr< distribution<double> > d1y, std::shared_ptr< distribution<double> > d1z,
                                  std::shared_ptr< distribution<double> > d2x, std::shared_ptr< distribution<double> > d2y, std::shared_ptr< distribution<double> > d2z, std::shared_ptr< distribution<double> > d3x, std::shared_ptr< distribution<double> > d3y, std::shared_ptr< distribution<double> > d3z, double p1, double p2, double p3 )
    : distribution(label), dist_1x(d1x), dist_1y(d1y), dist_1z(d1z), dist_2x(d2x), dist_2y(d2y), dist_2z(d2z), dist_3x(d3x), dist_3y(d3y), dist_3z(d3z), prob1(p1), prob2(p2), prob3(p3) {};
    ~sampleThreePoint_distribution() {};
    
    point sample();
};

class sampleHollowSphere_distribution : public distribution<point> {
private:
    double       inRadius, outRadius;
    const double twopi = 2.0 * std::acos(-1.0);
public:
    sampleHollowSphere_distribution( std::string label, double r1, double r2)
    : distribution(label), inRadius(r1), outRadius(r2) {};
    ~sampleHollowSphere_distribution() {};
    
    point sample();
};

class diskPositionZ_distribution : public distribution<point> {
private:
    double       oriX, oriY, oriZ, rad;
    const double twopi = 2.0 * std::acos(-1.0);
public:
    diskPositionZ_distribution( std::string label, double a, double b, double c, double d)
    : distribution(label), oriX(a), oriY(b), oriZ(c), rad(d) {};
    ~diskPositionZ_distribution() {};
    
    point sample();
};

class diskPositionX_distribution : public distribution<point> {
private:
    double       oriX, oriY, oriZ, rad;
    const double twopi = 2.0 * std::acos(-1.0);
public:
    diskPositionX_distribution( std::string label, double a, double b, double c, double d)
    : distribution(label), oriX(a), oriY(b), oriZ(c), rad(d) {};
    ~diskPositionX_distribution() {};
    
    point sample();
};


#endif
