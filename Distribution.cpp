#include <string>
#include <cmath>

#include "Distribution.h"
#include "Point.h"
#include "Random.h"

double uniform_distribution::sample() { return a + Urand() * ( b - a ); }

double linAnisotropicDirection_distribution::sample() {
    double b   = 3.0 * muBar;
    double mu  = (   std::sqrt( b*b + 4.0*b*Urand() -2.0*b + 1.0   ) - 1.0  ) / b;
    return mu;
}

double linear_distribution::sample() {
    double r1 = Urand(), r2 = Urand();
    double p  = 2.0 * std::fmin( fa, fb ) / ( fa + fb );
    if ( r1 < p ) { return a + r2 * ( b - a ); }
    else {
        if ( fb > fa ) { return a + ( b - a ) * std::sqrt( r2 ); }
        else           { return a + ( b - a ) * ( 1.0 - std::sqrt( r2 )); }
    }
}

double exponential_distribution::sample() { return -std::log( Urand() ) / lambda; }

double normal_distribution::sample()
{ return mu + sigma * std::sqrt( -2.0 * std::log( Urand() ) ) * std::cos( twopi * Urand() ); }

double HenyeyGreenstein_distribution::sample() {
    if ( a != 0.0 ) {
        return ( 1.0 + a*a - pow( ( 1.0 - a*a )/( 1.0 + a*(2.0*Urand() - 1.0) ), 2 ) ) / ( 2.0 * a );
    }
    else {
        return 2.0 * Urand() - 1.0;
    }
}

int meanMultiplicity_distribution::sample() {
    return (int) std::floor( nu + Urand() );
}

TerrellFission_distribution::TerrellFission_distribution( std::string label, double p1, double p2, double p3 )
: distribution(label), nubar(p1), sigma(p2), b(p3) {
    double c  = 0.0;
    double nu = 0.0;
    while ( c < 1.0 - 1.0e-12 ) {
        double a  = ( nu - nubar + 0.5 + b ) / sigma;
        c = 0.5 * ( 1 + erf( a / sqrt(2.0) ) ) ;
        
        cdf.push_back(c);
        nu += 1.0;
    }
    cdf.push_back(1.0);
}

int TerrellFission_distribution::sample() {
    double r  = Urand();
    double nu;
    for ( int i = 0 ; i < cdf.size() ; i++ ) {
        if ( r < cdf[i] ) {
            nu = (double) i;
            break;
        }
    }
    return nu;
}

point isotropicDirection_distribution::sample() {
    // sample polar cosine and azimuthal angle uniformly
    double mu  = 2.0 * Urand() - 1.0;
    double azi = twopi * Urand();
    
    // convert to Cartesian coordinates
    double c = std::sqrt( 1.0 - mu * mu );
    point p;
    p.x = std::cos( azi ) * c;
    p.y = std::sin( azi ) * c;
    p.z = mu;
    
    return p;
}

point forwardPeakDirection_distribution::sample() {
    // sample polar cosine with forward peaked distribution
    // choose mu randomly from -1 to 1
    // choose pTry randomly from 0 to 1.5 (highest pTemp is 1.5)
    double mu, pTemp, pTry;
    do{
        mu      = 2.0 * Urand() - 1.0;
        pTemp   = (1.0/6.0)*(1.0+ pow((mu+1.0),3) );
        pTry    = Urand()*1.5;
    }
    while(pTry > pTemp);
    
    //sample azimuthal angle uniformly
    double azi = twopi * Urand();
    
    // convert to Cartesian coordinates
    double c = std::sqrt( 1.0 - mu * mu );
    point p;
    p.x = std::cos( azi ) * c;
    p.y = std::sin( azi ) * c;
    p.z = mu;
    
    return p;
}

point anisotropicDirection_distribution::sample() {
    double mu  = dist_mu->sample();
    double azi = twopi * Urand();
    double cos_azi = std::cos(azi);
    double sin_azi = std::sin(azi);
    
    // rotate the local particle coordinate system aligned along the incident direction
    // to the global problem (x,y,z) coordinate system
    double sin_t0 = std::sqrt( 1.0 - mu * mu );
    double c = sin_t0 / sin_t;
    
    point p;
    p.x = axis.x * mu + ( axis.x * axis.z * cos_azi - axis.y * sin_azi ) * c;
    p.y = axis.y * mu + ( axis.y * axis.z * cos_azi + axis.x * sin_azi ) * c;
    p.z = axis.z * mu - cos_azi * sin_t0 * sin_t;
    return p;
}

point independentXYZ_distribution::sample() {
    return point( dist_x->sample(), dist_y->sample(), dist_z->sample() );
}

point sampleThreePoint_distribution::sample(){
    double rand = Urand();
    if(rand < prob1){
        return point(  dist_1x->sample(), dist_1y->sample(), dist_1z->sample()  );
    }
    else if(rand >= prob1 && rand < (prob1+prob2)  ){
        return point(  dist_2x->sample(), dist_2y->sample(), dist_2z->sample()   );
    }
    else{
        return point(  dist_3x->sample(), dist_3y->sample(), dist_3z->sample()   );
    }
}

point sampleHollowSphere_distribution::sample(){
    
    double phi      = Urand() * twopi;
    double cosTheta = ( 2.0* Urand() ) - 1.0;
    double theta    = std::acos(cosTheta);
    double r;
    while (1) {
        r  = outRadius * pow(Urand(), (1.0/3.0) );
        if(r > inRadius && r < outRadius){
            break;
        }
    }
    
    double x = r * std::sin(theta) * std::cos(phi);
    double y = r * std::sin(theta) * std::sin(phi);
    double z = r * cosTheta;
    
    return point(x, y, z);
}

point diskPositionZ_distribution::sample(){
    
    double theta = Urand() * twopi;
    double r     = rad * std::sqrt(Urand());
    
    double x = oriX + r * std::cos(theta);
    double y = oriY + r * std::sin(theta);
    double z = oriZ + 0.0;
    
    return point(x, y, z);
}

point diskPositionX_distribution::sample(){
    
    double theta = Urand() * twopi;
    double r     = rad * std::sqrt(Urand());
    
    double x = oriX + 0.0;
    double y = oriY + r * std::sin(theta);
    double z = oriZ + r * std::cos(theta);
    
    return point(x, y, z);
}







