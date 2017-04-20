#include <cmath>
#include <limits>
#include <cassert>

#include "Point.h"
#include "QuadSolver.h"
#include "Surface.h"

//PLANE----------------------------------------------------------------------------
// evaluates the surface equation w.r.t. to point p
double plane::eval( point p ) {
  return a * p.x  +  b * p.y  +  c * p.z  - d;
}

// determines the mininum positive distance to intersection for a ray r
// (returns a very large number if no intersection along ray for ease of calculation down the line)
double plane::distance( ray r ) {
  point p = r.pos;
  point u = r.dir;

  double denom = a * u.x  +  b * u.y  +  c * u.z;
  if ( std::fabs( denom ) > 100.0 * std::numeric_limits<double>::epsilon() ) {
    double dist = ( d - a * p.x - b * p.y - c * p.z ) / denom;
    if ( dist > 0.0 ) { return dist; }
    else { return std::numeric_limits<double>::max(); }
  }
  else {
    // moving in a direction that is (or is very close to) parallel to the surface
    return std::numeric_limits<double>::max();
  }
}

// get new reflected direction
point plane::reflect( ray r ) {
   assert( std::fabs( eval( r.pos ) ) < std::numeric_limits<float>::epsilon() );

   point  u = r.dir;
   double t = 2.0 * ( a * u.x  +  b * u.y  +  c * u.z ) / ( a*a + b*b + c*c );
   return point( u.x - a*t, u.y - b*t, u.z - c*t );
}

//SPHERE----------------------------------------------------------------------------
double sphere::eval( point p ) {
  return std::pow( p.x - x0, 2 ) + std::pow( p.y - y0, 2 ) + std::pow( p.z - z0, 2 )  - rad*rad;
}

double sphere::distance( ray r ) {
  point p = r.pos;
  point u = r.dir;

  // difference between each coordinate and current point
  point q( p.x - x0, p.y - y0, p.z - z0 );

  // put into quadratic equation form: a*s^2 + b*s + c = 0, where a = 1
  double b = 2.0 * ( q.x * u.x  +  q.y * u.y  +  q.z * u.z );
  double c = eval( p );

  return quad_solve( 1.0, b, c );
}

point sphere::reflect( ray r ) {
   assert( std::fabs( eval( r.pos ) ) < std::numeric_limits<float>::epsilon() );

   point p = r.pos;
   point u = r.dir;

   point q( p.x - x0, p.y - y0, p.z - z0 );

   double t = 2.0 * ( q.x * u.x  +  q.y * u.y  +  q.z * u.z ) / ( rad*rad );
   return point( u.x - q.x * t,  u.y - q.y * t,  u.z - q.z * t );
}

//CYLINDER X--------------------------------------------------------------------------
double cylinderx::eval( point p ) {
    return std::pow( p.y - y0, 2 ) + std::pow( p.z - z0, 2 )  - rad*rad;
}

double cylinderx::distance( ray r ) {
    point p = r.pos;
    point u = r.dir;
    point q( 0.0    , p.y - y0 , p.z - z0 );
    
    // put into quadratic equation form: a*s^2 + b*s + c = 0
    double a = u.y * u.y + u.z * u.z;
    double b = 2.0 * ( q.y * u.y  +  q.z * u.z );
    double c = q.y * q.y + q.z * q.z - rad * rad;
    
    return quad_solve( a, b, c );
}

point cylinderx::reflect( ray r ) {
    assert( std::fabs( eval( r.pos ) ) < std::numeric_limits<float>::epsilon() );
    
    point p = r.pos;
    point u = r.dir;
    
    point q( 0.0   , p.y - y0, p.z - z0 );
    
    double t = 2.0 * ( q.x * u.x  + q.y * u.y  +  q.z * u.z ) / ( rad*rad );
    return point( u.x - q.x * t ,  u.y - q.y * t,  u.z - q.z * t );

}

//CYLINDER Z--------------------------------------------------------------------------
double cylinderz::eval( point p ) {
    return std::pow( p.x - x0, 2 ) + std::pow( p.y - y0, 2 )  - rad*rad;
}

double cylinderz::distance( ray r ) {
    point p = r.pos;
    point u = r.dir;
    point q( p.x - x0  , p.y - y0 , 0.0 );
    
    // put into quadratic equation form: a*s^2 + b*s + c = 0
    double a = u.y * u.y + u.x * u.x;
    double b = 2.0 * ( q.y * u.y  +  q.x * u.x );
    double c = q.y * q.y + q.x * q.x - rad * rad;
    
    return quad_solve( a, b, c );
}

point cylinderz::reflect( ray r ) {
    assert( std::fabs( eval( r.pos ) ) < std::numeric_limits<float>::epsilon() );
    
    point p = r.pos;
    point u = r.dir;
    
    //normal vector
    point q( p.x - x0 , p.y - y0, 0.0 );
    
    double t = 2.0 * ( q.x * u.x  + q.y * u.y  +  q.z * u.z ) / ( rad*rad );
    return point( u.x - q.x * t ,  u.y - q.y * t,  u.z - q.z * t );
    
}

//CONE X-----------------------------------------------------------------------------
double conex::eval( point p ) {
    
    a2  = aDenom*aDenom;
    b2  = bDenom*bDenom;
    c2  = cDenom*cDenom;
    
    return std::pow( p.y - y0, 2 )/b2 + std::pow( p.z - z0, 2 )/c2   -  std::pow( p.x - x0, 2 )/a2 ;
}

double conex::distance( ray r ) {
    
    a2  = aDenom*aDenom;
    b2  = bDenom*bDenom;
    c2  = cDenom*cDenom;
    
    point p = r.pos;
    point u = r.dir;
    point q( p.x - x0  , p.y - y0 , p.z - z0 );
    
    // put into quadratic equation form: a*s^2 + b*s + c = 0
    double a = ( u.y * u.y / b2 ) + ( u.z * u.z / c2 ) - ( u.x * u.x /a2 );
    double b = 2.0 * ( (q.y * u.y /b2) + (q.z * u.z /c2) -  (q.x * u.x /a2) );
    double c = (q.y * q.y /b2) + (q.z * q.z /c2) - (q.x * q.x /a2);
    
    return quad_solve( a, b, c );
}

point conex::reflect( ray r ) { //need to check this, but it doesn't matter for this hw
    assert( std::fabs( eval( r.pos ) ) < std::numeric_limits<float>::epsilon() );
    
    a2  = aDenom*aDenom;
    b2  = bDenom*bDenom;
    c2  = cDenom*cDenom;
    
    point p = r.pos;
    point u = r.dir;
    
    //normal vector
    point q( -(p.x - x0)/a2 , (p.y - y0)/b2 , (p.z - z0)/c2 );
    
    double t = 2.0 * ( q.x * u.x  + q.y * u.y  +  q.z * u.z );
    return point( u.x - q.x * t ,  u.y - q.y * t,  u.z - q.z * t );
    
}






