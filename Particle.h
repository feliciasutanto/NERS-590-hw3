#ifndef _PARTICLE_HEADER_
#define _PARTICLE_HEADER_

#include <memory>

#include "Point.h"

class cell;

class particle {
private:
    point p_pos, p_dir;
    double p_energy;
    double p_speed;
    double p_wgt;
    double p_time;
    bool   exist;
    std::shared_ptr< cell > p_cell;
public:
    particle( point p, point d, double energy );
    ~particle() {};
    
    point   pos()      { return p_pos;       };       // return particle position
    point   dir()      { return p_dir;       };       // return particle direction
    double  wgt()      { return p_wgt;       };       // return particle weight
    double  time()     { return p_time;      };       // return distance travelled
    bool    alive()    { return exist;       };       // return particle state flag
    double  energy()   { return p_energy;    };       // return particle energy
    double  speed()    { return p_speed;     };       // return particle energy
    
    ray getRay() { return ray( p_pos, p_dir ); }
    
    std::shared_ptr< cell > cellPointer() { return p_cell; }
    
    void move( double s );
    void scatter( double mu0, double A );
    void kill();
    void setDirection( point p );
    void adjustWeight( double f );
    void recordCell( std::shared_ptr< cell > cel );
    void updateTime( double t);
    void setTime(double t);
    void setSpeed(double v);
};

#endif
