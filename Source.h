#ifndef _SOURCE_HEADER_
#define _SOURCE_HEADER_

#include <stack>
#include <memory>

#include "Point.h"
#include "Distribution.h"
#include "Particle.h"

class source {
private:
    std::shared_ptr< distribution<point> > dist_pos;
    std::shared_ptr< distribution<point> > dist_dir;
    std::shared_ptr< distribution<double> > dist_energy;
public:
    source( std::shared_ptr< distribution<point> > pos, std::shared_ptr< distribution<point> > dir, std::shared_ptr< distribution<double> > erg )
    : dist_pos(pos), dist_dir(dir), dist_energy(erg) {};
    ~source() {};
    
    std::stack<particle> sample();
};

#endif
