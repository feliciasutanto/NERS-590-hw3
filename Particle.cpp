#include <cmath>
#include <limits>
#include <iostream>
#include "Particle.h"
#include "Random.h"

// constructor for a new particle
particle::particle( point p, point d, double e ) : p_pos(p), p_dir(d), p_energy(e) {
    p_dir.normalize();
    exist      = true;
    p_wgt      = 1.0;
    p_time     = 0.0;
    p_cell     = nullptr;
    p_speed    = std::sqrt(2.0 * e / (1.0*931.493614838475)  ) * 29.9792458; //cm/ns
}

// move the particle along its current trajectory
void particle::move( double s ) {
    p_pos.x += s * p_dir.x;
    p_pos.y += s * p_dir.y;
    p_pos.z += s * p_dir.z;
}

// scatter the particle in center of mass given input direction cosine cos_t0 = mu0, and atomic mass of target A
void particle::scatter( double cos_t0, double A ) {
    
    point vlab(p_speed*p_dir.x, p_speed*p_dir.y, p_speed*p_dir.z); //cm/ns
    
    //calculate the u, velocity adjustment to go to center of mass frame
    point u(vlab.x / (1.0+A), vlab.y / (1.0+A), vlab.z / (1.0+A) ); //cm/ns
    
    //calculate velocity of neutron at center of mass frame
    point vc(vlab.x - u.x  ,  vlab.y - u.y ,  vlab.z - u.z); //cm.ns
    double speedc1 = std::sqrt( (vc.x*vc.x) + (vc.y*vc.y) + (vc.z*vc.z)   );
    
    //direction in center of mass
    point dirc(vc.x / speedc1, vc.y / speedc1, vc.z / speedc1); //no units
    
    // sample a random azimuthal angle uniformly
    double azi = 2.0 * std::acos(-1.0) * Urand();
    double cos_azi = std::cos(azi);
    double sin_azi = std::sin(azi);
    
    // rotate the local particle coordinate system aligned along the incident direction
    // to the global problem (x,y,z) coordinate system
    double sin_t  = std::sqrt( 1.0 - dirc.z  * dirc.z  );
    double sin_t0 = std::sqrt( 1.0 - cos_t0  * cos_t0 );
    
    point qc2; //no units
    if ( sin_t > std::numeric_limits<double>::epsilon() * 1000.0 ) {
        double c = sin_t0 / sin_t;
        qc2.x = dirc.x * cos_t0 + ( dirc.x * dirc.z * cos_azi - dirc.y * sin_azi ) * c;
        qc2.y = dirc.y * cos_t0 + ( dirc.y * dirc.z * cos_azi + dirc.x * sin_azi ) * c;
        qc2.z = dirc.z * cos_t0 - cos_azi * sin_t0 * sin_t;
    }
    else {
        // if incident direction along z, reorient axes to avoid division by zero
        sin_t  = std::sqrt( 1.0 -  dirc.y * dirc.y );
        double c = sin_t0 / sin_t;
        qc2.x = dirc.x * cos_t0 + ( dirc.x * dirc.y * cos_azi + dirc.z * sin_azi ) * c;
        qc2.y = dirc.y * cos_t0 - cos_azi * sin_t0 * sin_t;
        qc2.z = dirc.z * cos_t0 + ( dirc.y * dirc.z * cos_azi - dirc.x * sin_azi ) * c;
    }
    
    //calculate velocity of neutron at center of mass frame after collision
    point vc2(qc2.x * speedc1, qc2.y * speedc1, qc2.z * speedc1); //cm/ns
    
    //calculate velocity of neutron after coliision in lab frame
    point vlab2(vc2.x + u.x, vc2.y + u.y, vc2.z + u.z); //cm/ns
    double speedlab2 = std::sqrt( (vlab2.x*vlab2.x) + (vlab2.y*vlab2.y) + (vlab2.z*vlab2.z)  ); //cm/ns
    
    //update the new energy and speed
    double e2 = 0.5* (1.0*931.493614838475)  * pow((speedlab2*1.0/29.9792458), 2) ; // unit is MeV
    p_energy = e2;
    p_speed  = speedlab2;
    
    //update the new direction, no units
    p_dir.x = vlab2.x / speedlab2;
    p_dir.y = vlab2.y / speedlab2;
    p_dir.z = vlab2.z / speedlab2;
    
}

// set the particles life flag to false
void particle::kill() {
    exist = false;
}

// set the particle's direction
void particle::setDirection( point p ) {
    p_dir.x = p.x;
    p_dir.y = p.y;
    p_dir.z = p.z;
}

// adjust the weight by a factor f
void particle::adjustWeight( double f ) {
    p_wgt *= f;
}

// set the cell pointer for efficiency
void particle::recordCell( std::shared_ptr< cell > cel ) {
    p_cell = cel;
}

//update the time of the particle
void particle::updateTime( double t) {
    p_time += t ;
}

//set the time of the particle
void particle::setTime( double t) {
    p_time = t ;
}

//set the speed of the particle
void particle::setSpeed( double v) {
    p_speed = v ;
}




