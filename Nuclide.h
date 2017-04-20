#ifndef _NUCLIDE_HEADER_
#define _NUCLIDE_HEADER_

#include <vector>
#include <string>
#include <memory>

#include "Reaction.h"

class nuclide {
private:
    std::string nuclide_name;
    std::vector< std::shared_ptr< reaction > > rxn;
    double A;
public:
    nuclide( std::string label, double a ) : nuclide_name(label), A(a) {};
    ~nuclide() {};
    
    std::string name() { return nuclide_name; }
    std::vector< std::shared_ptr< reaction > > getReactions() {return rxn;} ;
    
    void        addReaction( std::shared_ptr< reaction > );
    double      total_xs( particle* p );
    double      cap_xs( particle* p );
    double      scat_xs( particle* p );
    double      getA(){ return A; };
    
    std::shared_ptr< reaction > sample_reaction( particle* p );
};


#endif
