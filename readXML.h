#ifndef _READXML_HEADER_
#define _READXML_HEADER_

#include <vector>
#include <iostream>
#include <cstring>
#include <memory>
#include <stack>
#include <cmath>
#include <sstream>

#include "pugixml.hpp"
#include "Distribution.h"
#include "Reaction.h"
#include "Nuclide.h"
#include "Material.h"
#include "Surface.h"
#include "Cell.h"
#include "Source.h"

// Function that returns an item from a vector of objects of type T by name provided
// the object class has a string and a method called name() allowing for it to be returned
template< typename T >
std::shared_ptr< T > findByName( std::vector< std::shared_ptr< T > > vec, std::string name ) {
    for ( auto v : vec ) {
        if ( v->name() == name ) { return v; }
    }
    return nullptr;
}

void readXML(
             unsigned long long&                                      nsamples,
             double&                                                  tCut,
             std::vector < std::shared_ptr< surface   > >&            surfaces,
             std::vector < std::shared_ptr< cell      > >&            cells,
             std::vector < std::shared_ptr< nuclide   > >&            nuclides,
             std::vector < std::shared_ptr< material  > >&            materials,
             std::vector < std::shared_ptr< estimator > >&            estimators,
             std::shared_ptr< source >                   &            src,
             std::vector < std::shared_ptr<distribution<double>> >&   double_distributions,
             std::vector < std::shared_ptr<distribution<int>>    >&   int_distributions,
             std::vector < std::shared_ptr<distribution<point>>  >&   point_distributions
) {
    
    // user enters the XML file name and pugixml will attempt to load
    std::string input_file_name;
    std::cout << " Enter XML input file name: " << std::endl;
    std::cin  >> input_file_name;
    
    pugi::xml_document input_file;
    pugi::xml_parse_result load_result = input_file.load_file( input_file_name.c_str() );
    
    // check to see if result failed and throw an exception if it did
    if ( ! load_result ) {
        std::cout << load_result.description() << std::endl;
        throw;
    }
    
    // distributuions
    pugi::xml_node input_distributions = input_file.child("distributions");
    
    // find total number of distributions
    int num_distributions = 0;
    for ( auto d : input_distributions ) { num_distributions++; }
    
    // since distributions may depend on other distributions, need to iterate
    int set_distributions = 0;
    while ( set_distributions < num_distributions ) {
        int previous_set_distributions = set_distributions;
        
        for ( auto d : input_distributions ) {
            std::string type = d.name();
            std::string name = d.attribute("name").value();
            std::string data = d.attribute("datatype").value();
            
            if ( data == "double" ) {
                // skip rest of loop if distribution already done
                if ( findByName( double_distributions, name ) ) { continue; }
                
                std::shared_ptr< distribution<double> > Dist;
                if ( type == "delta" ) {
                    double a = d.attribute("a").as_double();
                    Dist = std::make_shared< arbitraryDelta_distribution< double > > ( name, a );
                }
                else if ( type == "uniform" ) {
                    double a = d.attribute("a").as_double();
                    double b = d.attribute("b").as_double();
                    Dist = std::make_shared< uniform_distribution > ( name, a, b );
                }
                else if ( type == "linAnisotropic" ){
                    double muAve = d.attribute("muBar").as_double();
                    Dist = std::make_shared< linAnisotropicDirection_distribution > ( name, muAve );
                }
                else if ( type == "linear" ) {
                    double a  = d.attribute("a").as_double();
                    double b  = d.attribute("b").as_double();
                    double fa = d.attribute("fa").as_double();
                    double fb = d.attribute("fb").as_double();
                    Dist = std::make_shared< linear_distribution > ( name, a, b, fa, fb );
                }
                else if ( type == "henyeyGreenstein" ) {
                    double a = d.attribute("a").as_double();
                    Dist = std::make_shared< HenyeyGreenstein_distribution > ( name, a );
                }
                else {
                    std::cout << "unsupported distribution with data type " << data << std::endl;
                    throw;
                }
                double_distributions.push_back( Dist );
            }
            // integer-valued distributions
            else if ( data == "int" ) {
                // skip rest of loop if distribution already done
                if ( findByName( int_distributions, name ) ) { continue; }
                
                std::shared_ptr< distribution<int> > Dist;
                if ( type == "delta" ) {
                    double a = d.attribute("a").as_int();
                    Dist = std::make_shared< arbitraryDelta_distribution< int > > ( name, a );
                }
                else if ( type == "meanMultiplicity" ) {
                    double nubar = d.attribute("nubar").as_double();
                    Dist = std::make_shared< meanMultiplicity_distribution > ( name, nubar );
                }
                else if ( type == "terrellFission" ) {
                    double nubar = d.attribute("nubar").as_double();
                    double sigma = d.attribute("sigma").as_double();
                    double b     = d.attribute("b").as_double();
                    Dist = std::make_shared< TerrellFission_distribution > ( name, nubar, sigma, b );
                }
                else {
                    std::cout << "unsupported distribution with data type " << data << std::endl;
                    throw;
                }
                int_distributions.push_back( Dist );
            }
            else if ( data == "point" ) {
                // skip rest of loop if distribution already done
                if ( findByName( point_distributions, name ) ) { continue; }
                
                std::shared_ptr< distribution< point > > Dist;
                if ( type == "delta" ) {
                    double x = d.attribute("x").as_double();
                    double y = d.attribute("y").as_double();
                    double z = d.attribute("z").as_double();
                    Dist = std::make_shared< arbitraryDelta_distribution< point > > ( name, point( x, y, z ) );
                }
                else if ( type == "isotropic" ) {
                    Dist = std::make_shared< isotropicDirection_distribution > ( name );
                }
                else if ( type == "forwardPeak" ) {
                    Dist = std::make_shared< forwardPeakDirection_distribution > ( name );
                }
                else if ( type == "anisotropic" ) {
                    double u = d.attribute("u").as_double();
                    double v = d.attribute("v").as_double();
                    double w = d.attribute("w").as_double();
                    std::shared_ptr< distribution<double> > angDist =
                    findByName( double_distributions, d.attribute("distribution").value() );
                    
                    // in the angular distribution does not yet, skip to the end of the loop
                    if ( ! angDist ) { continue; }
                    
                    Dist = std::make_shared< anisotropicDirection_distribution > ( name, point( u, v, w ), angDist );
                }
                else if ( type == "independentXYZ" ) {
                    std::shared_ptr< distribution<double> > distX = findByName( double_distributions, d.attribute("x").value() );
                    std::shared_ptr< distribution<double> > distY = findByName( double_distributions, d.attribute("y").value() );
                    std::shared_ptr< distribution<double> > distZ = findByName( double_distributions, d.attribute("z").value() );
                    
                    // if any of these distributions have not yet been resolved, skip to the end of the loop
                    if ( !distX || !distY || !distZ ) { continue; }
                    
                    Dist = std::make_shared< independentXYZ_distribution > ( name, distX, distY, distZ );
                }
                else if ( type == "sampleThreePoint" ) {
                    std::shared_ptr< distribution<double> > sourceFirstX  = findByName( double_distributions, d.attribute("source1x").value() );
                    std::shared_ptr< distribution<double> > sourceFirstY  = findByName( double_distributions, d.attribute("source1y").value() );
                    std::shared_ptr< distribution<double> > sourceFirstZ  = findByName( double_distributions, d.attribute("source1z").value() );
                    std::shared_ptr< distribution<double> > sourceSecondX = findByName( double_distributions, d.attribute("source2x").value() );
                    std::shared_ptr< distribution<double> > sourceSecondY = findByName( double_distributions, d.attribute("source2y").value() );
                    std::shared_ptr< distribution<double> > sourceSecondZ = findByName( double_distributions, d.attribute("source2z").value() );
                    std::shared_ptr< distribution<double> > sourceThirdX  = findByName( double_distributions, d.attribute("source3x").value() );
                    std::shared_ptr< distribution<double> > sourceThirdY  = findByName( double_distributions, d.attribute("source3y").value() );
                    std::shared_ptr< distribution<double> > sourceThirdZ  = findByName( double_distributions, d.attribute("source3z").value() );
                    double probFirst = d.attribute("prob1").as_double();
                    double probSecond= d.attribute("prob2").as_double();
                    double probThird = d.attribute("prob3").as_double();
                    
                    // if any of these distributions have not yet been resolved, skip to the end of the loop
                    if ( !sourceFirstX || !sourceSecondX || !sourceThirdX || !sourceFirstY || !sourceSecondY || !sourceThirdY || !sourceFirstZ || !sourceSecondZ || !sourceThirdZ ) { continue; }
                    
                    Dist = std::make_shared< sampleThreePoint_distribution > ( name, sourceFirstX, sourceFirstY, sourceFirstZ, sourceSecondX, sourceSecondY, sourceSecondZ, sourceThirdX, sourceThirdY, sourceThirdZ, probFirst, probSecond, probThird );
                }
                else if ( type == "sampleHollowSphere" ) {
                    double rad1 = d.attribute("r1").as_double();
                    double rad2 = d.attribute("r2").as_double();
                    Dist = std::make_shared< sampleHollowSphere_distribution > ( name, rad1, rad2);
                }
                else if ( type == "diskPositionZ" ) {
                    double xOri = d.attribute("x0").as_double();
                    double yOri = d.attribute("y0").as_double();
                    double zOri = d.attribute("z0").as_double();
                    double rad  = d.attribute("r").as_double();
                    Dist = std::make_shared< diskPositionZ_distribution > ( name, xOri, yOri, zOri, rad);
                }
                else if ( type == "diskPositionX" ) {
                    double xOri = d.attribute("x0").as_double();
                    double yOri = d.attribute("y0").as_double();
                    double zOri = d.attribute("z0").as_double();
                    double rad  = d.attribute("r").as_double();
                    Dist = std::make_shared< diskPositionZ_distribution > ( name, xOri, yOri, zOri, rad);
                }
                
                else {
                    std::cout << "unsupported " << data << " distribution of type " << type << std::endl;
                    throw;
                }
                point_distributions.push_back( Dist );
            }
            else {
                std::cout << "unsupported distribution with data type " << data << std::endl;
                throw;
            }
            // if we reach here, assume distribution has been set
            set_distributions++;
        }
        // check to see if number of distributions has increased, if not, caught in an infinite loop
        if ( previous_set_distributions == set_distributions ) {
            std::cout << "distributions could not be resolved. " << std::endl;
            throw;
        }
    }
    
    // iterate over nuclides
    pugi::xml_node input_nuclides = input_file.child("nuclides");
    for ( auto n : input_nuclides ) {
        std::string name = n.attribute("name").value();
        double      A    = n.attribute("A").as_double();
        
        std::shared_ptr< nuclide > Nuc = std::make_shared< nuclide > ( n.attribute("name").value(), A );
        nuclides.push_back( Nuc );
        
        // iterate over its reactions
        for ( auto r : n.children() ) {
            std::shared_ptr< reaction > Rxn;
            std::string rxn_type = r.name();
            
            double xsa = r.attribute("xsa").as_double();
            double xsb = r.attribute("xsb").as_double();
            
            if ( rxn_type == "capture" ) {
                Nuc->addReaction( std::make_shared< capture_reaction > ( xsa, xsb ) );
            }
            else if ( rxn_type == "scatter" ) {
                std::string dist_name = r.attribute("distribution").value();
                std::shared_ptr< distribution<double> > scatterDist = findByName( double_distributions, dist_name );
                if ( scatterDist ) {
                    Nuc->addReaction( std::make_shared< scatter_reaction > ( xsa, xsb, scatterDist ) );
                }
                else {
                    std::cout << " unknown scattering distribution " << dist_name << " in nuclide " << name << std::endl;
                    throw;
                }
            }
            else if ( rxn_type == "fission" ) {
                std::string mult_dist_name = r.attribute("multiplicity").value();
                std::shared_ptr< distribution<int> > multDist = findByName( int_distributions, mult_dist_name );
                if ( multDist ) {
                    Nuc->addReaction( std::make_shared< fission_reaction > ( xsa, xsb, multDist ) );
                }
                else {
                    std::cout << " unknown multiplicity distribution " << mult_dist_name << " in nuclide " << name << std::endl;
                    throw;
                }
            }
            else {
                std::cout << "unknown reaction type " << rxn_type << std::endl;
                throw;
            }
        }
    }
    
    // iterate over materials
    pugi::xml_node input_materials = input_file.child("materials");
    for ( auto m : input_materials ) {
        std::string name = m.attribute("name").value();
        double      aden = m.attribute("density").as_double();
        
        std::shared_ptr< material > Mat = std::make_shared< material > ( name, aden );
        materials.push_back( Mat );
        
        // iterate over nuclides
        for ( auto n : m.children() ) {
            if ( (std::string) n.name() == "nuclide" ) {
                std::string nuclide_name = n.attribute("name").value();
                double      frac         = n.attribute("frac").as_double();
                
                
                Mat->addNuclide( findByName( nuclides, nuclide_name ), frac );
            }
        }
    }
    
    // iterate over surfaces
    pugi::xml_node input_surfaces = input_file.child("surfaces");
    for ( auto s : input_surfaces ) {
        std::string type = s.name();
        
        std::shared_ptr< surface > S;
        if ( type == "plane" ) {
            std::string name = s.attribute("name").value();
            double      a    = s.attribute("a").as_double();
            double      b    = s.attribute("b").as_double();
            double      c    = s.attribute("c").as_double();
            double      d    = s.attribute("d").as_double();
            S = std::make_shared< plane > ( name, a, b, c, d );
        }
        else if ( type == "sphere" ) {
            std::string name = s.attribute("name").value();
            double      x0   = s.attribute("x0").as_double();
            double      y0   = s.attribute("y0").as_double();
            double      z0   = s.attribute("z0").as_double();
            double      rad  = s.attribute("rad").as_double();
            S = std::make_shared< sphere > ( name, x0, y0, z0, rad );
        }
        else if ( type == "cylinderx" ) {
            std::string name = s.attribute("name").value();
            double      y0   = s.attribute("y0").as_double();
            double      z0   = s.attribute("z0").as_double();
            double      rad  = s.attribute("rad").as_double();
            S = std::make_shared< cylinderx > ( name, y0, z0, rad);
        }
        else if ( type == "cylinderz" ) {
            std::string name = s.attribute("name").value();
            double      x0   = s.attribute("x0").as_double();
            double      y0   = s.attribute("y0").as_double();
            double      rad  = s.attribute("rad").as_double();
            S = std::make_shared< cylinderz > ( name, x0, y0, rad);
        }
        else if ( type == "conex" ) {
            std::string name = s.attribute("name").value();
            double      x0   = s.attribute("x0").as_double();
            double      y0   = s.attribute("y0").as_double();
            double      z0   = s.attribute("z0").as_double();
            double      a    = s.attribute("a").as_double();
            double      b    = s.attribute("b").as_double();
            double      c    = s.attribute("c").as_double();
            S = std::make_shared< conex > ( name, x0, y0, z0, a, b, c);
        }
        else {
            std::cout << " unkown surface type " << type << std::endl;
            throw;
        }
        
        if ( (std::string) s.attribute("bc").value() == "reflect" ) {
            S->makeReflecting();
        }
        surfaces.push_back( S );
    }
    
    // iterate over cells
    pugi::xml_node input_cells = input_file.child("cells");
    for ( auto c : input_cells ) {
        std::string name = c.attribute("name").value();
        
        std::shared_ptr< cell > Cel = std::make_shared< cell > ( name );
        cells.push_back( Cel );
        
        // cell material
        if ( c.attribute("material") ) {
            std::shared_ptr< material > matPtr = findByName( materials, c.attribute("material").value() );
            if ( matPtr ) {
                Cel->setMaterial( matPtr );
            }
            else {
                std::cout << " unknown material " << c.attribute("material").value() << " in cell " << name << std::endl;
                throw;
            }
        }
        
        // cell importance
        if ( c.attribute("importance") ) {
            Cel->setImportance( c.attribute("importance").as_double() );
        }
        
        // iterate over surfaces
        for ( auto s : c.children() ) {
            if ( (std::string) s.name() == "surface" ) {
                std::string name  = s.attribute("name").value();
                int         sense = s.attribute("sense").as_int();
                
                std::shared_ptr< surface > SurfPtr = findByName( surfaces, name );
                if ( SurfPtr ) {
                    Cel->addSurface( findByName( surfaces, name ), sense );
                }
                else {
                    std::cout << " unknown surface with name " << name << std::endl;
                    throw;
                }
            }
            else {
                std::cout << " unknown data type " << s.name() << " in cell " << name << std::endl;
                throw;
            }
        }
    }
    
    // iterate over estimatators
    pugi::xml_node input_estimators = input_file.child("estimators");
    for ( auto e : input_estimators ) {
        std::string type = e.name();
        std::string name = e.attribute("name").value();
        
        std::shared_ptr< estimator > Est;
        if ( type == "current" ) {
            
            pugi::xml_node aaa = input_estimators.child("current");
            pugi::xml_node bbb = aaa.child("surface");
            pugi::xml_node ccc = aaa.child("lowBin");
            pugi::xml_node ddd = aaa.child("dBin");
            pugi::xml_node eee = aaa.child("numBin");
            
            std::string nameSurface = bbb.attribute("surface").value();
            double      lowBin      = ccc.attribute("lowBin").as_double();
            double      dBin        = ddd.attribute("dBin").as_double();
            int         numBin      = eee.attribute("numBin").as_int();
            
            
            Est = std::make_shared< surface_current_estimator > ( name, numBin, lowBin, dBin, numBin, nameSurface );
            
            // get the surfaces
            for ( auto s : e.children() ) {
                if ( (std::string) s.name() == "surface" ) {
                    std::string name = s.attribute("name").value();
                    std::shared_ptr< surface > SurfPtr = findByName( surfaces, name );
                    if ( SurfPtr ) {
                        SurfPtr->attachEstimator( Est );
                    }
                    else {
                        std::cout << " unknown surface label " << name << " in estimator " << e.attribute("name").value() << std::endl;
                    }
                }
            }
        }
        else if ( type == "countingSurface" ) {
            
            pugi::xml_node aaa = input_estimators.child("countingSurface");
            pugi::xml_node bbb = aaa.child("surface");
            pugi::xml_node ccc = aaa.child("lowBin");
            pugi::xml_node ddd = aaa.child("dBin");
            pugi::xml_node eee = aaa.child("numBin");
            
            std::string nameSurface = bbb.attribute("surface").value();
            double      lowBin      = ccc.attribute("lowBin").as_double();
            double      dBin        = ddd.attribute("dBin").as_double();
            int         numBin      = eee.attribute("numBin").as_int();
            
            Est = std::make_shared< counting_estimator > ( name, nameSurface ); //, lowBin, dBin, numBin );
            
            // get the surfaces
            for ( auto s : e.children() ) {
                if ( (std::string) s.name() == "surface" ) {
                    std::string name = s.attribute("name").value();
                    std::shared_ptr< surface > SurfPtr = findByName( surfaces, name );
                    if ( SurfPtr ) {
                        SurfPtr->attachEstimator( Est );
                    }
                    else {
                        std::cout << " unknown surface label " << name << " in estimator " << e.attribute("name").value() << std::endl;
                    }
                }
            }
        }
        else if( type == "pathLengthFlux"  ){
            
            pugi::xml_node aaa = input_estimators.child("pathLengthFlux");
            pugi::xml_node bbb = aaa.child("volume");
            pugi::xml_node ccc = aaa.child("lowBin");
            pugi::xml_node ddd = aaa.child("dBin");
            pugi::xml_node eee = aaa.child("numBin");
            pugi::xml_node fff = aaa.child("type");
            pugi::xml_node ggg = aaa.child("cell");
            
            double      volume      = bbb.attribute("vol").as_double();
            double      lowBin      = ccc.attribute("lowBin").as_double();
            double      dBin        = ddd.attribute("dBin").as_double();
            int         numBin      = eee.attribute("numBin").as_int();
            std::string type        = fff.attribute("type").value();
            std::string nameCell    = ggg.attribute("cell").value();
            
            Est = std::make_shared< cell_pathLengthFlux_estimator > ( name , numBin, volume, lowBin, dBin, type, nameCell );
            
            // get the surfaces
            for ( auto c : e.children() ) {
                if ( (std::string) c.name() == "cell" ) {
                    std::string name = c.attribute("name").value();
                    std::shared_ptr< cell > CellPtr = findByName( cells, name );
                    if ( CellPtr ) {
                        CellPtr->attachEstimator( Est );
                    }
                    else {
                        std::cout << " unknown cell label " << name << " in estimator " << e.attribute("name").value() << std::endl;
                    }
                }
            }
        }
        
        else {
            std::cout << "unknown estimator type " << name << std::endl;
            throw;
        }
        
        estimators.push_back( Est );
    }
    
    // create source
    pugi::xml_node input_source = input_file.child("source");
    pugi::xml_node input_source_position  = input_source.child("position");
    pugi::xml_node input_source_direction = input_source.child("direction");
    pugi::xml_node input_source_energy    = input_source.child("energy");
    
    std::string pos_dist_name    = input_source_position.attribute("distribution").value();
    std::string dir_dist_name    = input_source_direction.attribute("distribution").value();
    std::string erg_dist_name    = input_source_energy.attribute("distribution").value();
    
    std::shared_ptr< distribution< point > > posDist     = findByName( point_distributions, pos_dist_name );
    std::shared_ptr< distribution< point > > dirDist     = findByName( point_distributions, dir_dist_name );
    std::shared_ptr< distribution< double > >ergDist     = findByName( double_distributions,erg_dist_name );
    
    //std::shared_ptr< source > src;
    if ( posDist && dirDist && ergDist ) {
        src = std::make_shared< source > ( posDist, dirDist, ergDist );
    }
    else {
        if ( ! posDist )  { std::cout << " unknown position distribution "  << pos_dist_name   << " in source " << std::endl; }
        if ( ! dirDist )  { std::cout << " unknown direction distribution " << dir_dist_name   << " in source " << std::endl; }
        if ( ! ergDist )  { std::cout << " unknown energy distribution "    << erg_dist_name   << " in source " << std::endl; }
        throw;
    }
    
    // read number of histories
    pugi::xml_node input_numHistory          = input_file.child("simulation");
    pugi::xml_node input_numHistory_nsamples1 = input_numHistory.child("histories");
    pugi::xml_node input_numHistory_nsamples2 = input_numHistory.child("timeCut");
    nsamples = input_numHistory_nsamples1.attribute("end").as_int();
    tCut     = input_numHistory_nsamples2.attribute("tCut").as_double();
}


#endif
