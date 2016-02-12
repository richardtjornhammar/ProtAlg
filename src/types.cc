#include <cstdlib>
#include <vector>
#include <string>
#include <math.h>
#include "iofun.h"

using namespace rich;

void
math_helper::gsl_cross3D(gsl_vector *a, gsl_vector *b, gsl_vector *c) { // should rewrite in tensor notation
	double da[3];

	if( a->size!=b->size  || a->size!=c->size || a->size!=3)
		std::cout << "INFO:: CANNOT COMPUTE CROSS PRODUCT" << std::endl;

	da[0] = gsl_vector_get(a,1)*gsl_vector_get(b,2) - gsl_vector_get(a,2)*gsl_vector_get(b,1);
	da[1] = gsl_vector_get(a,2)*gsl_vector_get(b,0) - gsl_vector_get(a,0)*gsl_vector_get(b,2);
	da[2] = gsl_vector_get(a,0)*gsl_vector_get(b,1) - gsl_vector_get(a,1)*gsl_vector_get(b,0);

	gsl_vector_set(c,0,da[0]);
	gsl_vector_set(c,1,da[1]);
	gsl_vector_set(c,2,da[2]);
}

std::string
mmdb_helper::atom_type( char * aid ) {
	std::string atom_str(aid);
	std::size_t found0	= atom_str.find_last_of("/\\");
	std::size_t found1	= atom_str.find_first_of("[");
	std::string atype	= atom_str.substr(found0+1,found1-1-found0);
	return atype;
}

std::string
mmdb_helper::atom_symb( char * aid ) {
	std::string atom_str(aid);
	std::size_t found0	= atom_str.find_last_of ("[");
	std::size_t found1	= atom_str.find_first_of("]");
	std::string atype	= atom_str.substr(found0+1,found1-1-found0);
	return atype;
}

////fft_reso = clipper::Resolution(1.0/sqrt(fphi_.invresolsq_range().max()));
int
map_manager::assign_map( std::string inp_str ) {

	int verbose=0;
	map_io mio;
	mio.read_map_header( inp_str );

	clipper::String f_str=mio.get_ifil();
	clipper::String c_str=mio.get_isel();
	clipper::String t_str=mio.get_itot();

	clipper::CCP4MTZfile mtzin;
	if(verbose)
		std::cout << "INFO:: OPEN READ MTZ" << std::endl;
	mtzin.open_read( f_str );
	if(verbose)
		std::cout << "INFO:: DONE \nINFO:: OPEN READ CRYSTAL" << std::endl;
	mtzin.import_crystal ( xtl_ , t_str );
	if(verbose)
		std::cout << "INFO:: DONE \nINFO:: OPEN READ HKL" << std::endl;
  	mtzin.import_hkl_data( fphi_, set_, xtl_ , t_str );
	if(verbose)
		std::cout << "INFO:: DONE \nINFO:: CLOSE MTZ" << std::endl;
	mtzin.close_read();
	if(verbose)
		std::cout << "INFO:: DONE IO\nINFO:: SET GRID SAMPLING" << std::endl;
	m_s_rate_=1.5;
	clipper::Grid_sampling gs( fphi_.spacegroup(), fphi_.cell(), fphi_.resolution(), m_s_rate_ );
	xmap_.init	( fphi_.spacegroup(), fphi_.cell(), gs ); 		// INITIALIZE MAP
	if(verbose)
		std::cout << "INFO:: DONE GS\nINFO:: CALC FFT" << std::endl;
	xmap_.fft_from	( fphi_ );
	if(verbose)
		std::cout << "INFO:: DONE FFT" << std::endl;

	return 0;	
}
