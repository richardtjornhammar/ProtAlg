#ifndef MDIAG_H
#define MDIAG_H

#include <cstdlib>
#include <vector>
#include <string>

//	CLIPPER	STUFF
#include <clipper/clipper.h>

//// GSL STUFF
#include <gsl/gsl_matrix.h>

namespace rich {
	class calc_map {
		public:
			inline calc_map( ) { nbins_=100; }
			int set_nbins( int N )	{ nbins_=N; return N; };
			int get_nbins( void  )	{ return nbins_; };
			int bmp (	gsl_matrix * , std::string );
			int proj(	gsl_matrix * , gsl_matrix * , clipper::Xmap<float> ,
					gsl_matrix * , gsl_vector * , double , double,
					std::vector< std::pair<int,double> > * );
		private:
			int nbins_;
	};
}

#endif
