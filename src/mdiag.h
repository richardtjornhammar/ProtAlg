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
	class map_diag {
		public:
			map_diag( ) { nbins_=100; }	
			int set_nBins( int N ) { nbins_=N; return N;} ;
			int make_bmp(	gsl_matrix * , std::string );
			int calc_prj(	gsl_matrix * , clipper::Xmap<float> ,
					gsl_matrix * , double , double  );
		private:
			int nbins_;
	};
}

#endif
