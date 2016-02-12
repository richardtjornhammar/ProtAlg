#include <cstdlib>
#include <vector>
#include <string>
#include "mdiag.h"
#include "types.h"

//	CLIPPER	STUFF
#include <clipper/clipper.h>

using namespace rich;

int
calc_map::proj(	gsl_matrix *P  , clipper::Xmap<float> EDENS,
		gsl_matrix *OS , double RC , double ZC ) {

	if(P->size1!=P->size2 || P->size1 != nbins_ )
		return 1;
	if(OS->size1!=OS->size2 || OS->size1 != DIM )
		return 1;

	clipper::Xmap_base::Map_reference_index	ix;
	clipper::Coord_grid			cg0;

	for (ix = EDENS.first(); !ix.last(); ix.next()) {
		double v0	= EDENS[ix];
		cg0		= ix.coord();
	}

	return 0;
}
