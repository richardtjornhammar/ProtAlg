#include <cstdlib>
#include <vector>
#include <string>
#include <cmath>
#include "mdiag.h"
#include "types.h"

//	CLIPPER	STUFF
#include <clipper/clipper.h>

//	GSL	STUFF
#include <gsl/gsl_matrix.h>


using namespace rich;

int
calc_map::proj(		gsl_matrix *P , gsl_matrix *C , clipper::Xmap<float> EDENS, 
			gsl_matrix *OS, gsl_vector *origo,
			double RC , double ZC, std::vector< std::pair<int,double> > * theta, clipper::Grid_sampling gs ) {

	std::vector<double> fi_cnt;
	for(int i=0;i<theta->size();i++){
		(*(theta))[i].first	= i;
		(*(theta))[i].second	= 0.0;
		fi_cnt.push_back(1E-10);
	}

	int verbose = 0;
	gsl_vector *rvec	= gsl_vector_calloc(DIM);
	gsl_vector *xh		= gsl_vector_calloc(DIM);
	gsl_vector *yh		= gsl_vector_calloc(DIM);
	gsl_vector *zh		= gsl_vector_calloc(DIM);

	gsl_matrix_get_row( zh , OS , 0 );
	gsl_matrix_get_row( xh , OS , 1 );
	gsl_matrix_get_row( yh , OS , 2 );

	if(P->size1!=P->size2 || P->size1 != nbins_ )
		return 1;
	if(OS->size1!=OS->size2 || OS->size1 != DIM )
		return 1;

	clipper::Coord_orth o0( gsl_vector_get(origo,XX),
				gsl_vector_get(origo,YY), 
				gsl_vector_get(origo,ZZ) );

	clipper::Xmap_base::Map_reference_index	ix;
	clipper::Coord_grid	cg0;
	clipper::Coord_orth	co0;
	double Ntheta  = theta->size();
	double x = 0,y = 0,z = 0;
	for ( ix = EDENS.first(); !ix.last(); ix.next() ) {
		double v0	= EDENS[ix];
		cg0		= ix.coord();
		co0		= ix.coord_orth();
		clipper::Coord_orth rtmp;
		rtmp = co0-o0;
		gsl_vector_set( rvec, XX , rtmp.x() );
		gsl_vector_set( rvec, YY , rtmp.y() );
		gsl_vector_set( rvec, ZZ , rtmp.z() );
		gsl_blas_ddot ( rvec, xh , &x );
		gsl_blas_ddot ( rvec, yh , &y );
		gsl_blas_ddot ( rvec, zh , &z );
		double r2	= (x*x+y*y);	// PROJECTED
		if( r2 < RC*RC && z*z<ZC*ZC ) { // && z>0 THIS THING
			double	res   	 = atan2(y,x)*180.0/M_PI+180.0;
			int	ires 	 = floor(res);
			if(verbose) {
				std::cout << "   " << x 
					  << "   " << y 
					  << "   " << z 
					  << " | " << v0 << std::endl;
				std::cout << ires  << " "<< v0 << std::endl;
			}

			int itheta = floor(res/360.0*Ntheta);
			if( itheta>=0 && itheta<Ntheta ){ // &&r2>RC*RC*0.25)
				(*(theta))[itheta].second	+= sqrt(v0*v0);
				fi_cnt[itheta]+=1.0;
 			}
			int I = floor( 0.5*(x/RC+1.0)*nbins_ );
			int J = floor( 0.5*(y/RC+1.0)*nbins_ );
			if(I>=0&&I<nbins_&&J>=0&&J<nbins_) {
				double pval=0.0,cval=0.0;
				pval	 = gsl_matrix_get(P,I,J);
				cval	 = gsl_matrix_get(C,I,J);
				pval	+= v0*v0;
				cval	+= 1.0;
				gsl_matrix_set( P, I, J, pval );
				gsl_matrix_set( C, I, J, cval );
			}
		}
	}
	for(int i=0;i<theta->size();i++){
		(*(theta))[i].second/=fi_cnt[i];
	}

	gsl_vector_free( rvec );
	gsl_vector_free(  xh  );
	gsl_vector_free(  yh  );
	gsl_vector_free(  zh  );

	return 0;
}


int
calc_map::proj01(	gsl_matrix *P , gsl_matrix *C , clipper::Xmap<float> EDENS, 
			gsl_matrix *OS, gsl_vector *origo,
			double RC , double ZC, std::vector< std::pair<int,double> > * theta, clipper::Grid_sampling gs ) {

	std::vector<double> fi_cnt;
	for(int i=0;i<theta->size();i++){
		(*(theta))[i].first	= i;
		(*(theta))[i].second	= 0.0;
		fi_cnt.push_back(1E-10);
	}
	int Nt = theta->size();
	int Nz = ceil( sqrt((Nt*Nt)/(M_PI*RC*RC))*ZC );

	int verbose = 0;
	gsl_vector *rvec	= gsl_vector_calloc(DIM);
	gsl_vector *xh		= gsl_vector_calloc(DIM);
	gsl_vector *yh		= gsl_vector_calloc(DIM);
	gsl_vector *zh		= gsl_vector_calloc(DIM);

	gsl_matrix_get_row( zh , OS , 0 );
	gsl_matrix_get_row( xh , OS , 1 );
	gsl_matrix_get_row( yh , OS , 2 );

	if(P->size1!=P->size2 || P->size1 != nbins_ )
		return 1;
	if(OS->size1!=OS->size2 || OS->size1 != DIM )
		return 1;

	clipper::Coord_orth o0( gsl_vector_get(origo,XX),
				gsl_vector_get(origo,YY), 
				gsl_vector_get(origo,ZZ) );

	clipper::Xmap_base::Map_reference_index	ix;
	clipper::Coord_grid	cg0;
	clipper::Coord_orth	co0;
	double Ntheta  = theta->size();
	double x = 0, y = 0, z = 0;
/*
	for( int i=0 ; i<N ; i++ ) {
		;
	}

	for ( ix = EDENS.first(); !ix.last(); ix.next() ) {
		double v0	= EDENS[ix];
		cg0		= ix.coord();
		co0		= ix.coord_orth();
		clipper::Coord_orth rtmp;
		rtmp = co0-o0;
		gsl_vector_set( rvec, XX , rtmp.x() );
		gsl_vector_set( rvec, YY , rtmp.y() );
		gsl_vector_set( rvec, ZZ , rtmp.z() );
		gsl_blas_ddot ( rvec, xh , &x );
		gsl_blas_ddot ( rvec, yh , &y );
		gsl_blas_ddot ( rvec, zh , &z );
		double r2	= (x*x+y*y);	// PROJECTED
		if( r2 < RC*RC && z*z<ZC*ZC ) { // && z>0 THIS THING
			double	res   	 = atan2(y,x)*180.0/M_PI+180.0;
			int	ires 	 = floor(res);
			if(verbose) {
				std::cout << "   " << x 
					  << "   " << y 
					  << "   " << z 
					  << " | " << v0 << std::endl;
				std::cout << ires  << " "<< v0 << std::endl;
			}

			int itheta = floor(res/360.0*Ntheta);
			if( itheta>=0 && itheta<Ntheta ){ // &&r2>RC*RC*0.25)
				(*(theta))[itheta].second	+= sqrt(v0*v0);
				fi_cnt[itheta]+=1.0;
 			}
			int I = floor( 0.5*(x/RC+1.0)*nbins_ );
			int J = floor( 0.5*(y/RC+1.0)*nbins_ );
			if(I>=0&&I<nbins_&&J>=0&&J<nbins_) {
				double pval=0.0,cval=0.0;
				pval	 = gsl_matrix_get(P,I,J);
				cval	 = gsl_matrix_get(C,I,J);
				pval	+= v0*v0;
				cval	+= 1.0;
				gsl_matrix_set( P, I, J, pval );
				gsl_matrix_set( C, I, J, cval );
			}
		}
	}
*/
	for(int i=0;i<theta->size();i++){
		(*(theta))[i].second/=fi_cnt[i];
	}

	gsl_vector_free( rvec );
	gsl_vector_free(  xh  );
	gsl_vector_free(  yh  );
	gsl_vector_free(  zh  );

	return 0;
}

// 	BELOW IS FAST BUT BOGUS FOR SOME REASON
int
calc_map::proj00(	gsl_matrix *P , gsl_matrix *C , clipper::Xmap<float> EDENS, 
			gsl_matrix *OS, gsl_vector *origo,
			double RC , double ZC, std::vector< std::pair<int,double> > * theta, clipper::Grid_sampling gs) {

	std::vector<double> fi_cnt;
	for(int i=0;i<theta->size();i++){
		(*(theta))[i].first	= i;
		(*(theta))[i].second	= 0.0;
		fi_cnt.push_back(1E-10);
	}

	int verbose = 0;
	gsl_vector *rvec	= gsl_vector_calloc(DIM);
	gsl_vector *xh		= gsl_vector_calloc(DIM);
	gsl_vector *yh		= gsl_vector_calloc(DIM);
	gsl_vector *zh		= gsl_vector_calloc(DIM);

	gsl_matrix_get_row( zh , OS , 0 );
	gsl_matrix_get_row( xh , OS , 1 );
	gsl_matrix_get_row( yh , OS , 2 );

	if( P->size1 !=P->size2  || P->size1 != nbins_ )
		return 1;
	if( OS->size1!=OS->size2 || OS->size1 != DIM   )
		return 1;

	clipper::Coord_orth o0( gsl_vector_get(origo,XX),
				gsl_vector_get(origo,YY), 
				gsl_vector_get(origo,ZZ) );

	clipper::Xmap_base::Map_reference_index	ix;
	//clipper::Grid		grid;

	clipper::Coord_grid	cg0 = (o0.coord_frac( EDENS.cell() )).coord_grid( gs );
	clipper::Coord_grid	c_x, c_d;
	double maxdist = 2*sqrt( RC*RC + ZC*ZC );
	clipper::Skeleton_basic::Neighbours neighb(EDENS, 0.25, maxdist );
	int 	n_neighbs	= neighb.size();
	float	f_neig		= neighb.size();

//	HERE SWEEP CORRECT BIN + NEIGHBORS INSTEAD OF ENTIRE DENSITY
	for (int i = 0 ; i<n_neighbs ; i++ ) {
		c_x	 	 = cg0 + neighb[i];
		c_d	 	 = neighb[i]; 			// ORIGO MIRROR AND SYMMETRY ACCOUNTED d(u,v,w) 
		double v0	 = EDENS.get_data(c_x);		// THE RIGHT PLACE
		clipper::Xmap_base::Map_reference_coord mrc( EDENS , c_d );
		clipper::Coord_orth rtmp( mrc.coord_orth() );
		double Ntheta  = theta->size();

		double x = 0, y = 0, z = 0;

		gsl_vector_set( rvec, XX , rtmp.x() );
		gsl_vector_set( rvec, YY , rtmp.y() );
		gsl_vector_set( rvec, ZZ , rtmp.z() );
		gsl_blas_ddot ( rvec, xh , &x );
		gsl_blas_ddot ( rvec, yh , &y );
		gsl_blas_ddot ( rvec, zh , &z );
		double r2	= (x*x+y*y); 		 	// PROJECTED
		if( 1 ) { // r2 < RC*RC && z*z < ZC*ZC ) { 	// && z>0 IF AWAY FROM BB
			double	res   	 = atan2(y,x)*180.0/M_PI+180.0;
			int	ires 	 = floor(res);

			if(verbose) {
				std::cout << "   " << x 
					  << "   " << y 
					  << "   " << z 
					  << " | " << v0 << std::endl;
				std::cout << ires  << " "<< v0 << std::endl;
			}

			int itheta = floor(res/360.0*Ntheta);
			if( itheta>=0 && itheta<Ntheta ) { // && r2>RC*RC*0.25 
				(*(theta))[itheta].second	+= sqrt(v0*v0);
				fi_cnt[itheta]+=1.0;
 			}
			int I = floor( 0.5*(x/RC+1.0)*nbins_ );
			int J = floor( 0.5*(y/RC+1.0)*nbins_ );
			if(I>=0 && I<nbins_ && J>=0 && J<nbins_) {
				double pval = 0.0, cval = 0.0;
				pval	 = gsl_matrix_get(P,I,J);
				cval	 = gsl_matrix_get(C,I,J);
				pval	+= v0*v0;
				cval	+= 1.0;
				gsl_matrix_set( P, I, J, pval );
				gsl_matrix_set( C, I, J, cval );
			}
		}
	}
	//ix.coord().coord_frac(EDENS.grid_sampling()).coord_orth(EDENS.cell())

	for(int i=0;i<theta->size();i++){
		(*(theta))[i].second/=fi_cnt[i];
	}

	gsl_vector_free( rvec );
	gsl_vector_free(  xh  );
	gsl_vector_free(  yh  );
	gsl_vector_free(  zh  );

	return 0;
}

