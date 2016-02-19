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


double
math_helper::gsl_calc_orth(	gsl_vector *n2, gsl_vector *n1, gsl_vector *c2,
				gsl_vector *c1, gsl_matrix *OS ) { // should rewrite in tensor notation

	if( n1->size!=n2->size  || c1->size!=c2->size || OS->size1!=DIM && OS->size2!=DIM )
		std::cout << "INFO:: CANNOT COMPUTE CROSS PRODUCTS" << std::endl;

	gsl_vector *nh = gsl_vector_calloc(DIM);
	gsl_vector *ph = gsl_vector_calloc(DIM);
	gsl_vector *qh = gsl_vector_calloc(DIM);

	gsl_vector_sub(n2,n1);
	gsl_vector_sub(c2,c1);
	gsl_vector_memcpy(nh,n2);
	gsl_vector_memcpy(ph,c2);
	double nl = gsl_blas_dnrm2(nh);
	double pl = gsl_blas_dnrm2(ph);
	gsl_vector_scale(nh,1.0/nl); gsl_vector_scale(ph,1.0/pl);

	gsl_cross3D( nh, ph, qh );

	double ql = gsl_blas_dnrm2(qh);			
	gsl_vector_scale(qh, 1.0/ql );

	gsl_cross3D( nh, qh, ph );

	pl = gsl_blas_dnrm2(ph);
	gsl_vector_scale(ph, -1.0/pl);

	gsl_matrix_set_row(OS,0,nh);
	gsl_matrix_set_row(OS,1,ph);
	gsl_matrix_set_row(OS,2,qh);

	gsl_vector_free(nh);
	gsl_vector_free(ph);
	gsl_vector_free(qh);

	return nl;
}

int
mmdb_helper::check_clash( int imod, int icha, int iRes, CMMDBManager *MMDB, particles residue, double clashdist ) {

	PPCAtom		SelAtom;
	int		RC,selHnd,nSelAtoms,nAtoms;
 	int		nRes = MMDB->GetNumberOfResidues(imod,icha);

	nAtoms = residue.size();
 	selHnd = MMDB->NewSelection();
	for ( int k=0 ; k<nAtoms ; k++ ) {
		double x = gsl_vector_get( residue[k].second, XX );
		double y = gsl_vector_get( residue[k].second, YY );
		double z = gsl_vector_get( residue[k].second, ZZ );
		MMDB->Select (
			selHnd,
			STYPE_ATOM,
			imod,			// MODEL
			"*",			// CHAINS
     			ANY_RES	,"*",		// RESNR,
     			ANY_RES ,"*",		// INSCODE
     			"*",			// RESNAME
     			"*",			// ATOMNAME
     			"*",			// ATOMTYPES
     			"*",			// LOCINDICATOR
     			"*",			// SEGMENT
     			"*",			// CHARGES
     			-1.0,-1.0,		// OCCUPANCY (ANY)
     			x, y, z, clashdist,	// SPHERE
     			SKEY_OR      		// OR-selection
		);
	}
	MMDB->GetSelIndex ( selHnd, SelAtom, nSelAtoms );
	
	return nSelAtoms;
}


int
mmdb_helper::update_residue( int imod, int icha, int ires, CMMDBManager *mmdb0, particles residue ) {
	PPCAtom	atoms;
	int 	nAtoms;

	mmdb0->GetAtomTable    ( imod ,icha ,ires , atoms, nAtoms );

	if( nAtoms == residue.size() ) {
		for (int k=0;k<nAtoms;k++){
			atoms[k]->x=gsl_vector_get(residue[k].second,XX);
			atoms[k]->y=gsl_vector_get(residue[k].second,YY);
			atoms[k]->z=gsl_vector_get(residue[k].second,ZZ);
		}
		return 0;
	}else{
		return 1;
	}
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

int	
quaternion::assign_quaterion( gsl_vector *v, double angle ){ 
	// angle in radians
	if( v->size==DIM ){

		double norm=0.0;

		double vx = gsl_vector_get( v, XX );
		double vy = gsl_vector_get( v, YY );
		double vz = gsl_vector_get( v, ZZ );

  		norm   = 1.0/sqrt(vx*vx+vy*vy+vz*vz);

		double w = cos(angle*0.5);
		double x = vx*norm*sin(angle*0.5);
		double y = vy*norm*sin(angle*0.5);
		double z = vz*norm*sin(angle*0.5);

		gsl_vector_set ( q_ , XX , w );
		gsl_vector_set ( q_ , YY , x );
		gsl_vector_set ( q_ , ZZ , y );
		gsl_vector_set ( q_ ,DIM , z );

		gsl_matrix_set( R_, XX, XX, 1.0-2.0*(y*y+z*z)	);
		gsl_matrix_set( R_, XX, YY, 2*(x*y-z*w)		);
		gsl_matrix_set( R_, XX, ZZ, 2*(x*z+y*w)		);

		gsl_matrix_set( R_, YY, XX, 2*(x*y+z*w)		);
		gsl_matrix_set( R_, YY, YY, 1-2*(x*x+z*z)	);
		gsl_matrix_set( R_, YY, ZZ, 2*(y*z-x*w)		);

		gsl_matrix_set( R_, ZZ, XX, 2*(x*z-y*w)		);
		gsl_matrix_set( R_, ZZ, YY, 2*(y*z+x*w)		);
		gsl_matrix_set( R_, ZZ, ZZ, 1-2*(x*x+y*y)	);

		// R=[1-2*(y*y+z*z) 2*(x*y-z*w) 2*(x*z+y*w) ; 2*(x*y+z*w) 1-2*(x*x+z*z) 2*(y*z-x*w) ; 2*(x*z-y*w) 2*(y*z+x*w) 1-2*(x*x+y*y) ]

		bSet_=true;

		return 0;
	}else {
		return 1;
	}
}

int
quaternion::rotate_particles( particles ps ) {
	if( is_complete() ){
		for(int i=0;i<ps.size();i++){
			rotate_coord(ps[i].second);
		}
	}
	else {
		return 1;
	}

	return 0;
}

void
quaternion::print(){
	if( is_complete() ){
		tensorIO tIO;
		std::cout << "INFO::QUATERNION" << std::endl;
		tIO.output_vector( q_ );
		tIO.output_matrix( R_ );
	}
}

int
quaternion::rotate_particles( particles ps , gsl_vector *v0 ) {
	if( is_complete() ){
		gsl_vector *tmp = gsl_vector_calloc(DIM);
		for(int i=0;i<ps.size();i++){
			gsl_vector_sub( ps[i].second, v0 );
			rotate_coord  ( ps[i].second	 );
			gsl_vector_add( ps[i].second, v0 );
		}
		gsl_vector_free(tmp);
	}
	else {
		return 1;
	}

	return 0;
}

int
quaternion::rotate_coord( gsl_vector *x ) 
{
	double xX,yY,zZ;

	if( is_complete() && x->size==DIM ){

		if(0){
			double q[DIM+1],xo[DIM];
			q[XX] = gsl_vector_get(q_,XX); q[YY ] = gsl_vector_get(q_,YY );
			q[ZZ] = gsl_vector_get(q_,ZZ); q[DIM] = gsl_vector_get(q_,DIM);

			xo[XX] = gsl_vector_get(x, XX); 
			xo[YY] = gsl_vector_get(x, YY); 
			xo[ZZ] = gsl_vector_get(x, ZZ);
 
			xX = 	(q[0]*q[0] + q[1]*q[1]-q[2]*q[2]-q[3]*q[3] )	* xo[XX] + 
				(2*q[1]*q[2] - 2*q[0]*q[3])			* xo[YY] + 
				(2*q[1]*q[3] + 2*q[0]*q[2])			* xo[ZZ] ;

			yY = 	(2*q[1]*q[2] + 2*q[0]*q[3])			* xo[XX] + 
				(q[0]*q[0]-q[1]*q[1] + q[2]*q[2]-q[3]*q[3] )	* xo[YY] + 
				(2*q[2]*q[3]-2*q[0]*q[1])			* xo[ZZ] ;

			zZ = 	(2*q[1]*q[3] - 2*q[0]*q[2])			* xo[XX] + 
			 	(2*q[2]*q[3] + 2*q[0]*q[1])			* xo[YY] + 
			 	(q[0]*q[0]-q[1]*q[1]-q[2]*q[2] + q[3]*q[3] )	* xo[ZZ] ;

			gsl_vector_set(x,XX,xX);
			gsl_vector_set(x,YY,yY);
			gsl_vector_set(x,ZZ,zZ);
		}else{
			gsl_blas_dgemv (CblasNoTrans ,1.0, R_, x, 0.0, x);
		}

  		return 0;
	}else{
		return 1;
	}
}

