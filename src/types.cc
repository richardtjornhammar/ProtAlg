#include <cstdlib>
#include <vector>
#include <string>
#include <math.h>
#include <algorithm>
#include "iofun.h"
#include "mdiag.h"

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
mmdb_helper::copy_xtal( CMMDBManager* mmdb, CMMDBManager* mmdb_N ) {

	int RC;

	if (mmdb->isSpaceGroup() && !(mmdb_N->isSpaceGroup()) )  {
		// space group name is for demonstration only
		RC = mmdb_N->SetSpaceGroup ( mmdb->GetSpaceGroup () );
		if (RC!=SYMOP_Ok)   {
			switch (RC)  {
				case SYMOP_NoLibFile :
					printf ( " **** error: can't find symop.lib\n" );
					break;
				case SYMOP_UnknownSpaceGroup :
					printf ( " **** error: attempt to set up unknown space group\n" );
					break;
				case SYMOP_NoSymOps :
					printf ( " **** error: no symmetry operations found\n" );
					break;
				default :
					printf ( " **** error: unknown return code from "
					"CMMDBManager::SetSpaceGroup()\n" );
				}
			exit(2);
		}
		std::cout << "INFO:: ASSIGNED SYM. INFORMATION"<< std::endl;
	}

	if ( mmdb -> isCrystInfo() )  {
		double cell[8];
		int OC[0];
    		mmdb->GetCell ( cell[0], cell[1], cell[2], cell[3] , cell[4] , cell[5] , cell[6] , OC[0] );
		mmdb_N->SetCell ( cell[0], cell[1], cell[2], cell[3] , cell[4] , cell[5] , OC[0] );
		std::cout << "INFO:: ASSIGNED CELL INFORMATION"<< std::endl;
	}
	RC = mmdb_N->CrystReady();

	if ( RC>0 )  {
		if (RC & CRRDY_NotPrecise)
			std::cout << "WARNING:: 1 \n" ;
		if (RC & CRRDY_isTranslation)
			std::cout << "WARNING:: 2 \n" ;
		if (RC & CRRDY_NoOrthCode)
			std::cout << "WARNING:: 3 \n" ;
	}
//	LOOKS WIERD...
//	RC = mmdb_N->GenerateSymMates ( NULL );

	if (RC>0)  {
		std::cout << "WARNING:: " << RC <<" \n" ;
	}

	return RC;
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

std::string
mmdb_helper::atom_resn( char * aid ) {
	std::string atom_str(aid);
	std::size_t found0	= atom_str.find_last_of ("(");
	std::size_t found1	= atom_str.find_first_of(")");
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
		std::cout << "INFO:: DONE \nINFO:: SET GRID SAMPLING" << std::endl;
	m_s_rate_=1.5;
	clipper::Grid_sampling gs( fphi_.spacegroup(), fphi_.cell(), fphi_.resolution(), m_s_rate_ );
	xmap_.init	( fphi_.spacegroup(), fphi_.cell(), gs ); 		// INITIALIZE MAP
//	grid_.init	( fphi_.spacegroup(), fphi_.cell(), fphi_.resolution(), m_s_rate_ );
//	xmap_.init	( fphi_.spacegroup(), fphi_.cell(), grid_ ); 		// INITIALIZE MAP
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
quaternion::rotate_particles( particles ps , gsl_vector *v0 ,std::vector<bool> mask) {
	if( is_complete() ){
		gsl_vector *tmp = gsl_vector_calloc(DIM);
		for(int i=0;i<ps.size();i++){
			if(!mask[i]){
				gsl_vector_sub( ps[i].second, v0 );
				rotate_coord  ( ps[i].second	 );
				gsl_vector_add( ps[i].second, v0 );
			}
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

double
residue_helper::calc_fi( particles rfull, int I ) { 
	double fi = 0.0, cosfi=0.0;
	gsl_vector *r0 = gsl_vector_calloc(DIM);

	for( int i=0 ; i<rfull.size() ; i++ ){
		gsl_vector *r = gsl_vector_calloc(DIM);
		gsl_vector_memcpy( r, rfull[i].second );
		if(	(I==2) && (i==CA_||i==CB_||i<=CG_||i==N_||i==C_||i==O_) ||
			(I==1) && (i==CA_||i==CB_|| i==N_||i==C_||i==O_) )
			continue;
		gsl_vector_add( r0, r );
	}
	gsl_matrix *O = gsl_matrix_calloc(DIM,DIM);
	switch( I ) {
		case 2:
			copyO2(O);
			break;
		case 1:
			copyO1(O);
			break;
		default:
			copyOS(O);
	}
	double nrm	= gsl_blas_dnrm2( r0 );
	gsl_vector_scale( r0 , 1.0/nrm );
	gsl_vector *d 	= gsl_vector_calloc(DIM);
	gsl_matrix_get_row( d , O , 1 );
	gsl_blas_ddot(d,r0,&cosfi);
	fi = acos(cosfi);
	gsl_matrix_free(  O );
	gsl_vector_free(  d );
	gsl_vector_free( r0 );

	return fi;
}

std::vector<bool> 
residue_helper::get_mask( int sw ) {
	std::vector<bool> mask;
	switch(sw) {
		case 2:
			for(int i=0;i<nResAtoms_;i++) {
				mask.push_back(i==CA_||i==CB_||i<=CG_||i==N_||i==C_||i==O_);
			}
			break;
		default:
			for(int i=0;i<nResAtoms_;i++) {
				mask.push_back(i==CA_||i==CB_||i==N_||i==C_||i==O_);
			}
			break;
	}
	return mask;
}

int
residue_helper::calc_proj(	clipper::Xmap<float> density, clipper::Grid_sampling gs, 
				std::vector< std::pair<int,double> > *theta,
				int nb , int which ) {
	rich::calc_map cmap;
	cmap.set_nbins(nb);

	gsl_matrix *P	= gsl_matrix_calloc(	nb, nb	 );
	gsl_matrix *CN	= gsl_matrix_calloc(	nb, nb	 );
	switch( which ) {
		case 2:
			if( cmap.proj01( P , CN , density ,
					 O2_ , v0_ , rc_ , zc_,
					 theta ) ) {
				std::cout << "ERROR:: FAILED" << std::endl;
				exit(1);
			}
			break;
		case 1:
			if( cmap.proj01( P , CN , density ,
					 O1_ , v0_ , rc_ , zc_,
					 theta ) ) {
				std::cout << "ERROR:: FAILED" << std::endl;
				exit(1);
			}
			break;
		default:
			if( cmap.proj01( P , CN , density ,
				OS_ , v0_ , rc_ , zc_,
				theta ) ) {
				std::cout << "ERROR:: FAILED" << std::endl;
				exit(1);
			}
			break;
	}
	gsl_matrix_free(P);
	gsl_matrix_free(CN);

	return 0;
}

int
residue_helper::calc_proj(	clipper::Xmap<float> density, clipper::Grid_sampling gs, 
				std::vector< std::pair<int,double> > *theta, int nb , 
				int which, int ires ) {
	rich::calc_map cmap;
	cmap.set_nbins(nb);

	gsl_matrix *P	= gsl_matrix_calloc(	nb, nb	 );
	gsl_matrix *CN	= gsl_matrix_calloc(	nb, nb	 );
	switch( which ) {
		case 2:
			if( cmap.proj01( P , CN , density ,
					 O2_ , v0_ , rc_ , zc_,
					 theta ) ) {
				std::cout << "ERROR::FAILED" << std::endl;
				exit(1);
			}
			break;
		case 1:
			if( cmap.proj01( P , CN , density ,
					 O1_ , v0_ , rc_ , zc_,
					 theta ) ) {
				std::cout << "ERROR::FAILED" << std::endl;
				exit(1);
			}
			break;
		default:
			if( cmap.proj01( P , CN , density ,
					 OS_ , v0_ , rc_ , zc_,
					 theta ) ) {
				std::cout << "ERROR::FAILED" << std::endl;
				exit(1);
			}
			break;
	}
	rich::mat_io	mio;
	mio.write_gsl2datn(P, CN, "res"+std::to_string(which)+"P"+ std::to_string(ires) +".dat" );
	gsl_matrix_free(P);
	gsl_matrix_free(CN);

	return 0;
}

particles 
particle_helper::particles_memcpy( particles rfull ) {
	particles rblank;
	for(int i=0;i<rfull.size();i++){
		gsl_vector *r = gsl_vector_calloc(DIM);
		gsl_vector_memcpy( r, rfull[i].second );
		particle ptmp;
		ptmp.first  = rfull[i].first;
		ptmp.second = r;
		rblank.push_back(ptmp);
	}
	return rblank;
}

bool compare_pid ( std::pair< int , double > pid1 ,
		   std::pair< int , double > pid2 ) {
	return( pid1.second < pid2.second ) ; 
};

std::vector<double> 
residue_helper::prune_angles(std::vector< std::pair<int,double> > * theta, double TOL, int N){ 
	int verbose=0;
	std::sort( (*(theta)).begin(), (*(theta)).end(), compare_pid );
	double val = 100.0, alimit = TOL; 
	std::vector<double> v_ang;
	int Ith = (*(theta)).size();
	do {
		val		= (*(theta))[Ith-1].second;
		double ang	= (*(theta))[Ith-1].first  * 360.0	/ (float(N));
		if( verbose )
			std::cout << "INFO::MAX > " << Ith << " " << ang << " " << val << std::endl;
		if( val > alimit )
			v_ang.push_back(ang);
		Ith--;
	} while( val > alimit && Ith>=0 ) ;
	return v_ang;
}

double
residue_helper::calc_OS( void ) {
	double zc=0.0;
	if(bVecAs_){
		rich::math_helper mah;
		zc = mah.gsl_calc_orth( n2_, n1_, c2_, c1_, OS_ );
		zc_=zc;	// DO OVERWRITE ZC
	}
	return zc;
}

double
residue_helper::calc_O1( particles resatms ) {
	double zc=0.0;
	if(bVecAs_&&bHaveA_&&bHaveB_){
		rich::math_helper mah;
		gsl_vector_memcpy(cb_,resatms[CB_].second);
		gsl_vector_memcpy(ca_,resatms[CA_].second);
		gsl_vector_memcpy(v0_,ca_);
		zc = mah.gsl_calc_orth(  cb_, ca_, n2_, n1_, O1_ ); // DO NOT OVERWRITE OLD ZC
	}
	return zc;
}

double
residue_helper::calc_O2( particles resatms ) {
	double zc=0.0;
	if(bVecAs_&&bHaveB_&&bHaveG_){
		rich::math_helper mah;
		zc = mah.gsl_calc_orth(  cg_, cb_, n2_, n1_, O2_ );
	}
	return zc;
}


double 
residue_helper::analyze_stage1( int imod ,int icha ,int ires , CMMDBManager *mmdb, int nb, particles *residue_atoms ) {

	rich::mmdb_helper mmhelp;
	PPCAtom atoms_T,atoms_T2;

	gsl_vector *vt	= gsl_vector_calloc(DIM);
	int skip_res	= 0;
	int abg		= 0;
	int nResidues	= mmdb->GetNumberOfResidues(imod,icha);
	int nAtoms,nAtoms2;

	mmdb->GetAtomTable    ( imod ,icha ,ires , atoms_T, nAtoms );
	residue_atoms->clear(); // warning dealloc particles OK but particle is NOT OK

	char r_inf[256];
	atoms_T[0]->GetAtomID(r_inf);
	std::string rtype = mmhelp.atom_resn(r_inf);
	if( ires<nResidues-1 &&  rtype != "PRO" ) {
		mmdb->GetAtomTable ( imod ,icha ,ires+1 , atoms_T2, nAtoms2 );
		gsl_vector_set( n2_, XX, atoms_T2[0]->x );
		gsl_vector_set( n2_, YY, atoms_T2[0]->y );
		gsl_vector_set( n2_, ZZ, atoms_T2[0]->z );
		skip_res++;
	}

	if(!bHasAMAT_) {
		A_= gsl_matrix_calloc( nAtoms, DIM );
	} else if ( bHasAMAT_&&A_->size1!=nAtoms ) {
		gsl_matrix_free(A_);
		A_= gsl_matrix_calloc( nAtoms, DIM );
	}

	for (int iat = 0; iat<nAtoms ; iat++ ) {
		char a_inf[256];
		atoms_T[iat]->GetAtomID(a_inf);
		std::string atype = mmhelp.atom_type(a_inf);
		std::string etype = mmhelp.atom_symb(a_inf);

		gsl_vector_set(vt,XX,atoms_T[iat]->x);
		gsl_vector_set(vt,YY,atoms_T[iat]->y);
		gsl_vector_set(vt,ZZ,atoms_T[iat]->z);

		if(iat==0) // SHOULD BE N
			gsl_vector_memcpy(v0_,vt);

		rich::particle res_atom;
		res_atom.second = gsl_vector_alloc(DIM);
		res_atom.first  = etype; 
		gsl_vector_memcpy(res_atom.second, vt);
		residue_atoms->push_back(res_atom);

		if( atype=="C" && etype=="C" ) {	// CA
			gsl_vector_memcpy(c1_,vt);
			skip_res++;
			C_=iat;
		}
		if( atype=="O" && etype=="O" ) {	// CB NOT GLY
			gsl_vector_memcpy(c2_,vt);
			skip_res++;
			O_=iat;
		}
		if( atype=="N" && etype=="N" ) {
			gsl_vector_memcpy(n1_,vt);
			skip_res++;
			N_=iat;
		}
		if( ires>=nResidues-1 || rtype == "PRO" ) {
			if( atype=="C" && etype=="C" ) {
				gsl_vector_memcpy(n2_,vt);
				skip_res++;
			}
		}
		if( atype=="CA" && rtype != "PRO" ) {
			gsl_vector_memcpy(ca_,vt);
			abg++;
			bHaveA_ = true;
			CA_=iat;
		}
		if( atype=="CB" && rtype != "PRO" ) {
			gsl_vector_memcpy(cb_,vt);
			abg++;
			bHaveB_ = true;
			CB_=iat;
		}
		if( atype=="CG" && rtype != "PRO" ) {
			gsl_vector_memcpy(cg_,vt);
			abg++;
			bHaveG_ = true;
			CG_=iat;
		}
		gsl_matrix_set_row(A_,iat,vt);
	}
	if(skip_res!=4) {
		bSkip_ = true;
		std::cout << "INFO::SKIPPING RESIDUE "<< ires << " HAS " << skip_res << std::endl;
	}
//	CALCULATE LIMITS
	gsl_matrix *V	= gsl_matrix_calloc( DIM, DIM	 );
	gsl_vector *S	= gsl_vector_calloc(	DIM 	 );
	gsl_vector *wrk	= gsl_vector_calloc(	DIM	 );
	gsl_matrix *A	= gsl_matrix_calloc( nAtoms, DIM );
	gsl_matrix_memcpy( A , A_ );
	nResAtoms_=nAtoms;
	gsl_linalg_SV_decomp( A , V , S , wrk );
	double cutoff 	= sqrt(gsl_vector_get(S,0))*0.5;
	bVecAs_=true;
	rc_=cutoff;

	gsl_matrix_free(  A  );
	gsl_matrix_free(  V  );
	gsl_vector_free(  S  );
	gsl_vector_free( wrk );
	gsl_vector_free(  vt );

	return cutoff;
}
