#ifndef RTYPES_H
#define RTYPES_H

#include <cstdlib>
#include <vector>
#include <string>
#include <math.h>

//// GSL STUFF
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_permutation.h>

//	CLIPPER	STUFF
#include <clipper/clipper.h>
#include <clipper/clipper-ccp4.h>
#include <clipper/clipper-contrib.h>

//	MMDB	STUFF
#include "mmdb/mmdb_manager.h"

#define XX	0
#define YY	1
#define ZZ	2
#define	DIM	3
#define _VERBOSE_ 0

namespace rich {

	typedef std::pair<std::string, gsl_vector * >	particle;
	typedef std::vector< particle > 		particles;

	class particle_helper {
		public:		particles particles_memcpy( particles );
		private:
			;
	};

	class phys_constants {
		public:
			phys_constants() { // SI UNITS
				KBOL_ = 1.38064889E-23	;
				RGAS_ = 8.31446210	;
				NAVO_ = 6.022140857E23	;	
			};
		//! Linear algebra routines needed for fitting
			double	kB		( void ) { return KBOL_; };
			double  R		( void ) { return RGAS_; };
			double	Na		( void ) { return NAVO_; };
		private:
			double KBOL_ ;		//	m2 kg s-2 K-1
			double RGAS_ ;		//	J K−1 mol−1
			double NAVO_ ;
	};

	class map_manager {
		public:
			map_manager(){};
			int	assign_map( std::string );
			clipper::MTZdataset 		get_set(void) { return   set_; } ;
			clipper::MTZcrystal  		get_xtl(void) { return   xtl_; } ;
			clipper::Resolution		get_res(void) { return  reso_; } ;
			clipper::Grid_sampling		get_gsa(void) { return  grid_; } ;
			clipper::HKL_info		get_nfo(void) { return	 hkl_; } ;
			clipper::HKL_data< clipper::datatypes::F_phi<float> > get_fph() { return fphi_; } ;
			clipper::Xmap<float>		get_map(void) { return xmap_; } ;
		private:
			clipper::MTZdataset 	set_; 
			clipper::MTZcrystal  	xtl_;
			clipper::Resolution	reso_;
			clipper::Grid_sampling	grid_;
			clipper::HKL_info	hkl_;
			clipper::HKL_data< clipper::datatypes::F_phi<float> > fphi_;
			clipper::Xmap<float>	xmap_;
			float m_s_rate_;
	};

	class mmdb_helper {
		public:
			std::string	atom_type( char * );
			std::string	atom_symb( char * );
			std::string	atom_resn( char * );
			int	update_residue	( int, int, int, CMMDBManager*, particles );
			int	check_clash	( int, int, int, CMMDBManager*, particles , double );
			int	copy_xtal	( CMMDBManager*, CMMDBManager* );
		private:
			;
	};

	class math_helper {
		public:
			math_helper(){ bHelp_=false; };
			void	gsl_cross3D	( gsl_vector * , gsl_vector *, gsl_vector * );
			double	gsl_calc_orth	( gsl_vector * , gsl_vector *, gsl_vector *,
						  gsl_vector * , gsl_matrix * );
		private:
			bool bHelp_;
	};

	class quaternion {
		public:
			inline	quaternion() {	bSet_= false;  q_= gsl_vector_calloc(DIM+1); 
						R_= gsl_matrix_calloc(DIM,DIM); };
			void	clear(void)  {	gsl_vector_set_zero(q_); bSet_= false; };
			int	assign_quaterion( gsl_vector *x, double angle );
			int	rotate_coord	( gsl_vector *x );
			int	rotate_particles( particles );
			int	rotate_particles( particles , gsl_vector * );
			int	rotate_particles( particles , gsl_vector * , std::vector<bool> );
			void	print();
			bool	is_complete() { return bSet_; };
			~quaternion() {
    				gsl_vector_free( q_ );
				gsl_matrix_free( R_ );
			};
		private:
			bool bSet_;
			gsl_vector *q_; // quaternion 
			gsl_matrix *R_;
	};

	class residue_helper {
		public:
			residue_helper() {
				n2_=gsl_vector_calloc(DIM);	n1_=gsl_vector_calloc(DIM);
				c2_=gsl_vector_calloc(DIM);	c1_=gsl_vector_calloc(DIM);
				ca_=gsl_vector_calloc(DIM);	cb_=gsl_vector_calloc(DIM);
				cg_=gsl_vector_calloc(DIM);	v0_=gsl_vector_calloc(DIM);
				OS_=gsl_matrix_calloc(DIM,DIM);	O2_=gsl_matrix_calloc(DIM,DIM);
				O1_=gsl_matrix_calloc(DIM,DIM);
				isPRO_	= false; bOS_	= false; bHasAMAT_	= false;
				bHaveA_	= false; bHaveB_= false; bHaveG_	= false; 
				bVecAs_	= false; bSkip_	= false; nResAtoms_ 	= 0;
			};
			double	analyze_stage1( int , int , int , CMMDBManager * , int , particles * );
			double	calc_OS( void );
			double	calc_O1( particles  ); // particles only r no rw
			double	calc_O2( particles  ); // particles only r no rw
			std::vector<double> prune_angles(std::vector< std::pair<int,double> > *, double, int);
			int	calc_proj( clipper::Xmap<float>, clipper::Grid_sampling gs, std::vector< std::pair<int,double> > * , int, int );
			int	calc_proj( clipper::Xmap<float>, clipper::Grid_sampling gs, std::vector< std::pair<int,double> > * , int, int, int );
			void copyA (  gsl_matrix *A ){ if(A->size1==A_->size1&&A->size2==A_->size2){ gsl_matrix_memcpy(A,A_); }; };
			void copyv0( gsl_vector *v0 ){ if(v0->size==v0_->size){ gsl_vector_memcpy(v0,v0_); }; };
			void copyOS( gsl_matrix *OS ){ if(OS->size1==OS_->size1&&OS->size2==OS_->size2){ gsl_matrix_memcpy(OS,OS_); }; };
			void copyO1( gsl_matrix *O1 ){ if(O1->size1==O1_->size1&&O1->size2==O1_->size2){ gsl_matrix_memcpy(O1,O1_); }; };
			bool skip(void) { return bSkip_; };
			bool do2nd(void) { return bHaveA_&&bHaveB_; }
			bool do3nd(void) { return bHaveB_&&bHaveG_; }
			std::vector<bool> get_mask( int );
			void scale_c(double s){zc_*=s;rc_*=s;};
			void scale_z(double s){zc_*=s;};
			~residue_helper() {
				gsl_vector_free( n2_ ); gsl_vector_free( n1_ );
				gsl_vector_free( c2_ ); gsl_vector_free( c1_ );
				gsl_vector_free( ca_ ); gsl_vector_free( cb_ );
				gsl_vector_free( cg_ ); gsl_vector_free( v0_ );
				gsl_matrix_free( OS_ ); gsl_matrix_free( O2_ );
				gsl_matrix_free( O1_ );
			};
		private:
			bool bSkip_;
			char resT_;
			bool isPRO_, bHaveA_, bHaveB_, bHaveG_, bOS_;
			gsl_vector *n2_,*n1_,*c2_,*c1_,*ca_,*cb_,*cg_,*v0_;
			bool bVecAs_;
			gsl_matrix *OS_,*O1_,*O2_;
			gsl_matrix *A_;
			bool bHasAMAT_;
			double zc_,rc_;
			std::vector<int> baplus_;
			std::vector<int> gbplus_;
			std::vector< particles > confs0_;
			int N_, CA_, CB_, CG_, C_, O_;
			int nResAtoms_;
	};
}

#endif
