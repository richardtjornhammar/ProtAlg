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

#define XX	0
#define YY	1
#define ZZ	2
#define	DIM	3
#define _VERBOSE_ 0

namespace rich {

	typedef std::pair<std::string, gsl_vector * >	particle;
	typedef std::vector< particle > 		particles;

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
		private:
			;
	};

	class math_helper {
		public:
			math_helper(){ bHelp_=false; };
			void	gsl_cross3D( gsl_vector * , gsl_vector *, gsl_vector *);
		private:
			bool bHelp_;
	};

	class quaternion {
		public:
			inline	quaternion() { bSet_= false;  q_= gsl_vector_calloc(4); R_= gsl_matrix_calloc(DIM,DIM); };
			void	clear(void)  { gsl_vector_set_zero(q_); bSet_= false; };
			int	assign_quaterion( gsl_vector *x, double angle );
			int	rotate_coord	( gsl_vector *x );
			int	rotate_particles( particles );
			int	rotate_particles( particles , gsl_vector * );
			bool	is_complete(){ return bSet_; };
			~quaternion() {
    				gsl_vector_free( q_ );
				gsl_matrix_free( R_ );
			};
		private:
			bool bSet_;
			gsl_vector *q_; // quaternion 
			gsl_matrix *R_;
	};
}

#endif
