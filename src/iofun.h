#ifndef IOFUN_H
#define IOFUN_H
#include <cstdlib>
#include <string>

#include "types.h"

namespace rich {

	class helper_ops {
		public:
			std::string	split_string  ( const std::string& , const std::string& , int ) ;
			void 		split_filename( const std::string& );
		private:
			;
	};

	class tensorIO {
		public:
			inline tensorIO() {}  
			void output_matrix_label(gsl_matrix *, gsl_vector *v);	
			void output_matrix(gsl_matrix *O);
			void output_vector(gsl_vector *v);
			void output_geometry( particles px );
	};

	class map_io : public helper_ops {
		public:
			map_io(void) { bHaveHeader_=false; };
			void print_help( void ) { std::cout << "XMAP IO" << std::endl; };
			void set_path( std::string );
			void set_selection( std::string );
			void read_map_header( void );
			void read_map_header( std::string );
			clipper::String	get_ifil(void) { return    ipfile_; };
			clipper::String	get_isel(void) { return     ipcol_; };
			clipper::String	get_itot(void) { return  fphi_str_; };
		private:
			bool bHaveHeader_;
			std::string		path_;
			std::string		selection_;
			clipper::String		fphi_str_;
			clipper::CCP4MTZfile	mtzin_;
			clipper::String		ipfile_,ipcol_;
	};

	class mat_io {
		public:
			int write_gsl2dat  ( gsl_matrix *, std::string );
			int write_gsl2datn ( gsl_matrix *, gsl_matrix *, std::string );
			int write_vdbl2dat ( std::vector<double> , std::string );
		private:
			int N_;
	};
}

#endif
