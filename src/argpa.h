#ifndef ARGPAR_H
#define ARGPAR_H

#include <cstdlib>
#include <vector>
#include <string>

namespace rich {

	class arg_parser {
		public:
			arg_parser( void ) { bDefSet_=false; };

			void print_help( int , std::pair< int, std::vector< std::string > > );
			void print_help( int );
			void print_help( void );

			std::pair< int , std::vector< std::string > >  
				parse ( std::pair < int, char ** > mp ,
					std::vector < std::pair < int, std::pair<std::string, std::string > > > opts_n_defaults );
			std::pair< int , std::vector< std::string > >  
				parse ( std::pair < int, char ** > mp );

			void	set_options(std::vector < std::pair < int, std::pair<std::string, std::string > > > opts_n_def );

		private:
			std::vector < std::pair < int, std::pair<std::string, std::string > > > opts_n_def_;
			bool bDefSet_;
	};
}

#endif
