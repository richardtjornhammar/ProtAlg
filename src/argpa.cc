#include <iostream>
#include <algorithm>
#include <cstdlib>
#include <vector>
#include <string>
#include "argpa.h"
#include "types.h"

using namespace rich;

void	
arg_parser::set_options(std::vector < std::pair < int, std::pair<std::string, std::string > > > opts_n_def ) {
	opts_n_def_.clear();
	bDefSet_=false;
	for( int i=0 ; i<opts_n_def.size() ; i++ ) {
		if( !bDefSet_ ){
			opts_n_def_.push_back( opts_n_def[i] );
		}
	}
	bDefSet_=true;
}

void
arg_parser::print_help( int failed , std::pair< int, std::vector< std::string > > ret_args ){
	if(bDefSet_) {
		std::cout << "HELP:: VIABLE INPUT OPTIONS ARE:\nHELP:: [X] \t OPTION \t VALUE\n" << std::endl;
		bool bShow=ret_args.second.size()==opts_n_def_.size();
		for(int i=0;i<opts_n_def_.size();i++){
			std::cout << "HELP:: [" << opts_n_def_[i].first << "] \t " 
				<< opts_n_def_[i].second.first << " \t " ;
			if( bShow )
				std::cout << ret_args.second[i];
			std::cout << std::endl;
		}
		ret_args.first = 1;
		std::cout << std::endl;	
		if(failed==1) {
			std::cout << "ERROR::OBLIGATORY ARGUMENTS NOT SET ( MARKED [1] )\n" << std::endl;
		}

	} else {
		std::cout << "ERROR::ARGUMENT OPTIONS NOT SET\n" << std::endl;
	}
}

void
arg_parser::print_help( int failed ){
	if(bDefSet_) {
		std::cout << "HELP:: VIABLE INPUT OPTIONS ARE:\nHELP:: [X] \t OPTION \t VALUE\n" << std::endl;
		for(int i=0;i<opts_n_def_.size();i++){
			std::cout << "HELP:: [" << opts_n_def_[i].first << "] \t " 
				<< opts_n_def_[i].second.first << " \t " ;
			std::cout << std::endl;
		}
		std::cout << std::endl;	
		if (failed==1) {
			std::cout << "ERROR::OBLIGATORY ARGUMENTS NOT SET ( MARKED [1] )\n" << std::endl;
		}

	} else {
		std::cout << "ERROR::ARGUMENT OPTIONS NOT SET\n" << std::endl;
	}
}

void
arg_parser::print_help( void ){
	if(bDefSet_) {
		std::cout << "HELP:: VIABLE INPUT OPTIONS ARE:\nHELP:: [X] \t OPTION \t VALUE\n" << std::endl;
		for(int i=0;i<opts_n_def_.size();i++) {
			std::cout << "HELP:: [" << opts_n_def_[i].first << "] \t " 
				<< opts_n_def_[i].second.first << " \t " ;
			std::cout << std::endl;
		}
		std::cout << std::endl;	

	} else {
		std::cout << "ERROR::ARGUMENT OPTIONS NOT SET => PARSER NOT INITIALIZED\n" << std::endl;
	}
}

std::pair< int , std::vector< std::string > >
arg_parser::parse ( std::pair < int, char ** > mp ) {
	int failed=0;
	std::vector< std::string > args;
	std::pair  < int, std::vector< std::string > > ret_args;
	std::vector< int > v_set;

	for( int i=0 ; i<mp.first ; i++ ) {
		args.push_back			( mp.second[i] );
	}

	if( bDefSet_ ) {
		for( int i=0 ; i<opts_n_def_.size() ; i++ ) {
			ret_args.second.push_back ( opts_n_def_[i].second.second   );
			v_set.push_back(0);
		}
	} else {
		failed=1;
	}

	// COMMAND INPUT
	int arg 	= 0 ;
	ret_args.first 	= 0 ;
	int a_size = mp.first==1?2:mp.first;
	int help=0;
	while ( ++arg < a_size  ) {
		if(  args[arg] == "-h" ||  args[arg] =="--h" || args[arg] =="--help" || args[arg] =="-help" || mp.first==1) {
			help=1;
			break;
		}else{
			int ia=arg,iap=++arg;
			for(int i=0;i< opts_n_def_.size();i++) {
				if ( args[ia] == opts_n_def_[i].second.first ) {
					if(args[iap][0]=='-')
						failed=1;
					ret_args.second[i] = args[iap];
					v_set[i]=1;
				}
			}
		}
	}

	for(int i=0;i<v_set.size();i++) {
		if(  opts_n_def_[i].first > 0 && v_set[i] == 0 ) {
			ret_args.first = 1;
			failed=1;
			break;
		}
	}

	if ( help==1 || failed==1 ) {
		print_help( failed , ret_args );
	}

	return ret_args;
}

std::pair< int , std::vector< std::string > >
arg_parser::parse ( std::pair < int, char ** > mp, std::vector < std::pair < int, std::pair<std::string, std::string > > > opts_n_defaults ) {

	// ARGS OUTPUT ORDER CORRESPONDS TO OPTS_N_DEFAULTS FIRST VALUE
	std::vector< std::string > args;
	std::pair  < int, std::vector< std::string > > ret_args;
	std::vector< int > v_set;
	int failed=0;

	for( int i=0 ; i<mp.first ; i++ ) {
		args.push_back			( mp.second[i] );
	}
	for( int i=0 ; i<opts_n_defaults.size() ; i++ ) {
		ret_args.second.push_back	( opts_n_defaults[i].second.second );
		if( !bDefSet_ ){
			opts_n_def_.push_back( opts_n_defaults[i] );
		}
		v_set.push_back(0);
	}
	bDefSet_=true;

	// COMMAND INPUT
	int arg 	= 0 ;
	ret_args.first 	= 0 ;
	int a_size = mp.first==1?2:mp.first;
	int help=0;
	while ( ++arg < a_size  ) {
		if(  args[arg] == "-h" ||  args[arg] =="--h" || args[arg] =="--help" || args[arg] =="-help" || mp.first==1) {
			help=1;
			break;
		}else{
			int ia=arg,iap=++arg;
			for(int i=0;i<opts_n_defaults.size();i++) {
				if ( args[ia] == opts_n_defaults[i].second.first ) {
					if(args[iap][0]=='-')
						failed=1;
					ret_args.second[i] = args[iap];
					v_set[i]=1;
				}
			}
		}
	}

	for(int i=0;i<v_set.size();i++) {
		if( opts_n_defaults[i].first > 0 && v_set[i] == 0 ) {
			ret_args.first = 1;
			failed=1;
			break;
		}
	}


	if ( help==1 || failed==1 ) {
		print_help( failed , ret_args );
	}

	return ret_args;
}
