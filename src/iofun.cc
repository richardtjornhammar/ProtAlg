#include <cstdlib>
#include <string>
#include <iostream>
#include <algorithm>
#include <fstream>
#include <sstream>
#include "iofun.h"

using namespace rich;

void
map_io::set_path( std::string it_str ) {
	path_		= it_str;
}

void
map_io::set_selection( std::string it_str ) {
	selection_	= it_str;
}

void 
helper_ops::split_filename( const std::string& str)
{
  std::cout << "Splitting: " << str << '\n';
  std::size_t found = str.find_last_of("\\/");
  std::cout << " path: " << str.substr(0,found) << '\n';
  std::cout << " file: " << str.substr(found+1) << '\n';
}

std::string 
helper_ops::split_string ( const std::string& str, const std::string& split_str, int i)
{
	std::size_t found = str.find_last_of(split_str);
	std::string retstr;
	switch(i){
		case 1:
			retstr = str.substr(found+1);
			break;
		default:
			retstr = str.substr(0,found);
	}
	return retstr;	
}

void
map_io::read_map_header( std::string ifile ) {

	int		is_mtz_file=0;
	int		verbose=1;
	char		chr;
	clipper::String	ipcolf;
	clipper::String ipfile;
	ipfile = ifile;

	std::string	filename[3];

	std::cout << "INFO:: HAVE INPUT: \t" << ipfile << std::endl;
	std::string tmpline( ipfile );
	filename[0] = tmpline;
	std::size_t found = tmpline.find("mtz");
	if ( found!=std::string::npos ){
		std::cout << "INFO:: FOUND MTZ INPUT NAME" << std::endl;
		try { 
			mtzin_.open_read( ipfile );	
			is_mtz_file = 1;
		}
		catch ( ... ) {
			std::cout << "ERROR:: NOT A VALID MTZ FILE: " << filename[0] << std::endl;
			is_mtz_file = 0;
			exit(1);
		} 
		std::string label;
		std::string type;
		std::string mtypF("F"),mtypP("P");
		std::string mlab("WT");
		std::string selection("/[");
		std::string base_str;

		if (is_mtz_file) { 
			int nHit=0;
			mtzin_.set_column_label_mode( clipper::CCP4MTZfile::Legacy );
			std::vector<clipper::String> v = mtzin_.column_labels();
			if (v.size() > 1) { 
				int read_success = 1;
				for (unsigned int i=0; i<v.size(); i++) {
					if(verbose)
    						std::cout << i << " \t " << v[i] << "\n";
   						std::string::size_type ispace = v[i].find_last_of(" ");
   						if (ispace == std::string::npos) {
						std::cout <<  "WARNING:: uninterprettable label \"" 
						<< v[i] << "\" of "<< filename[0] << "\n";
					} else {
						label = v[i].substr(0, ispace);
						type  = v[i].substr(ispace+1);
						std::size_t found = label.find(mlab);
						if(found!=std::string::npos && nHit < 2) {
							nHit++;
							base_str 	 = split_string (label,"\\/",0);
							std::string toplabel_str = split_string (label,"\\/",1);
							selection	+= toplabel_str;
							selection	+= (nHit==1)?(","):("]");
						}
						if( verbose==2 ) {
							std::cout << "Got label :" << label 
							<< ": and type :" << type << ":\n";
							split_filename (label);	
						}
					}
				}
			}
		}
		mtzin_.close_read();
		path_		= base_str;
		selection_	= selection;
		fphi_str_= base_str+selection;
		std::cout << "INFO:: CAN GET DATA FROM: " << fphi_str_ << std::endl;	
		std::cout << "INFO:: DO YOU WANT TO USE THIS STRING? [Y/n]\t";
		while( true ) {
			std::cin >> chr;
			if( chr=='Y' || chr=='n' ) {
				break;
			} else {
				std::cout << "PLEASE ENTER AGAIN [Y/n]: \n";
			}
		}
		if( chr=='Y' )
			ipcolf=fphi_str_;
		if( chr=='n' ) {
			std::cout << "PLEASE ENTER COLUMN/S STRING (EX.: /FOO/BAR/[FWT,PHWT] ):\n";
			std::string colname;
			std::cin.ignore();
			std::getline(std::cin, colname);
			ipcolf	 = colname;
			fphi_str_ = colname;
		}
	}
	ipcol_ = ipcolf;
	ipfile_= ipfile;
	bHaveHeader_=true;

	if(verbose==2)
		std::cout << "INFO:: DONE"<< std::endl;
}


void 
tensorIO::output_vector(gsl_vector *v){
	std::cout <<"::INFO::VECTOR::"<< std::endl;
	for(int i=0;i<v->size;i++)
		std::cout << " " << gsl_vector_get(v,i);
	std::cout << std::endl;
}

void 
tensorIO::output_matrix_label(gsl_matrix *M, gsl_vector *v){
	std::cout <<"::INFO::MATRIX::"<< std::endl;
	for(int i=0;i<M->size2;i++){
		if(M->size2==v->size)
			std::cout << " " << gsl_vector_get(v,i);
		for(int j=0;j<M->size1;j++)
			std::cout << " " << gsl_matrix_get(M,j,i);
		std::cout << std::endl;	
	}
	std::cout << std::endl;
}

void 
tensorIO::output_matrix(gsl_matrix *M){
	std::cout <<"::INFO::MATRIX::"<< std::endl;
	for(int i=0;i<M->size2;i++){
		for(int j=0;j<M->size1;j++)
			std::cout << " " << gsl_matrix_get(M,j,i);
		std::cout << std::endl;	
	}
	std::cout << std::endl;
}



void 
tensorIO::output_geometry(particles px) {
	int n=px.size();
	std::cout << n << std::endl;
	std::cout << "FLUSHED COORDS" << std::endl;
	for(int i=0;i<n;i++){
		std::cout << " " << px[i].first << " " << gsl_vector_get(px[i].second,XX) << " " 
			<< gsl_vector_get(px[i].second,YY) << " " << gsl_vector_get(px[i].second,ZZ) << std::endl; 
	}
}

int
mat_io::write_gsl2dat ( gsl_matrix *A, std::string fstr ) {
	const char *c_filename = fstr.c_str();
	std::ofstream outf;
	outf.open(c_filename);
	for( int i=0 ; i<A->size1 ; i++ ) {
	for( int j=0 ; j<A->size2 ; j++ ) {
		outf	<< " " << gsl_matrix_get(A,i,j); 
	}
		outf	<< std::endl;
	}
	outf.close();
	return 0;
}

int
mat_io::write_vdbl2dat( std::vector<double> vdbl, std::string fstr ){
	const char *c_filename = fstr.c_str();
	std::ofstream outf;
	outf.open(c_filename);
	for( int i=0 ; i<vdbl.size() ; i++ ) {
		outf	<< " " << i << " " << vdbl[i] << std::endl; 
	}
	outf.close();
	return 0;
}

int
mat_io::write_gsl2datn( gsl_matrix *A, gsl_matrix *N, std::string fstr ) {
	const char *c_filename = fstr.c_str();
	std::ofstream outf;
	outf.open(c_filename);
	for( int i=0 ; i<A->size1 ; i++ ) {
	for( int j=0 ; j<A->size2 ; j++ ) {
		double cnt = gsl_matrix_get(N,i,j);
		if(cnt>0)
			outf	<< " " << sqrt(gsl_matrix_get(A,i,j))/cnt; 
		else
			outf	<< " " << 0.0; 
	}
		outf	<< std::endl;
	}
	outf.close();
	return 0;
}
