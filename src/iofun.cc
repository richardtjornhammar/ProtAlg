#include <cstdlib>
#include <string>
#include <iostream>
#include <algorithm>
#include <fstream>
#include <sstream>
#include <iomanip>
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


void 
fileIO::output_geometry( std::string filename, particles px ){
	const char *c_filename = filename.c_str();
	std::ofstream outp_coord;
	outp_coord.open(c_filename);

	int n=px.size();
	outp_coord << n << std::endl;
	outp_coord << "FLUSHED COORDS" << std::endl;
	for(int i=0;i<n;i++){
		outp_coord 	<< " " << px[i].first << " " << gsl_vector_get(px[i].second,XX) << " " 
				<< gsl_vector_get(px[i].second,YY) << " " << gsl_vector_get(px[i].second,ZZ) << std::endl; 
	}
	outp_coord.close();
}

void 
fileIO::output_geometry( std::string filename, particles px, std::string label ){
	const char *c_filename = filename.c_str();
	std::ofstream outp_coord;
	outp_coord.open(c_filename);

	int n=px.size();
	outp_coord << n << std::endl;
	outp_coord << label << std::endl;
	for(int i=0;i<n;i++){
		outp_coord	<< " " << px[i].first << " " << gsl_vector_get(px[i].second,XX) << " " 
				<< gsl_vector_get(px[i].second,YY) << " " << gsl_vector_get(px[i].second,ZZ) << std::endl; 
	}
	outp_coord.close();	
}

particles 
fileIO::read_xyz(std::string filename) {
	const char *c_filename = filename.c_str();
	std::ifstream inp_molecule;
	inp_molecule.open(c_filename);
	int Nparts = 0;
	std::string title;
	std::vector<double> BOX;
	BOX.push_back(1.0); BOX.push_back(1.0); BOX.push_back(1.0);
	particles ps;
	if(inp_molecule.is_open()){
		std::string iline;
		getline(inp_molecule,iline);
		std::stringstream data(iline);
		data >> Nparts;
		getline(inp_molecule, iline);
		title = iline;
		std::stringstream atad(iline);
		std::string wrd; atad >> wrd;
		if( wrd == "BOX" ){ // CAN SET A SCALE FOR COORDINATES (ORTHOGONAL SPACE)
			atad >> BOX[XX]; atad >> BOX[YY]; atad >> BOX[ZZ];
		}

		while( !(inp_molecule.eof()) ){
//	DECLARATIONS
			double r[3];
			particle sp;
			std::string atom;
			gsl_vector *rvec = gsl_vector_calloc(DIM);
//	ACTUAL IO
			getline(inp_molecule,iline);
			std::stringstream    data(iline);
			data >> atom >> r[0] >> r[1] >> r[2];
			gsl_vector_set(rvec,XX,r[XX]*BOX[XX]);
			gsl_vector_set(rvec,YY,r[YY]*BOX[YY]);
			gsl_vector_set(rvec,ZZ,r[ZZ]*BOX[ZZ]);
			sp.first = atom; sp.second = rvec;
			if(sp.first=="")
				continue;
			ps.push_back(sp);
		}
		inp_molecule.close();
	}else{
		std::cout << "FATAL:: COULD NOT OPEN FILE " << std::endl;
	}

	return ps;
}


void 
fileIO::output_pdb( std::string filename, particles px, std::vector<int> ndx ) {
//NOTE:: pymol can color the clusters with CMD:: util.cbc
	const char *c_filename = filename.c_str();
	std::ofstream outp_coord;
	outp_coord.open(c_filename);
	std::string title = "COMPND    CLUSTER";
	std::string autho = "AUTHOR    GENERATED BY RICHTOOL";
	std::string fin1  = "MASTER        0    0    0    0    0    0    0    0";
	std::stringstream    fins;
	fins << std::setw(5) << ((int)round(px.size())) << std::setw(5) << 0;
	fins << std::setw(5) << ((int)round(px.size())) << std::setw(5) << 0;

	if( ndx.size() == px.size() ){
		outp_coord << title << std::endl;
		outp_coord << autho << std::endl;
		for(int i=0;i<px.size();i++){
			outp_coord << std::right<< std::setw(6) << "HETATM" << std::setw(5) << i+1 << std::setw(3)
				<< px[i].first	<< std::setw(6) << "LIG" << std::setw(2) << ( (char) (65+ndx[i]) )
				<< std::setw(4) << 1 << "    "	<< std::setprecision(4) 
				<< std::setw(8) << gsl_vector_get(px[i].second,XX) 
				<< std::setw(8) << gsl_vector_get(px[i].second,YY)
				<< std::setw(8) << gsl_vector_get(px[i].second,ZZ) << std::setprecision(3)
				<< std::setw(6) << 1.00 << std::setw(6) << 0.00 << std::setw(12) << px[i].first << std::endl;
		}
	}
	outp_coord << fin1 << fins.str() << std::endl;
	outp_coord << "END" << std::endl;
	outp_coord.close();
}

void 
fileIO::output_pdb( std::string filename, particles px, std::vector<int> ndx , std::vector<int> order) {
//NOTE:: pymol can color the clusters with CMD:: util.cbc
	const char *c_filename = filename.c_str();
	std::ofstream outp_coord;
	outp_coord.open(c_filename);
	std::string title = "COMPND    CLUSTER";
	std::string autho = "AUTHOR    GENERATED BY RICHTOOL";
	std::string fin1  = "MASTER        0    0    0    0    0    0    0    0";
	std::stringstream    fins;
	fins << std::setw(5) << ((int)round(px.size())) << std::setw(5) << 0;
	fins << std::setw(5) << ((int)round(px.size())) << std::setw(5) << 0;

	if( ndx.size() == px.size() && order.size() == ndx.size() ){
		outp_coord << title << std::endl;
		outp_coord << autho << std::endl;
		for(int i=0;i<px.size();i++){
			outp_coord << std::right<< std::setw(6) << "HETATM" << std::setw(5) << i+1 << std::setw(3)
				<< px[order[i]].first	<< std::setw(6) << "LIG" << std::setw(2) << ( (char) (65+ndx[order[i]]) )
				<< std::setw(4) << 1 << "    "	<< std::setprecision(4) 
				<< std::setw(8) << gsl_vector_get(px[order[i]].second,XX) 
				<< std::setw(8) << gsl_vector_get(px[order[i]].second,YY)
				<< std::setw(8) << gsl_vector_get(px[order[i]].second,ZZ) << std::setprecision(3)
				<< std::setw(6) << 1.00 << std::setw(6) << 0.00 << std::setw(12) << px[order[i]].first << std::endl;
		}
	}
	outp_coord << fin1 << fins.str() << std::endl;
	outp_coord << "END" << std::endl;
	outp_coord.close();
}


void 
fileIO::output_pdb( std::string filename, particles px, std::vector<int> ndx , std::string alabel) {
//NOTE:: pymol can color the clusters with CMD:: util.cbc
	const char *c_filename = filename.c_str();
	std::ofstream outp_coord;
	outp_coord.open(c_filename);
	std::string title = "COMPND    CLUSTER";
	std::string autho = "AUTHOR    GENERATED BY RICHTOOL";
	std::string fin1  = "MASTER        0    0    0    0    0    0    0    0";
	std::stringstream    fins;
	fins << std::setw(5) << ((int)round(px.size())) << std::setw(5) << 0;
	fins << std::setw(5) << ((int)round(px.size())) << std::setw(5) << 0;

	if( ndx.size() == px.size() ){
		outp_coord << title << std::endl;
		outp_coord << autho << std::endl;
		for(int i=0;i<px.size();i++){
			outp_coord << std::right<< std::setw(6) << "HETATM" << std::setw(5) << i+1 << std::setw(3)
				<< alabel	<< std::setw(6) << "LIG" << std::setw(2) << ( (char) (65+ndx[i]) )
				<< std::setw(4) << 1 << "    "	<< std::setprecision(4) 
				<< std::setw(8) << gsl_vector_get(px[i].second,XX) 
				<< std::setw(8) << gsl_vector_get(px[i].second,YY)
				<< std::setw(8) << gsl_vector_get(px[i].second,ZZ) << std::setprecision(3)
				<< std::setw(6) << 1.00 << std::setw(6) << 0.00 << std::setw(12) << alabel << std::endl;
		}
	}
	outp_coord << fin1 << fins.str() << std::endl;
	outp_coord << "END" << std::endl;
	outp_coord.close();
}

void 
fileIO::output_pdb( std::string filename, particles px ) {
//NOTE:: pymol can color the clusters with CMD:: util.cbc
	const char *c_filename = filename.c_str();
	std::ofstream outp_coord;
	outp_coord.open(c_filename);
	std::string title = "COMPND    CLUSTER";
	std::string autho = "AUTHOR    GENERATED BY RICHTOOL";
	std::string fin1  = "MASTER        0    0    0    0    0    0    0    0";
	std::stringstream    fins;
	fins << std::setw(5) << ((int)round(px.size())) << std::setw(5) << 0;
	fins << std::setw(5) << ((int)round(px.size())) << std::setw(5) << 0;

	if(  px.size() > 0 ){
		outp_coord << title << std::endl;
		outp_coord << autho << std::endl;
		for(int i=0;i<px.size();i++){
			outp_coord << std::right<< std::setw(6) << "HETATM" << std::setw(5) << i+1 << std::setw(3)
				<< px[i].first	<< std::setw(6) << "LIG" << std::setw(2) << ( 'A' )
				<< std::setw(4) << 1 << "    "	<< std::setprecision(4) 
				<< std::setw(8) << gsl_vector_get(px[i].second,XX) 
				<< std::setw(8) << gsl_vector_get(px[i].second,YY)
				<< std::setw(8) << gsl_vector_get(px[i].second,ZZ) << std::setprecision(3)
				<< std::setw(6) << 1.00 << std::setw(6) << 0.00 << std::setw(12) << px[i].first << std::endl;
		}
	}
	outp_coord << fin1 << fins.str() << std::endl;
	outp_coord << "END" << std::endl;

	outp_coord.close();
}

void 
fileIO::output_pdb( std::string filename, particles px, gsl_vector *w ) {
//NOTE:: pymol can color the clusters with CMD:: util.cbc
	const char *c_filename = filename.c_str();
	std::ofstream outp_coord;
	outp_coord.open(c_filename);
	std::string title = "COMPND    CLUSTER";
	std::string autho = "AUTHOR    GENERATED BY RICHTOOL";
	std::string fin1  = "MASTER        0    0    0    0    0    0    0    0";
	std::stringstream    fins;
	fins << std::setw(5) << ((int)round(px.size())) << std::setw(5) << 0;
	fins << std::setw(5) << ((int)round(px.size())) << std::setw(5) << 0;

	if( w->size == px.size() ){
		outp_coord << title << std::endl;
		outp_coord << autho << std::endl;
		for(int i=0;i<px.size();i++){
			outp_coord << std::right<< std::setw(6) << "HETATM" << std::setw(5) << i+1 << std::setw(3)
				<< px[i].first	<< std::setw(6) << "LIG" << std::setw(2) << ( (char) (65+round(gsl_vector_get(w,i))) )
				<< std::setw(4) << 1 << "    "	<< std::setprecision(4) 
				<< std::setw(8) << gsl_vector_get(px[i].second,XX) 
				<< std::setw(8) << gsl_vector_get(px[i].second,YY)
				<< std::setw(8) << gsl_vector_get(px[i].second,ZZ) << std::setprecision(3)
				<< std::setw(6) << 1.00 << std::setw(6) << 0.00 << std::setw(12) << px[i].first << std::endl;
		}
	}
	outp_coord << fin1 << fins.str() << std::endl;
	outp_coord << "END" << std::endl;

	outp_coord.close();
}

void 
fileIO::output_pdb( std::string filename, gsl_matrix *M, gsl_vector *w ) {
//NOTE:: pymol can color the clusters with CMD:: util.cbc
	const char *c_filename = filename.c_str();
	std::ofstream outp_coord;
	outp_coord.open(c_filename);
	std::string title = "COMPND    CLUSTER";
	std::string autho = "AUTHOR    GENERATED BY RICHTOOL";
	std::string fin1  = "MASTER        0    0    0    0    0    0    0    0";
	std::stringstream    fins;
	double xX,yY,zZ;

	int n = M->size1, m = M->size2;
	int L = n==DIM?m:n;

	fins << std::setw(5) << L << std::setw(5) << 0;
	fins << std::setw(5) << L << std::setw(5) << 0;

	if( ((int)w->size) == L ){
		outp_coord << title << std::endl;
		outp_coord << autho << std::endl;
		for(int i=0; i<L; i++ ){

			xX=(n==DIM)?(gsl_matrix_get(M,XX,i)):(gsl_matrix_get(M,i,XX));
			yY=(n==DIM)?(gsl_matrix_get(M,YY,i)):(gsl_matrix_get(M,i,YY));
			zZ=(n==DIM)?(gsl_matrix_get(M,ZZ,i)):(gsl_matrix_get(M,i,ZZ));

			outp_coord << std::right<< std::setw(6) << "HETATM" << std::setw(5) << i+1 << std::setw(3)
				<< 'H'	<< std::setw(6) << "LIG" << std::setw(2) << ( (char) (65+round(gsl_vector_get(w,i))) )
				<< std::setw(4) << 1  << "    "	<< std::setprecision(4) 
				<< std::setw(8) << xX << std::setw(8) << yY << std::setw(8) << zZ 
				<< std::setprecision(3) << std::setw(6) << 1.00 
				<< std::setw(6) << 0.00 << std::setw(12) << 'H' << std::endl;
		}
	}
	outp_coord << fin1 << fins.str() << std::endl;
	outp_coord << "END" << std::endl;
	outp_coord.close();
}

void 
fileIO::output_pdb( std::string filename, gsl_matrix *M, std::vector<std::string > anames ) {
//NOTE:: pymol can color the clusters with CMD:: util.cbc
	const char *c_filename = filename.c_str();
	std::ofstream outp_coord;
	outp_coord.open(c_filename);
	std::string title = "COMPND    CLUSTER";
	std::string autho = "AUTHOR    GENERATED BY RICHTOOL";
	std::string fin1  = "MASTER        0    0    0    0    0    0    0    0";
	std::stringstream    fins;
	double xX,yY,zZ;

	int n = M->size1, m = M->size2;
	int L = n==DIM?m:n;

	fins << std::setw(5) << L << std::setw(5) << 0;
	fins << std::setw(5) << L << std::setw(5) << 0;

	if( (int)(anames.size()) == L ){
		outp_coord << title << std::endl;
		outp_coord << autho << std::endl;
		for(int i=0; i<L; i++ ){

			xX=(n==DIM)?(gsl_matrix_get(M,XX,i)):(gsl_matrix_get(M,i,XX));
			yY=(n==DIM)?(gsl_matrix_get(M,YY,i)):(gsl_matrix_get(M,i,YY));
			zZ=(n==DIM)?(gsl_matrix_get(M,ZZ,i)):(gsl_matrix_get(M,i,ZZ));

			outp_coord << std::right<< std::setw(6) << "HETATM" << std::setw(5) << i+1 << std::setw(3)
				<< anames[i]	<< std::setw(6) << "LIG" << std::setw(2) << 'A'
				<< std::setw(4) << 1  << "    "	<< std::setprecision(4) 
				<< std::setw(8) << xX << std::setw(8) << yY << std::setw(8) << zZ 
				<< std::setprecision(3) << std::setw(6) << 1.00 
				<< std::setw(6) << 0.00 << std::setw(12) << anames[i] << std::endl;
		}
	}
	outp_coord << fin1 << fins.str() << std::endl;
	outp_coord << "END" << std::endl;
	outp_coord.close();
}



