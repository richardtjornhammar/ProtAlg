//// C/C++ STUFF
#include <iostream>
#include <string>
#include <fstream>
#include <vector>
#include <sstream>
#include <math.h>
#include <cstdlib>
#include <cstring>
#include <random>

//	MMDB	STUFF
#include "mmdb/mmdb_manager.h"

//	GSL	STUFF
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_permutation.h>

//	CLIPPER	STUFF
#include <clipper/clipper.h>
#include <clipper/clipper-ccp4.h>
#include <clipper/clipper-contrib.h>

//	OWN STUFF
#include "types.h"
#include "argpa.h"	
#include "iofun.h"
#include "mdiag.h"

// rm src/*~
// COMPILE:: g++ -I/usr/local/include/clipper -std=c++11 src/* -L/usr/local/lib -lmmdb -lclipper-core -lclipper-contrib -lgsl -lblas -o rich_dyn
// COMPILE:: export LD_LIBRARY_PATH=/usr/local/lib 
// COMPILE:: g++ -std=c++11 src/* -lmmdb -lclipper-core -lclipper-ccp4 -lclipper-contrib -lgsl -lblas -o rich_dyn
// COMPILE:: export LD_LIBRARY_PATH=/home/richard/autobuild/Linux-Carmack-pre-release-gtk2-noguile/lib && g++ -std=c++11 src/* -lmmdb -lclipper-core -lclipper-ccp4 -lclipper-contrib -lgsl -lblas -o rich_dyn

void
fatal(void) {
	std::string	author("Richard Tjörnhammar (e-mail: richard.tjornhammar@gmail.com)");
	std::cout << "INFO:: PROGRAM FAILED" << std::endl;
	std::cout << "PLEASE CONTACT " << author << "\nWITH ANY QUESTIONS REGARDING THE FUNCTIONALITY" << std::endl;
	exit(1);
}

int main ( int argc, char ** argv ) {
	int verbose = 0;

//	PREPARING ARGPARSER
	std::pair < int, char ** >	margs;
	margs.first	= argc; 
	margs.second	= argv;

	std::vector < std::pair < int, std::pair<std::string, std::string > > >	opts_n_defaults;
	std::pair	< int, std::pair<std::string,std::string > > 		vipss;
	vipss.first=1;	vipss.second.first = "-mfile"; 
			vipss.second.second = "m_default.pdb";
	opts_n_defaults.push_back(vipss);
	vipss.first=0;	vipss.second.first = "-dfile"; 
			vipss.second.second = "d_default.mtz";
	opts_n_defaults.push_back(vipss);
	vipss.first=0;	vipss.second.first = "-model";
			vipss.second.second = "10";
	opts_n_defaults.push_back(vipss);
	vipss.first=0;	vipss.second.first = "-verbose";
			vipss.second.second = "0";
	opts_n_defaults.push_back(vipss);

//	EXECUTING ARGPARSE
	rich::arg_parser parser;
	std::pair<int, std::vector< std::string > > run_args = parser.parse( margs, opts_n_defaults );
	if(run_args.first) { // FAILED
		fatal();
	}

//	HERE WE SHOULD SET THE INPUT PARAMETERS
	verbose		= atoi(run_args.second[3].c_str());
	std::cout << "INFO:: SETTING VERBOSE LEVEL "<< verbose << std::endl;

	rich::map_manager mmap;
	int rv = mmap.assign_map( run_args.second[1] );
	clipper::Xmap<float> density = mmap.get_map();

//	THIS TOOL CURRENTLY ONLY FOR PDB
	CMMDBManager   mmdb;
	int itype, otype;
	itype = MMDB_FILE_PDB;

//	INITIALIZE MMDB
	InitMatType();

// 	SETTING SLACK FLAGS
	mmdb.SetFlag(	MMDBF_IgnoreDuplSeqNum	| 
                    	MMDBF_IgnoreBlankLines	|
                        MMDBF_IgnoreRemarks	|
                        MMDBF_IgnoreHash    	|
                        MMDBF_IgnoreNonCoorPDBErrors |
                        MMDBF_IgnoreSegID	|
                        MMDBF_IgnoreElement	|
                        MMDBF_IgnoreCharge	);

// 	READ COORDINATE FROM USING THE SUPPLIED NAME
	std::cout << "INFO:: WILL READ FILE " << run_args.second[0] << std::endl;
	int bytes = run_args.second[0].length();
	char * cpstr	= new char[bytes + 1];
	std::strcpy (cpstr, run_args.second[0].c_str() );
	int rval	= mmdb.ReadPDBASCII( cpstr );	//	READ PDB
	if(rval!=0) {
		fatal(  );
	}

//	START TO WORK WITH THE DATA (TABLES)
	PPCModel	model_T;
	PPCChain	chain_T;
	PPCResidue	resid_T;
	PPCAtom		atoms_T,atom_T2;
	int		nModels;
	int 		nChains,nResidues,nAtoms;
	int		imod,icha,ires,iat;

	mmdb.GetModelTable( model_T, nModels );
	std::cout << "INFO:: HAVE " << nModels << " MODELS" << std::endl;
	int nErr=0;

	rich::mmdb_helper mmhelp;

//	AMINO VECTORS
	gsl_vector *n1 = gsl_vector_calloc(DIM);
	gsl_vector *n2 = gsl_vector_calloc(DIM);
	gsl_vector *c1 = gsl_vector_calloc(DIM);
	gsl_vector *c2 = gsl_vector_calloc(DIM);
	gsl_vector *vt = gsl_vector_calloc(DIM);

//	ORTHONORMAL SYSTEM STORAGE
	gsl_vector *nh = gsl_vector_calloc(DIM);
	gsl_vector *ph = gsl_vector_calloc(DIM);
	gsl_vector *qh = gsl_vector_calloc(DIM);

	for ( imod=1 ; imod<=nModels ; imod++ ) {
		nChains = mmdb.GetNumberOfChains( imod ); 
		for ( icha = 0 ; icha < nChains ; icha++ ) { 
			nResidues = mmdb.GetNumberOfResidues( imod , icha ); 

			for ( ires = 0 ; ires < nResidues ; ires++ ) { 
				mmdb.GetAtomTable    ( imod ,icha ,ires , atoms_T, nAtoms );
				gsl_matrix *A	= gsl_matrix_calloc(nAtoms,DIM);
				gsl_matrix *V	= gsl_matrix_calloc(DIM,DIM);
				gsl_matrix *OS	= gsl_matrix_calloc(DIM,DIM);
				gsl_vector *S	= gsl_vector_calloc(DIM);
				gsl_vector *wrk	= gsl_vector_calloc(DIM);

				for (iat = 0; iat<nAtoms ; iat++ ) {
					char a_inf[256];
					atoms_T[iat]->GetAtomID(a_inf);
					std::string atype = mmhelp.atom_type(a_inf);
					std::string etype = mmhelp.atom_symb(a_inf);
					std::cout << "ATOM:: " << atype << " " << etype << std::endl;
					gsl_vector_set(vt,XX,atoms_T[iat]->x);
					gsl_vector_set(vt,YY,atoms_T[iat]->y);
					gsl_vector_set(vt,ZZ,atoms_T[iat]->z);
					if( atype=="CA" )
						gsl_vector_memcpy(c1,vt);
					if( atype=="CB" )
						gsl_vector_memcpy(c2,vt);
					if( atype=="N" && etype=="N" )
						gsl_vector_memcpy(n1,vt);	
					if( atype=="C" && etype=="C" )
						gsl_vector_memcpy(n2,vt);
					gsl_matrix_set_row(A,iat,vt);
				}
//			CALCULATE ORTHONORMAL SYSTEM
				rich::math_helper mah;
				gsl_vector_sub(n2,n1);
				gsl_vector_sub(c2,c1);
				gsl_vector_memcpy(nh,n2);
				gsl_vector_memcpy(ph,c2);
				double nl = gsl_blas_dnrm2(nh);
				double pl = gsl_blas_dnrm2(ph);
				gsl_vector_scale(nh,1.0/nl); gsl_vector_scale(ph,1.0/pl);
				mah.gsl_cross3D(nh,ph,qh);
				double ql = gsl_blas_dnrm2(qh);			
				gsl_vector_scale(qh,1.0/ql);

//			ASSIGN ORTHONORMAL SYSTEM
				gsl_matrix_set_row(OS,0,nh);
				gsl_matrix_set_row(OS,1,ph);
				gsl_matrix_set_row(OS,2,qh);

//			CALCULATE LIMITS
				gsl_linalg_SV_decomp(A,V,S,wrk);
				double rc = sqrt(gsl_vector_get(S,0))*0.5;
				double zc = nl;

//			CALULATE PROJECTION
				rich::calc_map cmap;
				int nb = cmap.get_nbins();
				gsl_matrix *P = gsl_matrix_calloc(nb,nb);
				cmap.proj(	P , density ,
						OS , n1 , rc , zc );
				if(verbose) {
					rich::tensorIO tIO;
					tIO.output_matrix(OS); tIO.output_matrix(V);
				}

				gsl_matrix_free( P );
				gsl_matrix_free( A );
				gsl_matrix_free( V );
				gsl_matrix_free( OS );
				gsl_vector_free( S );
				gsl_vector_free( wrk );
			}
		}
	}

	return 0;
}

