//// C/C++ STUFF
#include <iostream>
#include <string>
#include <fstream>
#include <vector>
#include <algorithm>
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

bool compare_pid ( std::pair< int , double > pid1 ,
		   std::pair< int , double > pid2 ) {
	return( pid1.second < pid2.second ) ; 
};

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
	vipss.first=0;	vipss.second.first = "-nbins";
			vipss.second.second = "24";
	opts_n_defaults.push_back(vipss);
	vipss.first=0;	vipss.second.first = "-verbose";
			vipss.second.second = "0";
	opts_n_defaults.push_back(vipss);

//	EXECUTING ARGPARSE
	rich::arg_parser parser;
	std::pair<int, std::vector< std::string > > run_args = parser.parse( margs, opts_n_defaults );
	if( run_args.first ) { 
		fatal();
	}

//	HERE WE SHOULD SET THE INPUT PARAMETERS
	verbose		= atoi(run_args.second[3].c_str());
	std::cout << "INFO:: SETTING VERBOSE LEVEL "<< verbose << std::endl;
	int NBINS_IO 	= atoi(run_args.second[2].c_str());
	std::cout << "INFO:: SETTING NBINS " << NBINS_IO << std::endl;

	rich::map_manager mmap;
	int rv = mmap.assign_map( run_args.second[1] );
	clipper::Xmap<float> density = mmap.get_map();

//	THIS TOOL CURRENTLY ONLY FOR PDB
	CMMDBManager   mmdb, mmdb_N;
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
		fatal( );
	}

//	START TO WORK WITH THE DATA (TABLES)
	PPCModel	model_T, model_N;
	PCModel		model;
	PPCChain	chain_T;
	PPCResidue	resid_T;
	PPCAtom		atoms_T,atoms_T2;
	int		nModels;
	int 		nChains,nResidues,nAtoms,nAtoms2;
	int		imod,icha,ires,iat;

	mmdb.GetModelTable( model_T, nModels );
	model = newCModel(	);
	model->Copy( model_T[0] );
	mmdb_N.AddModel(  model );
	int NM=1;

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

	std::vector< std::pair<int,double> > theta;
	int Ntheta=NBINS_IO;
	for ( int i=0 ; i<Ntheta ; i++ ) {
		std::pair< int, double > pid;
		pid.first=i; pid.second=0.0;
		theta.push_back(pid);
	}

	for ( imod=1 ; imod<=nModels ; imod++ ) {
		nChains = mmdb.GetNumberOfChains( imod ); 
		for ( icha = 0 ; icha < nChains ; icha++ ) { 
			nResidues = mmdb.GetNumberOfResidues( imod , icha ); 
			rich::particles residue_atoms;
			for ( ires = 0 ; ires < nResidues ; ires++ ) { 
				residue_atoms.clear();
				mmdb.GetAtomTable    ( imod ,icha ,ires , atoms_T, nAtoms );
				if( ires<nResidues-1 ) {
					mmdb.GetAtomTable    ( imod ,icha ,ires+1 , atoms_T2, nAtoms2 );
					gsl_vector_set(n2,XX,atoms_T2[0]->x);
					gsl_vector_set(n2,YY,atoms_T2[0]->y);
					gsl_vector_set(n2,ZZ,atoms_T2[0]->z);
				}

				gsl_matrix *A	= gsl_matrix_calloc( nAtoms,	DIM);
				gsl_matrix *V	= gsl_matrix_calloc( DIM,	DIM);
				gsl_matrix *OS	= gsl_matrix_calloc( DIM,	DIM);
				gsl_vector *S	= gsl_vector_calloc( DIM );
				gsl_vector *wrk	= gsl_vector_calloc( DIM );

				gsl_vector *v0 = gsl_vector_calloc(DIM);

				for (iat = 0; iat<nAtoms ; iat++ ) {

					char a_inf[256];
					atoms_T[iat]->GetAtomID(a_inf);
					std::string atype = mmhelp.atom_type(a_inf);
					std::string etype = mmhelp.atom_symb(a_inf);

					gsl_vector_set(vt,XX,atoms_T[iat]->x);
					gsl_vector_set(vt,YY,atoms_T[iat]->y);
					gsl_vector_set(vt,ZZ,atoms_T[iat]->z);

					if(iat==0)
						gsl_vector_memcpy(v0,vt);

					rich::particle res_atom;
					res_atom.second = gsl_vector_alloc(DIM);
					res_atom.first  = etype; 
					gsl_vector_memcpy(res_atom.second,vt);
					residue_atoms.push_back(res_atom);

					if( atype=="CA" )
						gsl_vector_memcpy(c1,vt);
					if( atype=="CB" )
						gsl_vector_memcpy(c2,vt);
					if( atype=="N" && etype=="N" )
						gsl_vector_memcpy(n1,vt);
					if( ires>=nResidues-1 )
						if( atype=="C" && etype=="C" )
							gsl_vector_memcpy(n2,vt);

					gsl_matrix_set_row(A,iat,vt);
				}

//			CALCULATE ORTHONORMAL SYSTEM
				rich::math_helper mah;
				double zc = mah.gsl_calc_orth( n2, n1, c2, c1, OS );

//			CALCULATE LIMITS
				gsl_linalg_SV_decomp(A,V,S,wrk);
				double rc = sqrt(gsl_vector_get(S,0))*0.5;

//			CALULATE PROJECTION
				rich::calc_map cmap;
				int nb = cmap.set_nbins(NBINS_IO);
				gsl_matrix *P  = gsl_matrix_calloc(nb,nb);
				gsl_matrix *CN = gsl_matrix_calloc(nb,nb);

				if( cmap.proj(	P , CN , density ,
						OS , v0 , rc , zc,
						&theta ) ) {
					std::cout << "ERROR::FAILED" << std::endl;
					fatal();
				}

				if( ires == 40 && icha==0 ) {		// TESTCASE
				//if( 1 ) {
					if( verbose ) {			// REALLY A DIAGNOSTICS TOOL, VERBOSE
						rich::mat_io mIO;
						mIO.write_gsl2datn( P, CN, "testproj.dat"  );
						mIO.write_vdbl2dat( theta, "testTheta.dat" );
					}
					std::sort( theta.begin(), theta.end(), compare_pid );
					double val=100.0, alimit=0.5;
					std::vector<double> v_ang;
					do {
						val		= theta.back().second		/ (float(NBINS_IO));
						double ang	= theta.back().first  * 360.0	/ (float(NBINS_IO));
						if(verbose)
							std::cout << "INFO::MAX > " << ang << " " << val << std::endl;
						if( val > alimit )
							v_ang.push_back(ang);
						theta.pop_back();
					} while( val > alimit ) ;
					rich::fileIO	fIO;
					gsl_matrix_get_row( nh, OS, 0 );

					if( v_ang.size() < 1)
						continue;

					for(int i_ang=0; i_ang < v_ang.size() ; i_ang++) {
						model = newCModel(	);
						mmdb_N.AddModel( model  ); 
						NM++;
						model->Copy( model_T[0] );

						rich::quaternion q;
						double fi = (i_ang==0)?( v_ang[i_ang] ):(v_ang[i_ang]-v_ang[i_ang-1]);
						q.assign_quaterion( nh , v_ang[i_ang] * M_PI/180.0 );
						q.rotate_particles( residue_atoms , v0 );

						if( mmhelp.update_residue( NM ,icha ,ires, &mmdb_N, residue_atoms ) ) {
							std::cout << "::ERROR::" << std::endl;
							fatal();
						}

						if( verbose ) {
							fIO.output_pdb( "rotres" + std::to_string(i_ang) + ".pdb" , residue_atoms );
						}

					}
					if ( verbose ) { 
						std::vector<std::string> vs;
						vs.push_back("Ga"); 
						vs.push_back("In"); 
						vs.push_back("Sn");
						fIO.output_pdb("system.pdb", OS, vs);
						rich::tensorIO tIO;
						rich::quaternion q_test;
						gsl_vector *rx = gsl_vector_calloc(DIM);
						gsl_vector *ax = gsl_vector_calloc(DIM);
						gsl_vector_set(rx,XX,1);
						gsl_vector_set(ax,YY,1);
						double fi=M_PI*0.5;
						tIO.output_vector( ax );
						q_test.assign_quaterion( ax , fi );
						q_test.print();
						tIO.output_vector(rx);
						q_test.rotate_coord(rx);
						tIO.output_vector(rx);
					}
				}

				if(verbose) {
					rich::tensorIO tIO;
					tIO.output_matrix(OS); tIO.output_matrix(V);
					tIO.output_vector(nh);
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

	mmdb_N.GetModelTable( model_T, nModels );
	std::cout << "INFO::GENERATED "<<nModels<<" MODELS"<< NM << std::endl;
	mmdb_N.FinishStructEdit();
	mmdb_N.WritePDBASCII( "multistate.pdb" );

	return 0;
}

