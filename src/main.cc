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
// COMPILE:: export LD_LIBRARY_PATH=/home/richard/autobuild/Linux-Carmack-pre-release-gtk2-noguile/lib
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
	vipss.first=1;	vipss.second.first = "-dfile"; 
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
	CMMDBManager	*mmdb_p;
	mmdb_p=&mmdb;
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
	int		nModels;
	int 		nChains,nResidues,nAtoms,nAtoms2;
	int		imod,icha,ires,iat;

	if( verbose )
		std::cout << "INFO:: BUILD HELP" << std::endl;
	rich::mmdb_helper mmhelp(&mmdb) ;
	if( verbose )
		std::cout << "INFO:: GEN SYM" << std::endl;
	mmhelp.generate_sym( );
	if( verbose )
		std::cout << "INFO:: GET MOL" << std::endl;
	mmhelp.get_mol( );
	if( verbose )
		std::cout << "INFO:: HAS MAN" << std::endl;
	mmhelp.has_manager( ); 
	if( verbose )
		std::cout << "INFO:: HAS UNI" << std::endl;
	mmhelp.has_unitcell( ); 
	if( verbose )
		std::cout << "INFO:: HAS CELL" << std::endl;
	mmhelp.has_cell( ); 
	if( verbose )
		std::cout << "INFO:: HAS CHEEZEBURGER" << std::endl;
	mmhelp.has_symmetry( );
	if( verbose )
		std::cout << "INFO:: SANE" << std::endl;
	mmhelp.remove_sym( );

	//mmdb.GetModelTable	( model_T, nModels );
	mmhelp.get_mol()->GetModelTable( model_T, nModels );
	if( verbose )
		std::cout << "INFO:: USED STRUCT" << std::endl;

	model = newCModel	(	);
	model -> Copy		( model_T[0] );
	mmdb_N.AddModel		(  model );
	int NM=1;

	// XTAL HEADER
	mmhelp.copy_xtal( &mmdb , &mmdb_N );

	if( verbose )
		std::cout << "INFO:: HAVE " << nModels << " MODELS" << std::endl;
	int nErr=0;

	// ORTHONORMAL SYSTEM VECTORS
	gsl_vector *nh = gsl_vector_calloc(DIM);
	gsl_vector *ph = gsl_vector_calloc(DIM);
	gsl_vector *qh = gsl_vector_calloc(DIM);

	std::vector< std::pair<int,double> > theta, chi;
	int Ntheta = NBINS_IO;

	for ( int i=0 ; i<Ntheta ; i++ ) {
		std::pair< int, double > pid;
		pid.first=i; pid.second=0.0;
		theta.push_back(pid);
		chi.push_back(pid);
	}

	int nb		= NBINS_IO;
	double TOL	= 30E-2;
	std::cout << "INFO:: HERE" << std::endl;

	for ( imod=1 ; imod<=nModels ; imod++ ) {
		nChains = mmdb.GetNumberOfChains( imod );
		if( verbose )
			std::cout << "INFO:: CHAIN LOOP " << nChains << std::endl; 
		for ( icha = 0 ; icha < nChains ; icha++ ) { 
			nResidues = mmdb.GetNumberOfResidues( imod , icha ); 
			rich::particles residue_atoms;
			if( verbose )
				std::cout << "INFO:: RES LOOP" << std::endl;
			for ( ires = 0 ; ires < nResidues ; ires++ ) { 

				rich::residue_helper rh;
				double rc = rh.analyze_stage1( imod, icha, ires, &mmdb, nb , &residue_atoms );
				double zc = rh.calc_OS();			// CALCULATE ORTHONORMAL SYSTEM
				if( rh.skip() )
					continue;
				if( (ires == 40 || ires==88) && icha==0 && verbose == 2)
					rh.calc_proj( density, mmap.get_gsa(), &theta , nb , 0 , ires );
				else
					rh.calc_proj( density, mmap.get_gsa(), &theta , nb , 0 );
				gsl_vector *v0  = gsl_vector_calloc(DIM);
				rh.copyv0(v0);

//			HERE WE ARE AT SPECIFIC PROJECTIONS PROBLEM
//			CALCULATE PROJECTION
				if( (ires == 40 || ires==88) && icha==0 ) {	// TESTCASE
					if( verbose ) {			
						rich::mat_io mIO;
						mIO.write_vdbl2dat( theta, "res"+std::to_string(ires) +
							 "NB" + std::to_string(NBINS_IO) + "Theta.dat" );
					}

					std::vector<double> v_ang = rh.prune_angles( &theta , TOL, nb ); 

					if( v_ang.size() < 1)
						continue;

					rich::fileIO	fIO;
					gsl_matrix *OS	= gsl_matrix_calloc( DIM, DIM );
					rh.copyOS(OS);
					gsl_matrix_get_row( nh, OS, 0 );
					rich::particle_helper parth;

					TOL=(ires==88)?(TOL*0.1):(TOL);

					for(int i_ang=0; i_ang < v_ang.size() ; i_ang++) {

						rich::particles rotres_atoms0;
						rich::quaternion q;
						q.assign_quaterion( nh , v_ang[i_ang] * M_PI/180.0 );
						rotres_atoms0 = parth.particles_memcpy( residue_atoms );
						q.rotate_particles( rotres_atoms0 , v0 );

						if( rh.do2nd() ) {

							gsl_vector *v00  = gsl_vector_calloc(DIM);
							rich::particles rotres_atoms1;
							rotres_atoms1 = parth.particles_memcpy( rotres_atoms0 ); 
							double zc_rr1 = rh.calc_O1( rotres_atoms1 );
							rh.copyv0(v00);

							rh.calc_proj( density, mmap.get_gsa(), &chi , nb , 1 , ires);

							std::vector<double> v_chi = rh.prune_angles( &chi , TOL, nb );
							gsl_matrix *O1	= gsl_matrix_calloc( DIM, DIM );
							rh.copyO1(O1);
							gsl_matrix_get_row( ph, O1, 0 );
							double chi0 = rh.calc_fi(rotres_atoms1,1); // HERE

							for(int i_chi=0; i_chi < v_chi.size() ; i_chi++) {
								rich::quaternion qq;
								qq.assign_quaterion( ph , (v_chi[i_chi]-chi0) * M_PI/180.0 );
								std::vector<bool> mask = rh.get_mask(0);
								if( mask.size()!=rotres_atoms1.size() )
									if(verbose)
										std::cout << "INFO:: WE HAVE A BAD MASK" << std::endl;		
								qq.rotate_particles( rotres_atoms1 , v00 , mask );
								if( mmhelp.check_clash_sym( 1 , icha, ires, &mmdb, rotres_atoms1, 1.2 )  ) {
									if(verbose)
										std::cout << "INFO:: WE HAVE CLASH" << std::endl;
								} else {
									model = newCModel(	);
									model->Copy( model_T[0] );
									NM++;
									mmdb_N.AddModel( model  );
									std::cout << "::GONNA UPDATE::" << std::endl;
									if( mmhelp.update_residue( NM ,icha ,ires, &mmdb_N, rotres_atoms1 ) ) {
										std::cout << "::ERROR::" << std::endl;
										fatal();
									}
									std::cout << "::DID IT::" << std::endl;
									if( verbose == 2 ) {
										char a_inf[256];
										PPCAtom		atoms_T, atoms_T2;
										mmdb_N.GetAtomTable    ( NM , icha ,ires , atoms_T, nAtoms );
										atoms_T[0]->GetAtomID(a_inf);
										std::cout << "INFO:1: " << a_inf << std::endl; 
									}
								}
							}
						}

						if( mmhelp.check_clash_sym( NM, icha, ires, &mmdb_N, rotres_atoms0, 1.0 )  ) {
							if(verbose)
								std::cout << "INFO:: WE HAVE CLASH" << std::endl;
						} else {
							model = newCModel(	);
							model->Copy( model_T[0] );
							NM++;
							mmdb_N.AddModel( model  ); 
							if( mmhelp.update_residue( NM , icha , ires, &mmdb_N, rotres_atoms0 ) ) {
								std::cout << "::ERROR::" << std::endl;
								fatal();
							}
							if( verbose == 2 ) {
								char a_inf[256];
								PPCAtom	 atoms_T,atoms_T2;
								mmdb_N.GetAtomTable    ( NM , icha ,ires , atoms_T, nAtoms );
								atoms_T[0]->GetAtomID(a_inf);
								std::cout << "INFO:0: " << a_inf << std::endl; 
							}	
						}
						
						if( verbose == 2 ) {
							fIO.output_pdb( "rotres" + std::to_string(i_ang) + ".pdb" , residue_atoms );
						}
					}

					if ( verbose == 2 ) { 
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
			}
		}
	}

	mmdb_N.GetModelTable( model_T, nModels );
	std::cout << "INFO:: GENERATED " << nModels << " MODELS (" << NM << ")" << std::endl;
	mmdb_N.FinishStructEdit();
	if(nModels>0)
		mmdb_N.WritePDBASCII( "state.pdb" );

	return 0;
}

