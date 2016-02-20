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

//	CRYSTAL CELL STUFF
	int RC;
	if (mmdb.isSpaceGroup() && !(mmdb_N.isSpaceGroup()) )  {
		// space group name is for demonstration only
		RC = mmdb_N.SetSpaceGroup ( mmdb.GetSpaceGroup () );
		if (RC!=SYMOP_Ok)   {
			switch (RC)  {
				case SYMOP_NoLibFile :
					printf ( " **** error: can't find symop.lib\n" );
					break;
				case SYMOP_UnknownSpaceGroup :
					printf ( " **** error: attempt to set up unknown space group\n" );
					break;
				case SYMOP_NoSymOps :
					printf ( " **** error: no symmetry operations found\n" );
					break;
				default :
					printf ( " **** error: unknown return code from "
					"CMMDBManager::SetSpaceGroup()\n" );
				}
			exit(2);
		}
		std::cout << "INFO::ASSIGNED SYM. INFORMATION"<< std::endl;
	}
	if (mmdb.isCrystInfo())  {
		// numerical values are for demonstration only
		double cell[8];
		int OC[0];
    		mmdb.GetCell ( cell[0], cell[1], cell[2], cell[3] , cell[4] , cell[5] , cell[6] , OC[0] );
		mmdb_N.SetCell ( cell[0], cell[1], cell[2], cell[3] , cell[4] , cell[5] , OC[0] );
		std::cout << "INFO::ASSIGNED CELL INFORMATION"<< std::endl;
	}
	RC = mmdb_N.CrystReady();
	if (RC>0)  {
		if (RC & CRRDY_NotPrecise)
			std::cout << "WARNING:: 1 \n" ;
		if (RC & CRRDY_isTranslation)
			std::cout << "WARNING:: 2 \n" ;
		if (RC & CRRDY_NoOrthCode)
			std::cout << "WARNING:: 3 \n" ;
	}
	RC = mmdb_N.GenerateSymMates ( NULL );
	if (RC>0)  {
		std::cout << "WARNING:: " << RC <<" \n" ;
	}
//	MOD MANAGER HAS CELL AND SYM

	std::cout << "INFO:: HAVE " << nModels << " MODELS" << std::endl;
	int nErr=0;

	rich::mmdb_helper mmhelp;

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
	int nb = NBINS_IO;
	for ( imod=1 ; imod<=nModels ; imod++ ) {
		nChains = mmdb.GetNumberOfChains( imod ); 
		for ( icha = 0 ; icha < nChains ; icha++ ) { 
			nResidues = mmdb.GetNumberOfResidues( imod , icha ); 
			rich::particles residue_atoms;
			for ( ires = 0 ; ires < nResidues ; ires++ ) { 

				rich::residue_helper rh;
				double rc = rh.analyze_stage1( imod, icha, ires, &mmdb, nb , &residue_atoms );
				double zc = rh.calc_OS();			// CALCULATE ORTHONORMAL SYSTEM
				gsl_matrix *OS	= gsl_matrix_calloc( DIM, DIM );
				rh.copyOS(OS);
				if(rh.skip())
					continue;
				rh.calc_proj( density, &theta , nb );
				gsl_vector *v0  = gsl_vector_calloc(DIM);
				rh.copyv0(v0);

//			HERE WE ARE AT SPECIFIC PROJECTIONS PROBLEM
//			CALULATE PROJECTION
				if( ires == 40 && icha==0 ) {	// TESTCASE
				//if( 1 ) {
					if( verbose ) {			// REALLY A DIAGNOSTICS TOOL, VERBOSE
						rich::mat_io mIO;
						mIO.write_vdbl2dat( theta, "testTheta.dat" );
					}
					std::sort( theta.begin(), theta.end(), compare_pid );
					double val = 100.0, alimit = 7E-2;
					std::vector<double> v_ang;
					int Ith = theta.size();
					do {
						val		= theta[Ith-1].second		/ sqrt(float(NBINS_IO));
						double ang	= theta[Ith-1].first  * 360.0	/ (float(NBINS_IO));
						if( verbose )
							std::cout << "INFO::MAX > " << Ith << " " << ang << " " << val << std::endl;
						if( val > alimit )
							v_ang.push_back(ang);
						Ith--;
					} while( val > alimit && Ith>=0 ) ;

					rich::fileIO	fIO;
					gsl_matrix_get_row( nh, OS, 0 );

					if( v_ang.size() < 1)
						continue;

					for(int i_ang=0; i_ang < v_ang.size() ; i_ang++) {
						model = newCModel(	);
						model->Copy( model_T[0] );

						rich::quaternion q;
						double fi = (i_ang==0)?( v_ang[i_ang] ):(v_ang[i_ang]-v_ang[i_ang-1]);
						q.assign_quaterion( nh , v_ang[i_ang] * M_PI/180.0 );
						q.rotate_particles( residue_atoms , v0 );

						if( mmhelp.check_clash( NM, icha, ires, &mmdb_N, residue_atoms, 1.0 ) > 1 ) {
							if(verbose)
								std::cout << "INFO::WE HAVE CLASH" << std::endl;
						} else {
							NM++;
							mmdb_N.AddModel( model  ); 
							if( mmhelp.update_residue( NM ,icha ,ires, &mmdb_N, residue_atoms ) ) {
								std::cout << "::ERROR::" << std::endl;
								fatal();
							}
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
			}
		}
	}

	mmdb_N.GetModelTable( model_T, nModels );
	std::cout << "INFO::GENERATED "<<nModels<<" MODELS"<< NM << std::endl;
	mmdb_N.FinishStructEdit();
	mmdb_N.WritePDBASCII( "multistate.pdb" );
	std::cout <<"INFO>> " << 0/1E-10 << std::endl;
	return 0;
}

