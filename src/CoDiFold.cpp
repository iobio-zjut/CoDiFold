#include<cstdlib>
#include<iostream>
#include<fstream>
#include<sstream>
#include<string>
#include<vector>
#include <algorithm>
#include <math.h>
#include<sys/stat.h>
#include<sys/types.h>
#include<iomanip>
#include<unistd.h>

using namespace std;


//sequence
string fasta;
//sequence route
string fasta_route;
//seqence length
int fasta_length;
//contact map route
string cmap1_route, cmap2_route;
//distance profile file route
string dp_route;
//native_pdb route
string native_route;
//fragment library route
string frag_lib_3_route;
string frag_lib_9_route; 

//check whether fastaFile is input
bool fFile;
//check whether contact1 and contact2 are input
bool cFile1,cFile2;
//check whether distance profile file is input
bool dpFile;
//check whether frag_library are input
bool f3File;
bool f9File;
//check whether ntive structure is input
bool nFile;


void getParameter(int argc, char ** argv, int & i);
void extractParameters(int argc, char ** argv);

void apply();


int main(int argc, char ** argv){
	cout << endl;
	cout << "###########################################################################" << endl;
	cout << "#                                CoDiFold                                 #" << endl;
	cout << "#                              Version: 1.0                               #" << endl;
	cout << "#          	a de novo protein structure prediction         		  #" << endl;
	cout << "#    	       by coupling contact with distance profile    		  #" << endl;
	cout << "#                         Copyright (C) Guijun Zhang                      #" << endl;
	cout << "#                     College of Information engineering                  #" << endl;
	cout << "#            Zhejiang university of technology, Hangzhou, China           #" << endl;
	cout << "###########################################################################" << endl;
	
	
	
	
	extractParameters(argc, argv);
	
	apply();
	
	
	cout << endl;
	cout << "###########################################################################" << endl;
	cout << "#                                  END                                    #" << endl;
	cout << "###########################################################################" << endl;
  
	return 0;
}


void extractParameters(int argc, char ** argv) {
    int i = 1;
    while (i < argc)
        getParameter(argc, argv, i);
    if (!fFile) {
        cout << endl;
        cout << "#  Error! seq.fasta File must be provided" << endl << endl;
        exit(0);
    }
    if (!cFile1) {
        cout << endl;
        cout << "#  Error! contact1 File (ResPRE) must be provided" << endl << endl;
        exit(0);
    }
    if (!cFile2) {
        cout << endl;
        cout << "#  Error! contact2 File (Raptor-X) must be provided" << endl << endl;
        exit(0);
    }
    if (!dpFile) {
        cout << endl;
        cout << "#  Error! distance profile File must be provided" << endl << endl;
        exit(0);
    }
    if (!f3File) {
        cout << endl;
        cout << "#  Error! 3-mer fragment library File must be provided" << endl << endl;
        exit(0);
    }
    if (!f9File) {
        cout << endl;
        cout << "#  Error! 9-mer fragment library File must be provided" << endl << endl;
        exit(0);
    }
}

void getParameter(int argc, char ** argv, int & i){
	string flag( argv[i] );
	if ( flag == "-f" ){
		if (argc < i + 2){
			cout << endl;
			cout << "#  Error! No fasta file provided" << endl << endl;
			cout << "   -f 			seq.fasta 			: protein sequence file  " << endl;
			cout << "   -c1 		c1-ResPRE   			: contact map file " << endl;
			cout << "   -c2 		c2-RaptorX   			: contact map file " << endl;
			cout << "   -dp 		distance profile   		: distance profile file " << endl;
			cout << "   -frag3	3mer_fragment_library  	: fragment library with fragment length 3" << endl;
			cout << "   -frag9	9mer_fragment_library  	: fragment library with fragment length 9" << endl;
			exit(0);
		}
		fasta_route = argv[++i];
		fFile = true;
	}
	else if( flag == "-c1" ){
		if (argc < i + 2){
			cout << endl;
			cout << "#  Error! No contact1(ResPRE) file provided" << endl << endl;
			cout << "   -f 			seq.fasta 			: protein sequence file  " << endl;
			cout << "   -c1 		c1-ResPRE   			: contact map file " << endl;
			cout << "   -c2 		c2-RaptorX   			: contact map file " << endl;
			cout << "   -dp 		distance profile   		: distance profile file " << endl;
			cout << "   -frag3	3mer_fragment_library  	: fragment library with fragment length 3" << endl;
			cout << "   -frag9	9mer_fragment_library  	: fragment library with fragment length 9" << endl;
			exit(0);
		}
		cmap1_route = argv[++i];
		cFile1 = true;
	}
	else if( flag == "-c2" ){
		if (argc < i + 2){
			cout << endl;
			cout << "#  Error! No contact2(RaptorX) file provided" << endl << endl;
			cout << "   -f 			seq.fasta 			: protein sequence file  " << endl;
			cout << "   -c1 		c1-ResPRE   			: contact map file " << endl;
			cout << "   -c2 		c2-RaptorX   			: contact map file " << endl;
			cout << "   -dp 		distance profile   		: distance profile file " << endl;
			cout << "   -frag3	3mer_fragment_library  	: fragment library with fragment length 3" << endl;
			cout << "   -frag9	9mer_fragment_library  	: fragment library with fragment length 9" << endl;
			exit(0);
		}
		cmap2_route = argv[++i];
		cFile2 = true;
	}
	else if( flag == "-dp" ){
		if (argc < i + 2){
			cout << endl;
			cout << "#  Error! No distance profile file provided" << endl << endl;
			cout << "   -f 			seq.fasta 			: protein sequence file  " << endl;
			cout << "   -c1 		c1-ResPRE   			: contact map file " << endl;
			cout << "   -c2 		c2-RaptorX   			: contact map file " << endl;
			cout << "   -dp 		distance profile   		: distance profile file " << endl;
			cout << "   -frag3	3mer_fragment_library  	: fragment library with fragment length 3" << endl;
			cout << "   -frag9	9mer_fragment_library  	: fragment library with fragment length 9" << endl;
			exit(0);
		}
		dp_route = argv[++i];
		dpFile = true;
	}
	else if( flag == "-frag3" ){
		if (argc < i + 2){
			cout << endl;
			cout << "#  Error! No 3mer_fragment_library file provided" << endl << endl;
			cout << "   -f 			seq.fasta 			: protein sequence file  " << endl;
			cout << "   -c1 		c1-ResPRE   			: contact map file " << endl;
			cout << "   -c2 		c2-RaptorX   			: contact map file " << endl;
			cout << "   -dp 		distance profile   		: distance profile file " << endl;
			cout << "   -frag3	3mer_fragment_library  	: fragment library with fragment length 3" << endl;
			cout << "   -frag9	9mer_fragment_library  	: fragment library with fragment length 9" << endl;
			exit(0);
		}
		frag_lib_3_route = argv[++i];
		f3File = true;
	}
	else if( flag == "-frag9" ){
		if (argc < i + 2){
			cout << endl;
			cout << "#  Error! No 9mer_fragment_library file provided" << endl << endl;
			cout << "   -f 			seq.fasta 			: protein sequence file  " << endl;
			cout << "   -c1 		c1-ResPRE   			: contact map file " << endl;
			cout << "   -c2 		c2-RaptorX   			: contact map file " << endl;
			cout << "   -dp 		distance profile   		: distance profile file " << endl;
			cout << "   -frag3	3mer_fragment_library  	: fragment library with fragment length 3" << endl;
			cout << "   -frag9	9mer_fragment_library  	: fragment library with fragment length 9" << endl;
			exit(0);
		}
		frag_lib_9_route = argv[++i];
		f9File = true;
	}
	++i;
}

void apply(){
	stringstream command;
	command << "~/rosetta_src_2018.33.60351_bundle/main/source/bin/AbinitioRelax.linuxgccrelease "
		<< "-in:file:fasta " << fasta_route << " "
		<< "-in:file:frag3 " << frag_lib_3_route << " "
		<< "-in:file:frag9 " << frag_lib_9_route << " "
		<< "-abinitio:increase_cycles " << 10 << " "
		<< "-nstruct " << 1 << " "
		<< "-out:pdb";
	string run_rosetta;
	getline(command, run_rosetta);
		
	system(run_rosetta.c_str());
	
	system("rm score.fsc default.out");	
	system("mv S_00000001.pdb lowestEtotal.pdb");
}

