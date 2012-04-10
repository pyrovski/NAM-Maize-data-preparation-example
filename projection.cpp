/*
check out lower_bound for binary search
 */
#include <fstream>
#include <iostream>
#include <stdlib.h>
#include <getopt.h>

using namespace std;

char *imputed_filename = 0, 
  *map_filename = 0,
  *fastphase_filename = 0,
  *pheno_filename = 0;

ifstream imputed_file, map_file, fastphase_file, pheno_file;

/*!
  
 */
 

/*!
  Project one row of the output matrix.

  Update non-zero entries of row @row in the output matrix.

  Output matrix rows are expected to be quite long, so 
  vector-processing an entire row is unwise.  
  Perhaps it would be better to process the row in short chunks.

  Is there a BLAS call for element-wise subtraction?
  gsl_vector_add_constant
 */
int projectRow(){

  /*!
    for each non-zero element (index j) of the output matrix, 
    
    start with the rightmost element of newmap(:, 4),
    and find elements < fasts(j).

    For the first non-matching element, ...
      If no such element is found, ...
    If all elements match, ...
    
   */
}

void usage(const char * name){
  cerr << name << " --imputed <> --map <> --fastphase <> --pheno <> [--verbose [<>]]" << endl;
}

int parse_inputs(int argc, char ** argv){
  int status, longindex;
  struct option options[] = {
    {"imputed", 1, 0, 0},
    {"map", 1, 0, 0},
    {"fastphase", 1, 0, 0},
    {"pheno", 1, 0, 0},
    {"verbose", 2, 0, 0},
    {0,0,0,0}
  };

  int exitVal = 0;

  while((status = getopt_long(argc, argv, "h", options, &longindex)) != -1){
    switch(status){
    case 0:
      //long options
      switch(longindex){
      case 0:
	imputed_filename = optarg;
	break;
      case 1:
	map_filename = optarg;
	break;
      case 2:
	fastphase_filename = optarg;
	break;
      case 3:
	pheno_filename = optarg;
	break;
      default:
	cerr << "unsupported long option: " << longindex << endl;
	exit(1);
      }
      break;
    default:
      cerr << "unsupported option: " << status << endl;
      exitVal = 1;
    case 'h':
      usage(argv[0]);
      exit(exitVal);
    }
  }

  if(!imputed_filename || !map_filename || 
     !fastphase_filename || !pheno_filename){
    cerr << "missing a filename!" << endl;
    exit(1);
  }
  return 0;
}

int open_inputs(){
  imputed_file.open(imputed_filename);
  map_file.open(map_filename);
  fastphase_file.open(fastphase_filename);
  pheno_file.open(pheno_filename);
  return 0;
}

/*!
  read imputed, map, fast, pheno files in their entirety
 */
int read_inputs(){
  
}

/*!
  Inputs: 
  - imputed marker
  - map
  - fastphase
  - phenotypes

  Outputs:
  - projected SNPs
 */
int main(int argc, char ** argv){

  parse_inputs(argc, argv);
  open_inputs();

  read_inputs();

  return 0;
}
