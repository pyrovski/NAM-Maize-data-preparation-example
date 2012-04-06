#include <fstream>
#include <iostream>
#include <stdlib.h>
#include <getopt.h>

using namespace std;

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
  char *imputed_filename = 0, 
    *map_filename = 0,
    *fastphase_filename = 0,
    *pheno_filename = 0;
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
  return 0;
}
