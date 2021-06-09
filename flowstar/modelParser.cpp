/*---
  Flow*: A Verification Tool for Cyber-Physical Systems.
  Authors: Xin Chen, Sriram Sankaranarayanan, and Erika Abraham.
  Email: Xin Chen <chenxin415@gmail.com> if you have questions or comments.
  
  The code is released as is under the GNU General Public License (GPL).
---*/

#include "unistd.h"
#include "modelParser.h"
#include "DNNResets.h"
#include "DNN.h"
#include "TaylorModel.h"
#include "NNTaylorModel.h"

extern int yyparse();

//std::string dnn::DNN_Filename;
extern std::vector<std::string> dnn::DNN_Filenames;
extern bool dnn::dnn_initialized;
extern std::vector<std::string> dnn::initialConds;
//extern bool dnn::load_reset;
//extern float dnn::totalNumBranches;
//extern float dnn::dnn_runtime;

bool dnn::plottingEnabled = false;
bool dnn::dumpingEnabled = false;

int main(int argc, char **argv)
{
	int num_threads = 1;
	int c;
	while( (c = getopt(argc, argv, "pdt:")) != -1 ) {
		switch (c) {
			case 'p':
				dnn::plottingEnabled = true;
				break;
			case 'd':
				dnn::dumpingEnabled = true;
				break;
			case 't':
				num_threads = atoi(optarg);
				break;
			case '?':
				if (optopt == 't') {
					fprintf (stderr, "Option -%c requires an argument.\n", optopt);
				} else if (isprint (optopt))
					fprintf (stderr, "Unknown option `-%c'.\n", optopt);
				else
					fprintf (stderr,
							"Unknown option character `\\x%x'.\n",
							optopt);
        		return 1;
      		default:
        		abort ();
		}
	}

	while (optind < argc) {
	        //dnn::DNN_Filename = argv[optind];

	        dnn::DNN_Filenames.push_back(argv[optind]);

		optind++;
	}

	if (dnn::DNN_Filenames.empty()) {
		fprintf( stderr, "Missing DNN filename. Expected usage: 'flowstar [-p, -d] <dnn filename>'\n" );
		//return(-1);
	}
	dnn::dnn_initialized = false;

	NNTaylorModelVec::nnpool = new ctpl::thread_pool(num_threads);
	
	yyparse();

	printf("total branches: %d\n", dnn::totalNumBranches);
	printf("dnn runtime: %f\n", dnn::dnn_runtime);
	//printf("DNN filename: %s\n", dnn::DNN_Filenames[0].c_str());

	printf("\nInitial conditions:\n");
	for(int i = 0; i < dnn::initialConds.size(); i++){
	        printf("%s\n", dnn::initialConds[i].c_str());
	}

	NNTaylorModelVec::nnpool->stop();
	delete NNTaylorModelVec::nnpool;	

	return 0;
}



