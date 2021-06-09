#ifndef NEURAL_NETWORK_H_
#define NEURAL_NETWORK_H_

#include "Hybrid.h"
#include "Continuous.h"
#include "Polynomial.h"
#include "NNTaylorModel.h"
#include <map>

using namespace flowstar;

namespace dnn {

enum activation {LINEAR = 0, SIGMOID = 1, SWISH = 2, RELU = 3, TANH = 4};
//extern std::string DNN_Filename;
extern std::vector<std::string> DNN_Filenames;

extern std::map<int, Continuous_Reachability_Setting> dnn_crs;
extern std::map<int, std::vector<iMatrix>> dnn_reset_weights;
extern std::map<int, std::vector<iMatrix>> dnn_reset_biases;
extern std::map<int, std::vector<activation>> dnn_activations;
extern std::map<int, std::vector<std::string>> augmentedVarNames;
extern std::vector<std::string> curAugmentedVarNames;
extern std::map<int, Variables> augmentedStateVars;
extern std::map<int, NNTaylorModelVec> activation_reset_map;
 
extern bool dnn_initialized;
extern std::map<int, int> branch_origin;
extern std::map<int, Flowpipe> saved_plant_states;
extern std::map<int, TaylorModelVec> saved_plant_tmv;
extern int totalNumBranches;
extern int curBranchId;
extern float dnn_runtime;
extern bool storedInitialConds;
extern std::vector<std::string> initialConds;

extern bool plottingEnabled;
extern bool dumpingEnabled;

void act_reset(NNTaylorModel &tmReset, const Interval intC, const int varInd, const int numVars, const bool tanh_act, const int order);
 
void sig_reset(TaylorModel &tmReset, const Interval intC, const int varInd, const int numVars);
  
void tanh_reset(TaylorModel &tmReset, const Interval intC, const int varInd, const int numVars);

void convert_TM_dimension(NNTaylorModel &tmNew, const NNTaylorModel &tmOld,
			  const int dimNew, const int varInd, const std::vector<std::string> & varNames,
			  const std::map<int, int> &indexMap=std::map<int,int>());

void get_state_to_f_map(std::map<int, int> &stateToF, std::map<int, int> &fToState,
			const std::vector<std::string> &stateVarNames, const TaylorModelVec &tmv);
 
void load_dnn(std::vector<iMatrix> &resets, std::vector<iMatrix> &reset_biases,
	      std::vector<dnn::activation> &activations, Variables &augmentedStateVars,
	      std::vector<std::string> &augmentedVarNames, const Variables &vars,
	      const std::vector<std::string> &stateVarNames, const std::string filename);

}

#endif /* NEURAL_NETWORK_H_ */
