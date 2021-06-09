#ifndef NEURAL_NETWORK__REACH_H_
#define NEURAL_NETWORK__REACH_H_

#include "Hybrid.h"

using namespace flowstar;

namespace dnn_reachability {

  void compute_dnn_reachability(Flowpipe &result, const TaylorModelVec tmvImage,
				const Flowpipe fpAggregation, const std::string modeName,
				const std::vector<std::string> stateVarNames, const Variables realStateVars, const bool bPrint);

}

#endif /* NEURAL_NETWORK__REACH_H_ */
