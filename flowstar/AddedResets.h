#include "Hybrid.h"

using namespace flowstar;

void div_reset(TaylorModel &tmReset, const Interval intC, const int varStoreInd, const int varDenInd, const int numVars);

void sec_reset(TaylorModel &tmReset, const Interval intC, const int varStoreInd, const int varInputInd, const int numVars, const TaylorModelVec tmvAggregation, const std::vector<Interval> doAggregation);

void arc_reset(TaylorModel &tmReset, const Interval intC, const int varStoreInd, const int varInputInd, const int numVars);

void tan_reset(TaylorModel &tmReset, const Interval intC, const int varStoreInd, const int varInputInd, const int numVars);

void cos_reset(TaylorModel &tmReset, const Interval intC, const int varStoreInd, const int varInputInd, const int numVars);

void sin_reset(TaylorModel &tmReset, const Interval intC, const int varStoreInd, const int varInputInd, const int numVars);

void sqrt_reset(TaylorModel &tmReset, const Interval intC, const int varStoreInd, const int varInputInd, const int numVars);
