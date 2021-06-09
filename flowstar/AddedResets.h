#include "Hybrid.h"

using namespace flowstar;

enum reset_types {SEC = 0, ARC = 1, COS = 2, TAN = 3, SIN = 4, DIV = 5};

void div_reset(TaylorModel &tmReset, const Interval intC, const int varStoreInd, const int varDenInd, const int numVars);

void sec_reset(TaylorModel &tmReset, const Interval intC, const int varStoreInd, const int varInputInd, const int numVars, const TaylorModelVec tmvAggregation, const std::vector<Interval> doAggregation);

void arc_reset(TaylorModel &tmReset, const Interval intC, const int varStoreInd, const int varInputInd, const int numVars);

void tan_reset(TaylorModel &tmReset, const Interval intC, const int varStoreInd, const int varInputInd, const int numVars);

void cos_reset(TaylorModel &tmReset, const Interval intC, const int varStoreInd, const int varInputInd, const int numVars);

void sin_reset(TaylorModel &tmReset, const Interval intC, const int varStoreInd, const int varInputInd, const int numVars);

void sqrt_reset(TaylorModel &tmReset, const Interval intC, const int varStoreInd, const int varInputInd, const int numVars);
