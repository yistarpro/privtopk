#ifndef EIF_TESTCODE_H
#define EIF_TESTCODE_H

#include "openfhe.h"
#include <iostream>

using namespace lbcrypto;
using namespace std;

namespace ckkseif {

    string statTime(const vector<double> times, const usint iteration);

	void SortTest(const uint32_t scaleModSize, const usint size, const int32_t arraybound, const usint levelBudgetElmt, const usint logbatchsize, usint iteration);

	void TopkTest_MT_opt(const uint32_t scaleModSize, const usint size, const int32_t arraybound, const usint numk, const usint levelBudgetElmt, usint iteration, const int32_t mtenvironment = 0, const int32_t initial_subsize_bound = 64, const int32_t init_numct_bound = 128, const bool optimized = false, const int32_t custom_mergecrit = 1);

	void nexusTest(const uint32_t scaleModSize, const usint size, const int32_t arraybound,const usint levelBudgetElmt, usint iteration);


}
#endif
