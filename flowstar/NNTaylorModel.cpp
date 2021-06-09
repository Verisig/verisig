/*---
  Verisig: A Verification Tool for Autonomous Systems with Neural Network Components.
  Authors: Radoslav Ivanov, Taylor Carpenter, James Weimer, Rajeev Alur, George Pappas, Insup Lee.
  Email: Radoslav Ivanov <rivanov@seas.upenn.edu> if you have questions or comments.

  The code is released as is under the GNU General Public License (GPL).
---*/

#include "NNTaylorModel.h"
#include "expression.h"
#include <chrono>

using namespace flowstar;
using namespace std;
ctpl::thread_pool* NNTaylorModelVec::nnpool = NULL;

double NNTaylorModel::mulTime = 0;
double NNTaylorModel::remTime = 0;
double NNTaylorModel::addAssignTime = 0;
double NNTaylorModel::cutoffTime = 0;
/*----------------------------------NNTaylorModel constructors----------------------------------*/

NNTaylorModel::NNTaylorModel()
{
}

NNTaylorModel::NNTaylorModel(const NNTaylorModel & tm)
{
        expansion = NNPolynomial(tm.expansion);
	remainder = tm.remainder;
}

NNTaylorModel::NNTaylorModel(const TaylorModel tm, const std::vector<std::string> & varNames)
{
        expansion = NNPolynomial(tm.expansion, varNames);
	remainder = tm.remainder;
}

NNTaylorModel::NNTaylorModel(const NNPolynomial & polyExp, const Interval & I)
{
        expansion = NNPolynomial(polyExp);
	remainder = I;
}

NNTaylorModel::NNTaylorModel(const string & strPolynomial, const Variables & vars)
{
	parsePolynomial.clear();

        string prefix(str_prefix_multivariate_polynomial);
        string suffix(str_suffix);

	parsePolynomial.strPolynomial = prefix + strPolynomial + suffix;
	parsePolynomial.variables = vars;

	parseMultivariatePolynomial();

	this->expansion = NNPolynomial(parsePolynomial.result, vars.varNames);
}

NNTaylorModel::NNTaylorModel(const Interval & I, const int numVars)
{
	NNPolynomial polyTemp(I, numVars);
	
	expansion = polyTemp;
	remainder.set(0.0, 0.0);
}

NNTaylorModel::NNTaylorModel(const iMatrix coefficients, const int rowIndex, const bool noTime)
{

        std::vector<Interval> new_coefficients;
	Interval intZero;

	if (noTime){
		new_coefficients.push_back(intZero); // add time
	}

	int cols = coefficients.cols();

	for(int i = 0; i < cols; i++){
	        new_coefficients.push_back(coefficients.getDataAt(rowIndex * cols + i));
	}

	NNPolynomial polyTemp(new_coefficients);
	expansion = polyTemp;
	remainder.set(0.0, 0.0);	

}

/*----------------------------------NNTaylorModelVec constructors----------------------------------*/

NNTaylorModelVec::NNTaylorModelVec()
{
}

NNTaylorModelVec::NNTaylorModelVec(const std::vector<Interval> & coefficients)
{
	int numVars = coefficients.size() + 1;
	for(int i=0; i<coefficients.size(); ++i)
	{
		NNTaylorModel tmTemp(coefficients[i], numVars);
		tmTemp.expansion.mul_assign(i+1, 1);
		tms.push_back(tmTemp);
	}
}

NNTaylorModelVec::NNTaylorModelVec(const std::vector<Interval> & constants, const int numVars)
{
	for(int i=0; i<constants.size(); ++i)
	{
		NNTaylorModel tmTemp(constants[i], numVars);
		tms.push_back(tmTemp);
	}
}

NNTaylorModelVec::NNTaylorModelVec(const iMatrix & coefficients, bool noTime)
{

        int cols = coefficients.cols();

	if(noTime) cols ++; //add an extra column for time
	
	int rows = coefficients.rows();
	
	for(int i=0; i<rows; ++i)
	{

	        NNTaylorModel tmTemp(coefficients, i, noTime);
		tms.push_back(tmTemp);
	}
}

/*----------------------------------NNTaylorModel methods----------------------------------*/

void NNTaylorModel::add(NNTaylorModel & result, const NNTaylorModel & tm) const
{
	result.expansion = expansion + tm.expansion;
	result.remainder = remainder + tm.remainder;
}

void NNTaylorModel::add_assign(const NNTaylorModel & tm)
{
        expansion.add_assign(tm.expansion);
  
	remainder += tm.remainder;
}

void NNTaylorModel::add_assign(const std::vector<std::shared_ptr<NNTaylorModel>> &tms)
{

	for(int i = 0; i < tms.size(); i++){
	        expansion.add_assign(tms[i]->expansion);
	        remainder += tms[i]->remainder;
	}
}

void NNTaylorModel::sub(NNTaylorModel & result, const NNTaylorModel & tm) const
{
	result.expansion = expansion - tm.expansion;
	result.remainder = remainder - tm.remainder;
}

void NNTaylorModel::intEvalNormal(Interval & result, const vector<Interval> & step_exp_table) const
{
	expansion.intEvalNormal(result, step_exp_table);
	result += remainder;
}

void NNTaylorModel::insert(NNTaylorModel & result, const NNTaylorModelVec & vars,
			   const vector<Interval> & varsPolyRange, const vector<Interval> & step_exp_table,
			   const bool ctrunc, const bool cutoff,
			   const Interval & cutoff_threshold, const int order) const
{
        
	if(vars.tms.size() == 0)
	{
	        result = *this;

		for(auto iter = result.expansion.monomials_map.begin(); iter != result.expansion.monomials_map.end(); )
		{
		        if( ((iter->second)->degree() - (iter->second)->getDegree(0)) > 0 )
			{
				iter = result.expansion.monomials_map.erase(iter);
			}
			else
			{
				++iter;
			}
		}
	}
	else
	{

	        NNTaylorModel::mulTime = 0;
		NNTaylorModel::remTime = 0;
		NNTaylorModel::addAssignTime = 0;
		NNTaylorModel::cutoffTime = 0;

		auto start = chrono::high_resolution_clock::now(); 		
	        shared_ptr<NNHornerForm> hf(new NNHornerForm());
		expansion.toHornerForm(hf);
		auto stop = chrono::high_resolution_clock::now();
		auto duration = chrono::duration_cast<chrono::nanoseconds>(stop - start);
		double toHFtime = duration.count();
		
		//hf->insertNN(result, vars, varsPolyRange, step_exp_table, ctrunc, cutoff, cutoff_threshold, order);

		NNPolynomial result_poly;
		Interval result_rem;

		start = chrono::high_resolution_clock::now();

		hf->insertNN(result_poly, result_rem, vars, varsPolyRange, step_exp_table,
			     ctrunc, cutoff, cutoff_threshold, order);
		stop = chrono::high_resolution_clock::now();
		duration = chrono::duration_cast<chrono::nanoseconds>(stop - start);	

		result.expansion = result_poly;
		result.remainder = result_rem;		
		
		result.remainder += remainder;
	}
}

void NNTaylorModel::polyRange(Interval & result, const std::vector<Interval> & step_exp_table) const
{
	expansion.intEvalNormal(result, step_exp_table);
}

void NNTaylorModel::constant(Interval & result) const
{
	expansion.constant(result);
}

void NNTaylorModel::addConstant(const Interval & constant, const int num_vars)
{

        expansion.addConstant(constant, num_vars);
  
}

void NNTaylorModel::mul(NNTaylorModel & result, const Interval & I) const
{
	expansion.mul(result.expansion, I);
	result.remainder.set(remainder * I);
}

void NNTaylorModel::mul(shared_ptr<NNTaylorModel> & result, const Interval & I) const
{
	expansion.mul(result->expansion, I);
	result->remainder.set(remainder * I);
}

void NNTaylorModel::mul_assign(const int varIndex, const int degree)
{
	expansion.mul_assign(varIndex, degree);
}

void NNTaylorModel::mul_assign(const Interval & I)
{
        expansion.mul_assign(I);
	remainder *= I;
}

void NNTaylorModel::mul_insert(NNTaylorModel & result, const NNTaylorModel & tm, const Interval & tmPolyRange,
			       const vector<Interval> & step_exp_table, const bool ctrunc, const bool cutoff, 
			       const Interval & cutoff_threshold, const int order) const
{
        Interval temp, temp_res;
	expansion.get_remainder_from_mul(result.remainder, remainder, tm.remainder,
					 tmPolyRange, temp_res, temp, step_exp_table);

	result.expansion = expansion * tm.expansion;

	if(ctrunc){
	        result.ctrunc(step_exp_table, order);
	}
	if(cutoff){
	        result.cutoff(step_exp_table, cutoff_threshold);
	}
}

void NNTaylorModel::mul_insert_assign(const NNTaylorModel & tm, const Interval & tmPolyRange,
				      const vector<Interval> & step_exp_table, const bool ctrunc, const bool cutoff, 
				      const Interval & cutoff_threshold, const int order)
{
        Interval temp, temp_res;
	expansion.get_remainder_from_mul(remainder, remainder, tm.remainder,
					 tmPolyRange, temp_res, temp, step_exp_table);

        expansion = expansion * tm.expansion;  

	if(ctrunc){
	        this->ctrunc(step_exp_table, order);
	}
	if(cutoff){
	        this->cutoff(step_exp_table, cutoff_threshold);
	}
}

void NNTaylorModel::rmConstant()
{
	expansion.rmConstant();
}

void NNTaylorModel::ctrunc(const std::vector<Interval> & step_exp_table, const int order)
{
	Interval I;
	expansion.ctrunc(I, step_exp_table, order);

	remainder += I;
}

void NNTaylorModel::cutoff(const vector<Interval> & step_exp_table, const Interval & cutoff_threshold)
{
	Interval I;
	expansion.cutoff(I, step_exp_table, cutoff_threshold);

	remainder += I;
}

void NNTaylorModel::clear()
{
	expansion.clear();
	remainder.set(0.0, 0.0);
}

/*----------------------------------NNTaylorModelVec methods----------------------------------*/

void NNTaylorModelVec::add(NNTaylorModelVec & result, const NNTaylorModelVec & tmv) const
{
	result.clear();

	if(tms.size() != tmv.tms.size())
	{
		printf("Dimensions do not coincide.\n");
		return;
	}

	for(int i=0; i<tms.size(); ++i)
	{
		NNTaylorModel tmTemp;
		tms[i].add(tmTemp, tmv.tms[i]);
		result.tms.push_back(tmTemp);
	}
}

void NNTaylorModelVec::add_assign(const NNTaylorModelVec & tmv)
{

	if(tms.size() != tmv.tms.size())
	{
		printf("Dimensions do not coincide.\n");
		return;
	}

	for(int i=0; i<tms.size(); ++i)
	{
		tms[i].add_assign(tmv.tms[i]);
	}	
}

void NNTaylorModelVec::sub(NNTaylorModelVec & result, const NNTaylorModelVec & tmv) const
{
	result.clear();

	if(tms.size() != tmv.tms.size())
	{
		printf("Dimensions do not coincide.\n");
		return;
	}

	for(int i=0; i<tms.size(); ++i)
	{
		NNTaylorModel tmTemp;
		tms[i].sub(tmTemp, tmv.tms[i]);
		result.tms.push_back(tmTemp);
	}
}

void NNTaylorModelVec::insert(NNTaylorModelVec & result, const NNTaylorModelVec & vars,
			      const std::vector<Interval> & varsPolyRange,
			      const std::vector<Interval> & step_exp_table, const bool ctrunc, const bool cutoff,
			      const Interval & cutoff_threshold, const int order) const
{
	result.clear();

	std::vector<std::future<NNTaylorModel>> futures;

	for(int i=0; i<tms.size(); ++i)
	{
		auto tm = tms[i];
		futures.push_back(nnpool->push(
					       [tm, vars, varsPolyRange, step_exp_table,
						ctrunc, cutoff, cutoff_threshold, order](int id) {
					       NNTaylorModel tmTemp;
					       tm.insert(tmTemp, vars, varsPolyRange, step_exp_table,
							 ctrunc, cutoff, cutoff_threshold, order);

					       return tmTemp;
					     }
					     ));
	}

	for(int i=0; i<tms.size(); ++i) {
		result.tms.push_back(futures[i].get());
	}
}

void NNTaylorModelVec::linearTrans(NNTaylorModelVec & result, const iMatrix & A) const
{
	result.clear();
	
	if(tms.size() != A.cols())
	{
		printf("Dimensions do not coincide.\n");
		return;
	}

	Interval I;

	int rows = A.rows();
	int cols = A.cols();
	for(int i=0; i<rows; ++i)
	{		

	  	NNTaylorModel tm1;

		for(int j=0; j<cols; ++j)
		{

			A.getDataAt(I, j + i * cols);
			if(I.subsetZero()) continue;

		        shared_ptr<NNTaylorModel> tm2(new NNTaylorModel());			

			tms[j].mul(tm2, I);
			tm1.add_assign(*tm2);
		}
		
		result.tms.push_back(tm1);
	}
}

void NNTaylorModelVec::linearTrans_par(NNTaylorModelVec & result, const iMatrix & A) const
{
	result.clear();
	
	if(tms.size() != A.cols())
	{
		printf("Dimensions do not coincide.\n");
		return;
	}

	int rows = A.rows();
	int cols = A.cols();

	std::vector<std::future<NNTaylorModel>> futures;
	
	for(int i=0; i<rows; ++i)
	{
	  
	        futures.push_back(nnpool->push(
					       [cols, i, A, this](int id){
						       vector<shared_ptr<NNTaylorModel>> all_tms;
						       Interval intZero;
						       for(int j=0; j<cols; ++j)
						       {
							       shared_ptr<NNTaylorModel> tm2(new NNTaylorModel());

							       Interval I = A.getDataAt(j + i * cols);
							       if(I.subseteq(intZero)) continue;

							       tms[j].mul(tm2, I);
							       all_tms.push_back(tm2);
						       }

						       NNTaylorModel tm1;
						       tm1.add_assign(all_tms);
						       return tm1;
					       }));
	}

	for(int i=0; i<tms.size(); ++i) {
		result.tms.push_back(futures[i].get());
	}	
}

void NNTaylorModelVec::addConstant(const iMatrix & A)
{
	
	if(A.cols() != 1)
	{
		printf("Can only add constants through a vector.\n");
		return;
	}

	int rows = A.rows();
	for(int i=0; i<rows; ++i)
	{
		Interval I = A.getDataAt(i);

		if(I.subsetZero()) continue;

		tms[i].addConstant(I, tms.size());
	}
}

void NNTaylorModelVec::polyRange(std::vector<Interval> & result, const std::vector<Interval> & step_exp_table) const
{
	result.clear();

	for(int i=0; i<tms.size(); ++i)
	{
		Interval intTemp;
		tms[i].polyRange(intTemp, step_exp_table);
		result.push_back(intTemp);
	}
}

void NNTaylorModelVec::intEvalNormal(std::vector<Interval> & result, const std::vector<Interval> & step_exp_table) const
{
	result.clear();
	for(int i=0; i<tms.size(); ++i)
	{
		Interval I;
		tms[i].intEvalNormal(I, step_exp_table);
		result.push_back(I);
	}
}

void NNTaylorModelVec::constant(std::vector<Interval> & result) const
{
	result.clear();

	for(int i=0; i<tms.size(); ++i)
	{
		Interval I;
		tms[i].constant(I);
		result.push_back(I);
	}
}

void NNTaylorModelVec::constant(iMatrix & result) const
{
        result = iMatrix(tms.size(), 1);

	for(int i=0; i<tms.size(); ++i)
	{
		Interval I;
		tms[i].constant(I);
		result.setDataAt(i, 0, I);
	}
}

void NNTaylorModelVec::rmConstant()
{
	for(int i=0; i<tms.size(); ++i)
		tms[i].rmConstant();
}

void NNTaylorModelVec::scale_assign(const std::vector<Interval> & S)
{
	for(int i=0; i<tms.size(); ++i)
	{
		tms[i].mul_assign(S[i]);
	}
}

void NNTaylorModelVec::cutoff(const std::vector<Interval> & step_exp_table, const Interval & cutoff_threshold)
{
	for(int i=0; i<tms.size(); ++i)
	{
		tms[i].cutoff(step_exp_table, cutoff_threshold);
	}
}

void NNTaylorModelVec::clear()
{
	tms.clear();
}
