/*---
  Verisig: A Verification Tool for Autonomous Systems with Neural Network Components.
  Authors: Radoslav Ivanov, Taylor Carpenter, James Weimer, Rajeev Alur, George Pappas, Insup Lee.
  Email: Radoslav Ivanov <rivanov@seas.upenn.edu> if you have questions or comments.

  The code is released as is under the GNU General Public License (GPL).
---*/

#ifndef NNTAYLORMODEL_H_
#define NNTAYLORMODEL_H_

#include "Polynomial.h"
#include "NNPolynomial.h"

#include "ctpl/ctpl_stl.h"

namespace flowstar
{  
  
class NNTaylorModelVec;
  
class NNTaylorModel
{
public:
	NNPolynomial expansion;	// Taylor expansion
	Interval remainder;		// remainder interval

	static double addAssignTime;
	static double mulTime;
	static double remTime;
	static double cutoffTime;
public:
	NNTaylorModel();			// empty Taylor model.
        NNTaylorModel(const NNTaylorModel & tm);
	NNTaylorModel(const TaylorModel tm, const std::vector<std::string> & varNames);
	NNTaylorModel(const NNPolynomial & polyExp, const Interval & I);	// Taylor model (P,I)
	NNTaylorModel(const Interval & I, const int numVars);				// constant
	NNTaylorModel(const iMatrix coefficients, const int rowIndex, const bool noTime); // added by Rado
	NNTaylorModel(const std::string & strPolynomial, const Variables & vars);	
  
	void constant(Interval & result) const;      // Return the constant part of the expansion.
	void addConstant(const Interval & constant, const int num_vars);		// add a constant
	void rmConstant();	

	void insert(NNTaylorModel & result, const NNTaylorModelVec & vars, const std::vector<Interval> & varsPolyRange,
		    const std::vector<Interval> & step_exp_table, const bool ctrunc, const bool cutoff,
		    const Interval & cutoff_threshold, const int order) const;
	
	void polyRange(Interval & result, const std::vector<Interval> & step_exp_table) const;
	void intEvalNormal(Interval & result, const std::vector<Interval> & step_exp_table) const;	

	void add(NNTaylorModel & result, const NNTaylorModel & tm) const;			// addition	
	void add_assign(const NNTaylorModel & tm);

	void add_assign(const std::vector<std::shared_ptr<NNTaylorModel>> &tms); // a more efficient way to add a vector of TMs

	void sub(NNTaylorModel & result, const NNTaylorModel & tm) const;			// subtraction	
	
	void mul(NNTaylorModel & result, const Interval & I) const;
	void mul(std::shared_ptr<NNTaylorModel> & result, const Interval & I) const;
	void mul_assign(const Interval & I);
	void mul_assign(const int varIndex, const int degree);

	void mul_insert(NNTaylorModel & result, const NNTaylorModel & tm, const Interval & tmPolyRange,
			const std::vector<Interval> & step_exp_table, const bool ctrunc, const bool cutoff,
			const Interval & cutoff_threshold, const int order) const;
	void mul_insert_assign(const NNTaylorModel & tm, const Interval & tmPolyRange,
			       const std::vector<Interval> & step_exp_table, const bool ctrunc, const bool cutoff,
			       const Interval & cutoff_threshold, const int order);

	void cutoff(const std::vector<Interval> & step_exp_table, const Interval & cutoff_threshold);
	
	void ctrunc(const std::vector<Interval> & step_exp_table, const int order);		// conservative truncation	

	void clear();	


	friend class NNTaylorModelVec;
	friend class Monomial;
	friend class Polynomial;
};

class NNTaylorModelVec			// Taylor models: R^n -> R^m
{
public:
	static ctpl::thread_pool* nnpool;
public:
	std::vector<NNTaylorModel> tms;
public:
	NNTaylorModelVec();
	NNTaylorModelVec(const std::vector<Interval> & coefficients);
	NNTaylorModelVec(const std::vector<Interval> & constants, const int numVars);
	NNTaylorModelVec(const iMatrix & coefficients, bool noTime); //added by Rado

	void constant(std::vector<Interval> & result) const;
	void constant(iMatrix & result) const;

	void rmConstant();

	void insert(NNTaylorModelVec & result, const NNTaylorModelVec & vars,
		    const std::vector<Interval> & varsPolyRange,
		    const std::vector<Interval> & step_exp_table, const bool ctrunc, const bool cutoff, 
		    const Interval & cutoff_threshold, const int order) const;
	void polyRange(std::vector<Interval> & result, const std::vector<Interval> & step_exp_table) const;
	void intEvalNormal(std::vector<Interval> & result, const std::vector<Interval> & step_exp_table) const;	

	void add(NNTaylorModelVec & result, const NNTaylorModelVec & tmv) const;		// addition	
	void add_assign(const NNTaylorModelVec & tmv);

	void sub(NNTaylorModelVec & result, const NNTaylorModelVec & tmv) const;		// subtraction

	void scale_assign(const std::vector<Interval> & S);

	void addConstant(const iMatrix & A);		// add a vector of constants
	void linearTrans(NNTaylorModelVec & result, const iMatrix & A) const;		// linear transformation
	void linearTrans_par(NNTaylorModelVec & result, const iMatrix & A) const;		// linear transformation (parallelized)

	void cutoff(const std::vector<Interval> & step_exp_table, const Interval & cutoff_threshold);	

	void clear();
};

}

#endif /* NNTAYLORMODEL_H_ */
