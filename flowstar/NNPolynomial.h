/*---
  Verisig: A Verification Tool for Autonomous Systems with Neural Network Components.
  Authors: Radoslav Ivanov, Taylor Carpenter, James Weimer, Rajeev Alur, George Pappas, Insup Lee.
  Email: Radoslav Ivanov <rivanov@seas.upenn.edu> if you have questions or comments.

  The code is released as is under the GNU General Public License (GPL).
---*/

#ifndef NNPOLYNOMIAL_H_
#define NNPOLYNOMIAL_H_

#include "Monomial.h"
#include "Polynomial.h"
#include "Matrix.h"
#include <memory>

namespace flowstar
{

class NNMonomial;
class NNPolynomial;
class NNTaylorModel;
class NNTaylorModelVec;

class NNHornerForm    // c + (...)*x1 + (...)*x2 + ... + (...)*xn
{
protected:
        std::shared_ptr<Interval> constant;  // constant part
	std::vector<std::shared_ptr<NNHornerForm>> hornerForms; // other parts
public:
	NNHornerForm();

	void intEval(Interval & result, const std::vector<Interval> & domain) const;		// interval evaluation of the Horner form

	void insertNN(NNTaylorModel & result, const NNTaylorModelVec & vars, const std::vector<Interval> & varsPolyRange,
		      const std::vector<Interval> & step_exp_table, const bool ctrunc, const bool cutoff,
		      const Interval & cutoff_threshold, const int order) const;

	void insertNN(NNPolynomial &poly, Interval &remainder,
		      const NNTaylorModelVec & vars, const std::vector<Interval> & varsPolyRange,
		      const std::vector<Interval> & step_exp_table, const bool ctrunc, const bool cutoff,
		      const Interval & cutoff_threshold, const int order) const;

	void clear();

	friend class NNTaylorModel;
	friend class NNTaylorModelVec;
	friend class NNPolynomial;
};

class NNPolynomial				// polynomials in monomial form
{
public:
  //        std::list<std::shared_ptr<NNMonomial>> monomials;
	std::map<std::string, std::shared_ptr<NNMonomial>> monomials_map;
public:
	NNPolynomial();														// empty polynomial
	NNPolynomial(const Polynomial &poly, const std::vector<std::string> & varNames); // convert between the two polynomial classes
	NNPolynomial(const NNPolynomial &polynomial);
	NNPolynomial(const Interval & constant, const int numVars);  // constant polynomial where dim is the number of the variables
	NNPolynomial(const std::vector<Interval> & coefficients);   // linear polynomial with the given coefficients, the input matrix is a row vector
	NNPolynomial(const std::list<std::shared_ptr<NNMonomial>> & monos);
	NNPolynomial(const std::shared_ptr<NNMonomial> & monomial);   // polynomial with one monomial	

	void intEval(Interval & result, const std::vector<Interval> & domain) const;	// interval evaluation of the polynomial	
	void intEvalNormal(Interval & result, const std::vector<Interval> & step_exp_table) const;	// fast evaluation over normalized domain
	void intEvalNormalTemp(Interval & result, Interval & temp, const std::vector<Interval> & step_exp_table) const;	// fast evaluation over normalized domain		

	void inv(NNPolynomial & result) const;   // additive inverse
	void inv_assign();

	void add_assign(const NNPolynomial &polynomial);   // add a polynomial
	void add_assign(const std::shared_ptr<NNMonomial> & monomial);   // add a monomial
	void add_assign(const std::shared_ptr<NNMonomial> & monomial, const std::vector<std::string> & varNames);   // add a monomial; this should only be used when converting between NNTM and TM
	
	void add_assign(const std::list<std::shared_ptr<NNMonomial>> & monomials);   // add a list of monomial
	//void add_to_end(const std::shared_ptr<NNMonomial> & monomial);

	void get_remainder_from_mul(Interval &result_rem, const Interval &remainder,
				    const Interval &remainder_other, const Interval & tmPolyRange,
				    Interval &temp_res, Interval &temp, const std::vector<Interval> & step_exp_table) const;

	void mul(NNPolynomial & result, const Interval & I) const;
	//void mul_no_add(NNPolynomial &tempPoly, const NNPolynomial & polynomial) const;	
	void mul(NNPolynomial result, const int varIndex, const int degree) const;
	void mul_assign(const Interval & I);   // multiplied by an interval	
	void mul_assign(const int varIndex, const int degree);   // multiplied by a term x^d
	void mul_assign(const NNPolynomial & p, std::shared_ptr<NNMonomial> &mono_temp);   // multiplied by a polynomial

	void ctrunc(Interval & remainder, const std::vector<Interval> & step_exp_table, const int order);	// conservative truncation	
	void cutoff(Interval & intRem, const std::vector<Interval> & step_exp_table, const Interval & cutoff_threshold);

	void toHornerForm(std::shared_ptr<NNHornerForm> & result) const;          // transform the polynomial into a Horner form	

	void constant(Interval & result) const;   // constant part of the polynomial
	void addConstant(const Interval & constant, const int num_vars);		// add a vector of constants
	void rmConstant();				// remove the constant part

	void toString(std::string & result, const std::vector<std::string> & varNames) const;
	//	void toStringMap(std::string & result, const std::vector<std::string> & varNames) const;
	
	//void reorder();   // sort the monomials.
	void clear();	

	NNPolynomial & operator += (const NNPolynomial & polynomial);
	NNPolynomial & operator -= (const NNPolynomial & polynomial);
	NNPolynomial & operator *= (const NNPolynomial & polynomial);
	NNPolynomial & operator *= (const Interval & I);	
	
	NNPolynomial operator + (const NNPolynomial & polynomial) const;
	NNPolynomial operator - (const NNPolynomial & polynomial) const;
	NNPolynomial operator * (const NNPolynomial & polynomial) const;
	NNPolynomial operator * (const Interval & I) const;	
};
  
}


#endif /* NNPOLYNOMIAL_H_ */
