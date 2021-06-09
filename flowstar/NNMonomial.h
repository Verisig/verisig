/*---
  Verisig: A Verification Tool for Autonomous Systems with Neural Network Components.
  Authors: Radoslav Ivanov, Taylor Carpenter, James Weimer, Rajeev Alur, George Pappas, Insup Lee.
  Email: Radoslav Ivanov <rivanov@seas.upenn.edu> if you have questions or comments.

  The code is released as is under the GNU General Public License (GPL).
---*/

#ifndef NN_MONOMIAL_H_
#define NN_MONOMIAL_H_

#include "Interval.h"
#include <memory>

namespace flowstar
{

class Monomial;

class NNMonomial
{
protected:
        std::shared_ptr<Interval> coefficient;		// the coefficient of the monomial
	std::map<int, int> degrees_map;
	int d;			        		// the degree of the monomial, it is the sum of the values in degrees.
	int num_vars;
public:
	NNMonomial(); // empty monomial.
	NNMonomial(const Interval & I, const std::vector<int> & degs);	
	NNMonomial(const Monomial &monomial); // convert between the two types	
	NNMonomial(const Interval & I, const int numVars);  // a constant
	NNMonomial(const NNMonomial & monomial);
	NNMonomial(const std::shared_ptr<NNMonomial> & monomial);
	NNMonomial(const std::shared_ptr<NNMonomial> & monomial, const bool shallow); // option for a shallow copy
	//NNMonomial(const Interval & I, const std::map<int, int> & degs, int numVars);

	int degree() const; // degree of the monomial
	Interval getCoefficient() const;
	void setDegree(const int index, const int degree);
	void addDegree(const int index, const int degree);
	int getDegree(const int index);
	int getNumVars();

	bool isLinear(int & index) const;

	void intEval(Interval & result, const std::vector<Interval> & domain) const;	// interval evaluation of the monomial
	void intEvalNormal(Interval & result, const std::vector<Interval> & step_exp_table) const;

	void add_assign(const std::shared_ptr<NNMonomial> &monomial);	

	void mul(std::shared_ptr<NNMonomial> &result, const std::shared_ptr<NNMonomial> &monomial);

	NNMonomial & operator += (const NNMonomial & monomial);			// we assume the two monomials can be added up
	NNMonomial & operator *= (const NNMonomial & monomial);
	const NNMonomial operator + (const NNMonomial & monomial) const;
	const NNMonomial operator * (const NNMonomial & monomial) const;
	
	bool operator < (const NNMonomial & b) const;  // Define a partial order over the monomials
	bool operator == (const NNMonomial & b) const;

	int cutoff(std::shared_ptr<NNMonomial> & monoRem, const Interval & cutoff_threshold);

	void toString(std::string & result, const std::vector<std::string> & varNames) const;
	void toStringNoCoef(std::string & result, const std::vector<std::string> & varNames) const;	

	friend class Monomial;
	friend class NNPolynomial;	
	friend class NNTaylorModel;
	friend class NNTaylorModelVec;
};
 
}

bool compare(const std::shared_ptr<flowstar::NNMonomial>& lhs, const std::shared_ptr<flowstar::NNMonomial>& rhs);

#endif /* NN_MONOMIAL_H_ */
