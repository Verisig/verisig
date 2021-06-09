/*---
  Verisig: A Verification Tool for Autonomous Systems with Neural Network Components.
  Authors: Radoslav Ivanov, Taylor Carpenter, James Weimer, Rajeev Alur, George Pappas, Insup Lee.
  Email: Radoslav Ivanov <rivanov@seas.upenn.edu> if you have questions or comments.

  The code is released as is under the GNU General Public License (GPL).
---*/

#include "Monomial.h"
#include "NNMonomial.h"

using namespace flowstar;
using namespace std;

/*----------------------------------NNMonomial constructors----------------------------------*/

NNMonomial::NNMonomial()
{
	d = 0;
	num_vars = 0;
}

NNMonomial::NNMonomial(const Monomial &monomial)
{
        coefficient = make_shared<Interval>(monomial.coefficient);
	d = monomial.d;

	for(int i=0; i<monomial.degrees.size(); ++i)
	{
		if(monomial.degrees[i] > 0){
		        degrees_map[i] = monomial.degrees[i];
		}
	}
	num_vars = monomial.degrees.size();
}

NNMonomial::NNMonomial(const Interval & I, const std::vector<int> & degs):d(0)
{
	for(int i=0; i<degs.size(); ++i)
	{
		d += degs[i];

		if(degs[i] > 0){
		        degrees_map[i] = degs[i];
		}
	}
	num_vars = degs.size();
	coefficient = make_shared<Interval>(I);
}

NNMonomial::NNMonomial(const Interval & I, const int numVars):d(0)
{

	coefficient = make_shared<Interval>(I);
	num_vars = numVars;
}

NNMonomial::NNMonomial(const NNMonomial & monomial): d(monomial.d)
{
	degrees_map = monomial.degrees_map;
	num_vars = monomial.num_vars;
	coefficient = make_shared<Interval>(*monomial.coefficient);
}

NNMonomial::NNMonomial(const shared_ptr<NNMonomial> & monomial)
{
        coefficient = make_shared<Interval>(*monomial->coefficient);
	d = monomial->d;
	degrees_map = monomial->degrees_map;
	num_vars = monomial->num_vars;	
}

NNMonomial::NNMonomial(const shared_ptr<NNMonomial> & monomial, const bool shallow)
{
        if(shallow)
	        coefficient = monomial->coefficient;
	else
	        coefficient = make_shared<Interval>(*monomial->coefficient);
  
	d = monomial->d;
	degrees_map = monomial->degrees_map;
	num_vars = monomial->num_vars;	
}

/*----------------------------------NNMonomial methods----------------------------------*/

NNMonomial & NNMonomial::operator += (const NNMonomial & monomial)
{
        (*coefficient) += (*monomial.coefficient);
	return *this;
}

NNMonomial & NNMonomial::operator *= (const NNMonomial & monomial)
{
        (*coefficient) *= (*monomial.coefficient);

	for(auto i=monomial.degrees_map.begin(); i != monomial.degrees_map.end(); ++i){
	        if(degrees_map.find(i->first) == degrees_map.end()){
		        degrees_map[i->first] = 0;
		}

		degrees_map[i->first] += i->second;
	}	

	d += monomial.d;
	return *this;
}

const NNMonomial NNMonomial::operator + (const NNMonomial & monomial) const
{
	NNMonomial result = *this;
	result += monomial;
	return result;
}

const NNMonomial NNMonomial::operator * (const NNMonomial & monomial) const
{
	NNMonomial result = *this;
	result *= monomial;
	return result;
}

void NNMonomial::setDegree(const int index, const int degree)
{

        if(degrees_map.find(index) == degrees_map.end())
	        degrees_map[index] = 0;
	
        d += (degree - degrees_map[index]);
	
	degrees_map[index] = degree;

	if(degrees_map[index] <= 0)
	        degrees_map.erase(index);
}

void NNMonomial::addDegree(const int index, const int degree)
{
        d += degree;
        
	if(degrees_map.find(index) == degrees_map.end()){
	        degrees_map[index] = 0;
	}
	
	degrees_map[index] += degree;

	if(degrees_map[index] <= 0)
	        degrees_map.erase(index);
}

int NNMonomial::getDegree(const int index)
{

        if(degrees_map.find(index) == degrees_map.end())
	        return 0;

	return degrees_map[index];
}

int NNMonomial::getNumVars()
{
        return num_vars;
}

void NNMonomial::add_assign(const std::shared_ptr<NNMonomial> &monomial)
{
        (*coefficient) += (*monomial->coefficient);
}

void NNMonomial::mul(std::shared_ptr<NNMonomial> &result, const std::shared_ptr<NNMonomial> &monomial)
{
        //result->coefficient = make_shared<Interval>((*coefficient) * (*monomial->coefficient));
  
        //two operations to avoid a copy
        result->coefficient->set((*coefficient));
	(*result->coefficient) *= (*monomial->coefficient);

	result->degrees_map = degrees_map;
        for(auto i=monomial->degrees_map.begin(); i != monomial->degrees_map.end(); ++i){
	        if(result->degrees_map.find(i->first) == result->degrees_map.end()){
		        result->degrees_map[i->first] = 0;
		}

		result->degrees_map[i->first] += monomial->degrees_map[i->first];
	}

	result->num_vars = num_vars;	
	result->d = d + monomial->d;
}

bool NNMonomial::operator < (const NNMonomial & b) const
{
	if(d < b.d)
	{
		return true;
	}
	else if(d > b.d)
	{
		return false;
	}
	else	// a.d == b.d
	{
	        auto A=degrees_map.begin();
		auto B=b.degrees_map.begin();
	        for(; A != degrees_map.end(); ++A, ++B)	
		{
		        if(A->first < B->first || (A->first == B->first && A->second < B->second))
				return true;
			if(B->first < A->first || (A->first == B->first && B->second < A->second))
				return false;
		}
	}

	return false;	// a == b
}

bool NNMonomial::operator == (const NNMonomial & b) const
{
	if (d == b.d)
	{
	        for(auto i=degrees_map.begin(); i != degrees_map.end(); ++i)		  
		{
		        if(b.degrees_map.find(i->first) == b.degrees_map.end() || i->second != b.degrees_map.find(i->first)->second)
				return false;
		}

		return true;	// The two monomials are identical without considering the coefficients.
	}
	else
		return false;
}

int NNMonomial::degree() const
{
	return d;
}

Interval NNMonomial::getCoefficient() const
{

        return (*coefficient);
}

bool NNMonomial::isLinear(int & index) const
{
	if(d == 1)
	{
	        for(auto i=degrees_map.begin(); i != degrees_map.end(); ++i) 		  
		{
			if(i->second == 1)
			{
				index = i->first;
				return true;
			}
		}
	}

	return false;
}

void NNMonomial::intEval(Interval & result, const std::vector<Interval> & domain) const
{
        result = (*coefficient);

	for(auto i=degrees_map.begin(); i != degrees_map.end(); ++i)	  
	{
		Interval tmpI(1,1);
		for(int j=0; j<i->second; ++j)
		{
			tmpI *= i->first;
		}

		result *= tmpI;
	}
}

void NNMonomial::intEvalNormal(Interval & result, const std::vector<Interval> & step_exp_table) const
{

        result.set(0.0, 0.0);

	if(num_vars == 0 || coefficient->subsetZero())
		return;

	result.set(*coefficient);

	bool bSet = false;

	for(auto i=degrees_map.begin(); i != degrees_map.end(); ++i)  
	{
		if(i->second%2 == 0)	// degree is an even number
		{
			if(!bSet)
			{
				bSet = true;
			}
		}
		else			// degree is an odd number
		{
		        result.mul_m11();
		        return;
		}
	}

	if(bSet)
	        result.mul_01();
}

int NNMonomial::cutoff(std::shared_ptr<NNMonomial> & monoRem, const Interval & cutoff_threshold)
{
	Interval M;
	coefficient->midpoint(M);

	if(M.subseteq(cutoff_threshold))
	{
		return 2;
	}
	else if(coefficient->width() >= MAX_WIDTH)
	{
	        monoRem = shared_ptr<NNMonomial>(new NNMonomial(*this));
		monoRem->coefficient->remove_midpoint();
		coefficient->set(M);
		return 1;
	}
	else
	{
		return 0;
	}
}

void NNMonomial::toString(string & result, const vector<string> & varNames) const
{
        string strMono;

	strMono += '(';

	std::string strInt;
	coefficient->toString(strInt);
	strMono += strInt;

	for(auto i=degrees_map.begin(); i != degrees_map.end(); ++i)
	{

	        if(i->second == 1)
		{
		        strMono += ' ';
			strMono += '*';
			strMono += ' ';
			strMono += varNames[i->first];
		}
		else
		{
		        strMono += ' ';
			strMono += '*';
			strMono += ' ';
			strMono += varNames[i->first];
			strMono += '^';

			char strNum[NUM_LENGTH];
			sprintf(strNum, "%d", i->second);
			string num(strNum);
			strMono += num;
		}
		
	}

	strMono += ')';

	result = strMono;
}

void NNMonomial::toStringNoCoef(string & result, const vector<string> & varNames) const
{
        string strMono;

	std::string strInt;

	bool isFirst = true;

	for(auto i=degrees_map.begin(); i != degrees_map.end(); ++i)
	{
	        if(i->second == 1)
		{
		        if(isFirst){
			        isFirst = false;
			}
			else{
			        strMono += '*';
			}
	
			strMono += varNames[i->first];
		}
		else
		{
		        if(isFirst){
			        isFirst = false;
			}
			else{
			        strMono += '*';
			}
			strMono += varNames[i->first];
			strMono += '^';

			char strNum[NUM_LENGTH];
			sprintf(strNum, "%d", i->second);
			string num(strNum);
			strMono += num;
		}
	}

	result = strMono;
}

bool compare(const std::shared_ptr<NNMonomial>& lhs, const std::shared_ptr<NNMonomial>& rhs) {
        return *lhs < *rhs;
}
