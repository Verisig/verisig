/*---
  Verisig: A Verification Tool for Autonomous Systems with Neural Network Components.
  Authors: Radoslav Ivanov, Taylor Carpenter, James Weimer, Rajeev Alur, George Pappas, Insup Lee.
  Email: Radoslav Ivanov <rivanov@seas.upenn.edu> if you have questions or comments.

  The code is released as is under the GNU General Public License (GPL).
---*/

#include "DNNResets.h"
#include "NNMonomial.h"
#include "NNPolynomial.h"
#include "NNTaylorModel.h"
#include <chrono>

using namespace flowstar;
using namespace std;

extern std::vector<std::string> dnn::curAugmentedVarNames;

/*----------------------------------NNPolynomial constructors----------------------------------*/

NNPolynomial::NNPolynomial()
{
}

NNPolynomial::NNPolynomial(const NNPolynomial &polynomial)
{
	
	for(auto i=polynomial.monomials_map.begin(); i != polynomial.monomials_map.end(); ++i){
	        monomials_map[i->first] = make_shared<NNMonomial>(i->second);
	}
}

NNPolynomial::NNPolynomial(const Polynomial &poly, const std::vector<std::string> & varNames)
{

        list<Monomial>::const_iterator iter;
	string mono_string;
	
	for(iter = poly.monomials.begin(); iter != poly.monomials.end(); ++iter){

		shared_ptr<NNMonomial> p(new NNMonomial((*iter)));
		p->toStringNoCoef(mono_string, varNames);
		monomials_map[mono_string] = p;
	}

}

NNPolynomial::NNPolynomial(const Interval & constant, const int numVars)
{

	if(!constant.subsetZero())
	{
		string mono_string;
		shared_ptr<NNMonomial> p(new NNMonomial(constant, numVars));
		p->toStringNoCoef(mono_string, dnn::curAugmentedVarNames);
		monomials_map[mono_string] = p;
	}
}

NNPolynomial::NNPolynomial(const std::vector<Interval> & coefficients)
{
	int numVars = coefficients.size();
	string mono_string;
	
	for(int i=0; i<numVars; ++i)
	{
		if(coefficients[i].subsetZero())		// the coefficient is zero
			continue;

		shared_ptr<NNMonomial> p(new NNMonomial(coefficients[i], numVars));
		p->setDegree(i, 1);
		p->toStringNoCoef(mono_string, dnn::curAugmentedVarNames);
		monomials_map[mono_string] = p;
	}
}

NNPolynomial::NNPolynomial(const list<shared_ptr<NNMonomial>> & monos)
{
        string mono_string;
	for(auto iter = monos.begin(); iter != monos.end(); ++iter){
		(*iter)->toStringNoCoef(mono_string, dnn::curAugmentedVarNames);
		monomials_map[mono_string] = make_shared<NNMonomial>(*iter);
	}
	
}

NNPolynomial::NNPolynomial(const std::shared_ptr<NNMonomial> & monomial)
{
	string mono_string;
	monomial->toStringNoCoef(mono_string, dnn::curAugmentedVarNames);
	monomials_map[mono_string] = make_shared<NNMonomial>(monomial);	
}

/*----------------------------------NNHornerForm constructors----------------------------------*/

NNHornerForm::NNHornerForm()
{
}

/*----------------------------------NNPolynomial methods----------------------------------*/

NNPolynomial & NNPolynomial::operator += (const NNPolynomial & polynomial)
{

	Interval intTemp;

	for(auto iterM = polynomial.monomials_map.begin(); iterM != polynomial.monomials_map.end(); iterM++){
	  
	        string mono_string;
		(iterM->second)->toStringNoCoef(mono_string, dnn::curAugmentedVarNames);

		if(monomials_map.find(mono_string) == monomials_map.end()){
		        monomials_map[mono_string] = make_shared<NNMonomial>(iterM->second);
		}
		else{
		        (*monomials_map[mono_string]) += (*iterM->second);
		}
	}

	return *this;
}

NNPolynomial NNPolynomial::operator + (const NNPolynomial & polynomial) const
{
	NNPolynomial result = *this;
	result += polynomial;
	return result;
}

NNPolynomial & NNPolynomial::operator -= (const NNPolynomial & polynomial)
{
	NNPolynomial polyTemp = polynomial;
	polyTemp.inv_assign();
	*this += polyTemp;

	return *this;
}

NNPolynomial NNPolynomial::operator - (const NNPolynomial & polynomial) const
{
	NNPolynomial result = *this;
	result -= polynomial;
	return result;
}

NNPolynomial & NNPolynomial::operator *= (const NNPolynomial & polynomial)
{
	if((monomials_map.size() == 0) || (polynomial.monomials_map.size() == 0))
	{
		this->clear();
		return *this;
	}

	map<string, shared_ptr<NNMonomial>> monomials_map_new;

	shared_ptr<NNMonomial> temp_mono(new NNMonomial());
	temp_mono->coefficient = make_shared<Interval>();
	string mono_string;
	
	for(auto iterA = monomials_map.begin(); iterA != monomials_map.end(); ++iterA){
	        for(auto iterB = polynomial.monomials_map.begin(); iterB != polynomial.monomials_map.end(); ++iterB){

		        (iterA->second)->mul(temp_mono, (iterB->second));

			temp_mono->toStringNoCoef(mono_string, dnn::curAugmentedVarNames);

			if(monomials_map_new.find(mono_string) == monomials_map_new.end()){
			        monomials_map_new[mono_string] = make_shared<NNMonomial>(temp_mono);
			}
			else{
			        (*monomials_map_new[mono_string]) += (*temp_mono);
			}

		}
	}

	monomials_map = monomials_map_new;

	return *this;
}

void NNPolynomial::mul_assign(const NNPolynomial & polynomial, shared_ptr<NNMonomial> &temp_mono)
{
	if((monomials_map.size() == 0) || (polynomial.monomials_map.size() == 0))
	{
		this->clear();
		return;
	}

	map<string, shared_ptr<NNMonomial>> monomials_map_new;

	string mono_string;
	
	for(auto iterA = monomials_map.begin(); iterA != monomials_map.end(); ++iterA){
	        for(auto iterB = polynomial.monomials_map.begin(); iterB != polynomial.monomials_map.end(); ++iterB){

		        (iterA->second)->mul(temp_mono, (iterB->second));

			temp_mono->toStringNoCoef(mono_string, dnn::curAugmentedVarNames);

			if(monomials_map_new.find(mono_string) == monomials_map_new.end()){
			        monomials_map_new[mono_string] = make_shared<NNMonomial>(temp_mono);
			}
			else{
			        (*monomials_map_new[mono_string]) += (*temp_mono);
			}

		}
	}

	monomials_map = monomials_map_new;
}

NNPolynomial & NNPolynomial::operator *= (const Interval & I)
{

	if(I.subsetZero())	// the interval is zero
	{
		clear();
	}
	else
	{
		for(auto i=monomials_map.begin(); i != monomials_map.end(); ++i){
		        (*monomials_map[i->first]->coefficient) *= I;
		}
	}

	return *this;
}

NNPolynomial NNPolynomial::operator * (const NNPolynomial & polynomial) const
{
	NNPolynomial result = *this;
	result *= polynomial;
	return result;
}

NNPolynomial NNPolynomial::operator * (const Interval & I) const
{
	NNPolynomial result = *this;
	result *= I;
	return result;
}

void NNPolynomial::add_assign(const shared_ptr<NNMonomial> & monomial)
{
	bool bAdded = false;

	if(monomial->coefficient->subsetZero())
	{
		return;
	}

	string mono_string;
	monomial->toStringNoCoef(mono_string, dnn::curAugmentedVarNames);

	//NB: this is doing a shallow copy
	if(monomials_map.find(mono_string) == monomials_map.end()){
	        monomials_map[mono_string] = monomial;
	}
	else{
	        (*monomials_map[mono_string]) += (*monomial);
	}
}

void NNPolynomial::add_assign(const shared_ptr<NNMonomial> & monomial, const std::vector<std::string> & varNames)
{
	bool bAdded = false;

	if(monomial->coefficient->subsetZero())
	{
		return;
	}

	string mono_string;
	monomial->toStringNoCoef(mono_string, varNames);

	//NB: this is doing a shallow copy
	if(monomials_map.find(mono_string) == monomials_map.end()){
	        monomials_map[mono_string] = monomial;
	}
	else{
	        (*monomials_map[mono_string]) += (*monomial);
	}
}

void NNPolynomial::add_assign(const std::list<std::shared_ptr<NNMonomial>> & monomials_other)
{

        if(monomials_other.size() == 0) return;

	for(auto iter = monomials_other.begin(); iter != monomials_other.end(); iter++){
	        add_assign(*iter);
	}
}

void NNPolynomial::add_assign(const NNPolynomial &polynomial)
{
	if(polynomial.monomials_map.size() == 0) return;
  
	for(auto iterM = polynomial.monomials_map.begin(); iterM != polynomial.monomials_map.end(); iterM++){
	  
	        string mono_string;
		(iterM->second)->toStringNoCoef(mono_string, dnn::curAugmentedVarNames);

		//NB: this is doing a shallow copy
		if(monomials_map.find(mono_string) == monomials_map.end()){
		        monomials_map[mono_string] = iterM->second;
		}
		else{
		        (*monomials_map[mono_string]) += (*iterM->second);
		}
	}	
}

void NNPolynomial::mul(NNPolynomial & result, const Interval & I) const
{
        result = NNPolynomial(*this);
	result.mul_assign(I);
}

void NNPolynomial::mul(NNPolynomial result, const int varIndex, const int degree) const
{
        result = NNPolynomial(*this);
	result.mul_assign(varIndex, degree);
}

void NNPolynomial::mul_assign(const Interval & I)
{

	if(I.subsetZero())	// the interval is zero
	{
		clear();
	}
	else
	{
		for(auto i=monomials_map.begin(); i != monomials_map.end(); ++i){
		        (*(i->second)->coefficient) *= I;
		}
	}
}

void NNPolynomial::mul_assign(const int varIndex, const int degree)
{

	for(auto i=monomials_map.begin(); i != monomials_map.end(); ++i){
	        (i->second)->addDegree(varIndex, degree);
	}
}

void NNPolynomial::inv(NNPolynomial & result) const
{
	result = *this;
	result.inv_assign();
}

void NNPolynomial::inv_assign()
{

	for(auto i=monomials_map.begin(); i != monomials_map.end(); ++i){
	        (i->second)->coefficient->inv_assign();
	}
}

void NNPolynomial::intEval(Interval & result, const std::vector<Interval> & domain) const
{
        shared_ptr<NNHornerForm> hf(new NNHornerForm());

	toHornerForm(hf);
	hf->intEval(result, domain);

}

void NNPolynomial::intEvalNormal(Interval & result, const std::vector<Interval> & step_exp_table) const
{
        result.set(0.0, 0.0);

	Interval intTemp;
	
	for(auto i=monomials_map.begin(); i != monomials_map.end(); i++){
	  
	        (i->second)->intEvalNormal(intTemp, step_exp_table);
		result += intTemp;
	}
}

void NNPolynomial::intEvalNormalTemp(Interval & result, Interval & intTemp, const std::vector<Interval> & step_exp_table) const
{
        result.set(0.0, 0.0);

	for(auto i=monomials_map.begin(); i != monomials_map.end(); i++){
	  
	        (i->second)->intEvalNormal(intTemp, step_exp_table);
		result += intTemp;
	}
}

// temp is just a placeholder, in order to avoid creating unnecessary Interval objects
void NNPolynomial::get_remainder_from_mul(Interval &result_rem, const Interval &remainder,
					  const Interval &remainder_other, const Interval & tmPolyRange,
					  Interval &temp_res, Interval &temp, const std::vector<Interval> & step_exp_table) const
{
	result_rem.set(0.0, 0.0);
	
	// these booleans are used to avoid unnecessary computations
	bool isZero, isZeroOther;
	if(!remainder_other.subsetZero())
	{
	        intEvalNormalTemp(temp_res, temp, step_exp_table);
		temp_res *= remainder_other;
		result_rem += temp_res;
		isZeroOther = false;
	}
	else{
	        isZeroOther = true;
	}

	if(!remainder.subsetZero())
	{
	        //having these two lines avoids creating a new interval
	        //to compute temp = tmPolyRange * remainder
	        temp.set(tmPolyRange);
		temp *= remainder;

		result_rem += temp;
		isZero = false;		
	}
	else{
	        isZero = true;
	}

	if(!isZero && !isZeroOther){
	        temp.set(remainder);
		temp *= remainder_other;
		
		result_rem += temp;
	}
}

void NNPolynomial::toHornerForm(shared_ptr<NNHornerForm> & result) const
{

        result->clear();

	if(monomials_map.size() == 0)
		return;	

	int numVars = (monomials_map.begin()->second)->getNumVars();

	list<shared_ptr<NNMonomial>> lstMono;

	for(auto iter_mono = monomials_map.begin(); iter_mono != monomials_map.end(); iter_mono++){
	        shared_ptr<NNMonomial> p(new NNMonomial(iter_mono->second, true)); //shallow copy
	        lstMono.push_back(p);
	}

        list<shared_ptr<NNMonomial>>::const_iterator iter = lstMono.begin();

	if((*iter)->degree() == 0)
	{
	        result->constant = make_shared<Interval>(*(*iter)->coefficient);
		iter = lstMono.erase(iter);

		if(lstMono.size() == 0)
			return;
	}

	vector<list<shared_ptr<NNMonomial>> > vlMono;
	
	for(int i=0; i<numVars; ++i)
	{
	        list<shared_ptr<NNMonomial>> lst_ith;
		vlMono.push_back(lst_ith);

		for(iter = lstMono.begin(); iter != lstMono.end();)
		{
		        if((*iter)->getDegree(i) > 0)
			{
			        (*iter)->addDegree(i, -1);
			        vlMono[i].push_back(*iter);
				iter = lstMono.erase(iter);
			}
			else
			{
				++iter;
			}
		}
	}	

	for(int i=0; i<numVars; ++i)
	{
		NNPolynomial polyTemp(vlMono[i]);
		shared_ptr<NNHornerForm> hf(new NNHornerForm());
		polyTemp.toHornerForm(hf);
		result->hornerForms.push_back(hf);
	}
}

void NNPolynomial::ctrunc(Interval & remainder, const vector<Interval> & step_exp_table, const int order)
{
	NNPolynomial polyTemp;
	string mono_string;
	
	for(auto i=monomials_map.begin(); i != monomials_map.end(); ){
	        if((i->second)->degree() > order){

			(i->second)->toStringNoCoef(mono_string, dnn::curAugmentedVarNames);
			polyTemp.monomials_map[mono_string] = i->second;
			i = monomials_map.erase(i);
		}
		else{
		        i++;
		}
	}

	polyTemp.intEvalNormal(remainder, step_exp_table);
}

void NNPolynomial::cutoff(Interval & intRem, const vector<Interval> & step_exp_table, const Interval & cutoff_threshold)
{
	NNPolynomial polyTemp;

	for(auto i=monomials_map.begin(); i != monomials_map.end(); ){
	        shared_ptr<NNMonomial> monoTemp;
		int res = (i->second)->cutoff(monoTemp, cutoff_threshold);

		switch(res)
		{
		case 0:
			++i;
			break;
		case 1:
		        polyTemp.add_assign(monoTemp);
			++i;
			break;
		case 2:
		        polyTemp.add_assign(i->second);
			i = monomials_map.erase(i);
			break;
		}		
	}
	polyTemp.intEvalNormal(intRem, step_exp_table);
}

void NNPolynomial::constant(Interval & result) const
{

	string empty = "";
	if(monomials_map.find(empty) == monomials_map.end()){
	        result.set(0.0, 0.0);	  
	}
	else{
	        result = (monomials_map.find(empty)->second)->getCoefficient();
	}
}

void NNPolynomial::rmConstant()
{

	string empty = "";
	if(monomials_map.find(empty) != monomials_map.end())
	      monomials_map.erase(empty);
}

void NNPolynomial::addConstant(const Interval & constant, const int num_vars)
{
        string empty = "";
	if(monomials_map.find(empty) == monomials_map.end()){
	        shared_ptr<NNMonomial> new_mono(new NNMonomial(constant, num_vars));
		monomials_map[empty] = new_mono;
	}
	else{
	        (*monomials_map[empty]->coefficient) += constant;
	}
}  

void NNPolynomial::clear()
{
	monomials_map.clear();
}

void NNPolynomial::toString(string & result, const vector<string> & varNames) const
{
	std::string strPoly;

	if(monomials_map.size() == 0)
	{
		strPoly = "(0)";
		return;
	}

	strPoly += '(';

	for(auto i=monomials_map.begin(); i != monomials_map.end(); ){

	        string strTemp;
		(i->second)->toString(strTemp, varNames);
		strPoly += strTemp;
		
		i++;

		if(i == monomials_map.end()){
		        strPoly += ')';
		}
		else{
		        strPoly += ' ';
			strPoly += '+';
			strPoly += ' ';
		}
	}

	result = strPoly;
}

/*----------------------------------NNHornerForm methods----------------------------------*/

void NNHornerForm::intEval(Interval & result, const std::vector<Interval> & domain) const
{
        result.set(*constant);		

	for(int i=0; i<hornerForms.size(); ++i)
	{
		Interval intHF;
		hornerForms[i]->intEval(intHF, domain);
		intHF *= domain[i];
		result += intHF;

	}
}

void NNHornerForm::insertNN(NNTaylorModel & result, const NNTaylorModelVec & vars,
			    const std::vector<Interval> & varsPolyRange,
			    const vector<Interval> & step_exp_table, const bool ctrunc, const bool cutoff,
			    const Interval & cutoff_threshold, const int order) const
{
	int numVars = vars.tms.size() + 1; // the extra 1 is for time

	result.clear();

	if(!constant->subsetZero())
	{
		NNTaylorModel tmConstant(*constant, numVars);
		result = NNTaylorModel(tmConstant);
	}

	if(hornerForms.size() > 0)						// the first variable is t
	{
		NNTaylorModel tmTemp;

		// NB: killing all time operations since there is no time in the DNN
		
		// hornerForms[0]->insertNN(tmTemp, vars, varsPolyRange, domain, step_exp_table, cutoff_threshold, ctrunc, order);
		// tmTemp.expansion.mul_assign(0,1);			// multiplied by t
		// tmTemp.remainder *= domain[0];
		// result.add_assign(tmTemp);
		
		for(int i=1; i<hornerForms.size(); ++i)
		{
		        hornerForms[i]->insertNN(tmTemp, vars, varsPolyRange, step_exp_table,
						 ctrunc, cutoff,
						 cutoff_threshold, order);	// recursive call
			
			tmTemp.mul_insert_assign(vars.tms[i-1], varsPolyRange[i-1],
						 step_exp_table, ctrunc, cutoff, cutoff_threshold, order);
			
			result.add_assign(tmTemp);
			
		}
	}
}

void NNHornerForm::insertNN(NNPolynomial &poly, Interval &remainder,
	      const NNTaylorModelVec & vars, const vector<Interval> & varsPolyRange,
	      const vector<Interval> & step_exp_table, const bool ctrunc, const bool cutoff,
	      const Interval & cutoff_threshold, const int order) const
{

	poly.clear();

	// declaring these at the top in order to avoid allocating too many mpfr variables
	Interval new_remainder, result_rem, temp_res, temp;
	int numVars = vars.tms.size() + 1; // the extra 1 is for time
	shared_ptr<NNMonomial> mono_temp(new NNMonomial(temp, numVars));

	if(constant && !constant->subsetZero())
	{
		shared_ptr<NNMonomial> new_mono(new NNMonomial(*constant, numVars));
		
		string empty = "";
		poly.monomials_map[empty] = new_mono;
	}

	if(hornerForms.size() > 0) // the first variable is t
	{
		
		for(int i=1; i<hornerForms.size(); ++i)
		{

		        NNPolynomial new_poly;
			new_remainder.set(0.0, 0.0);
		  
		        // recursive call
		        hornerForms[i]->insertNN(new_poly, new_remainder, vars, varsPolyRange,
						 step_exp_table, ctrunc, cutoff, cutoff_threshold, order);	

			auto start = chrono::high_resolution_clock::now();
			new_poly.get_remainder_from_mul(result_rem, new_remainder, vars.tms[i-1].remainder,
			 				varsPolyRange[i-1], temp_res, temp, step_exp_table);

			auto stop = chrono::high_resolution_clock::now();
			auto duration = chrono::duration_cast<chrono::nanoseconds>(stop - start);
			NNTaylorModel::remTime += duration.count();

		        start = chrono::high_resolution_clock::now();
			new_poly.mul_assign(vars.tms[i-1].expansion, mono_temp);
			stop = chrono::high_resolution_clock::now();
			duration = chrono::duration_cast<chrono::nanoseconds>(stop - start);
			NNTaylorModel::mulTime += duration.count();

			start = chrono::high_resolution_clock::now(); 
			poly.add_assign(new_poly);
			stop = chrono::high_resolution_clock::now();
			duration = chrono::duration_cast<chrono::nanoseconds>(stop - start);
			NNTaylorModel::addAssignTime += duration.count();

			remainder += result_rem;
			
		}

		if(cutoff){
			auto start = chrono::high_resolution_clock::now(); 
			poly.cutoff(temp, step_exp_table, cutoff_threshold);
			auto stop = chrono::high_resolution_clock::now();
			auto duration = chrono::duration_cast<chrono::nanoseconds>(stop - start);
			NNTaylorModel::cutoffTime += duration.count();
			remainder += temp;
		}
		if(ctrunc){
		        poly.ctrunc(temp, step_exp_table, order);
			remainder += temp;
		}
	}

}

void NNHornerForm::clear()
{
        if(constant)
	        constant->set(0,0);
	hornerForms.clear();
}
