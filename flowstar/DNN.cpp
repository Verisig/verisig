#include "DNN.h"
#include <sstream>
#include <map>
#include <yaml-cpp/yaml.h>
#include <bits/stdc++.h> 

using namespace YAML;
using namespace flowstar;
using namespace std;

Real sigmoid(Real input){
        Real sig = Real(-1) * Real(input);
	sig.exp_assign();
	
        return Real(1) / (Real(1) + sig);
}

Real tanh(Real input){
        Real posExp = Real(input);
	posExp.exp_assign();

	Real negExp = Real(Real(-1) * input);
	negExp.exp_assign();
	
        return (posExp - negExp)/(posExp + negExp);
}

Real tanhinv(Real input){

        Real r1 = Real(1) - input;
	Real r2 = Real(1) + input;

	r1.log_assign();
	r2.log_assign();

	return Real(0.5) * (r1 - r2);
}

Real swish(Real input){
        return input * sigmoid(input);
}

Real swishTen(Real input){
        return input * sigmoid(Real(10) * input);
}

Real swishHundred(Real input){
        return input * sigmoid(Real(100) * input);
}

Real lambda(int lorder, int rorder, int order, Real input){
        if (lorder + rorder == order + 1){

	        Real sig = sigmoid(input);
		sig.pow_assign(lorder);

		Real oneMsig = Real(1) - sigmoid(input);
		oneMsig.pow_assign(rorder);

		return sig * oneMsig;
	    
	  
	        //return pow(sigmoid(input), lorder) * pow(1 - sigmoid(input), rorder);
	}

	return Real(lorder) * lambda(lorder, rorder + 1, order, input) -
	  Real(rorder) * lambda(lorder + 1, rorder, order, input);
}

Real tanhlambda(int lorder, int rorder, int rec, int order, Real input){

        Real output;

        if (rec == order - 1){

	        output = Real(1) - tanh(input) * tanh(input);
		output.pow_assign(rorder);
	  
	        //output = pow(1 - tanh(input) * tanh(input), rorder);
        
		if (lorder > 0){
		        Real recReal = tanh(input);
			recReal.pow_assign(lorder);

			output = output * recReal;
			
		        //output = output * pow(tanh(input), lorder);
		}

	}
        
	else{                
	        output =  Real(-1) * Real(2) * Real(rorder) * tanhlambda(lorder + 1, rorder, rec + 1, order, input);
        
		if (lorder > 0)
		        output = output + Real(lorder) * tanhlambda(lorder - 1, rorder + 1, rec + 1, order, input);            
	}

	return output;

}

Real tanhDer(int order, Real input){
        return tanhlambda(0, 1, 0, order, input);
}

Real sigDer(int order, Real input){
        return lambda(1, 1, order, input);
}

Real swish1stDer(Real input){
        int derOrder = 1;
	return sigmoid(input) + input * lambda(1, 1, derOrder, input);
}

Real swish2ndDer(Real input){
        int derOrder = 1;
	return Real(2) * lambda(1, 1, derOrder, input) + input * lambda(1, 1, derOrder + 1, input);
}

Real swish3rdDer(Real input){
        int derOrder = 2;
	return Real(3) * lambda(1, 1, derOrder, input) + input * lambda(1, 1, derOrder + 1, input);
}

Real swishTen1stDer(Real input){
        int derOrder = 1;
	return sigmoid(Real(10) * input) + Real(10) * input * lambda(1, 1, derOrder, Real(10) * input);
}

Real swishTen2ndDer(Real input){
        int derOrder = 1;
	return Real(20) * lambda(1, 1, derOrder, Real(10) * input) + Real(100) * input * lambda(1, 1, derOrder + 1, Real(10) * input);
}

Real swishTen3rdDer(Real input){
        int derOrder = 2;
	return Real(300) * lambda(1, 1, derOrder, Real(10) * input) + Real(1000) * input * lambda(1, 1, derOrder + 1, Real(10) * input);
}

Real swishHundred1stDer(Real input){
        int derOrder = 1;
	return sigmoid(Real(100) * input) + Real(100) * input * lambda(1, 1, derOrder, Real(100) * input);
}

Real swishHundred2ndDer(Real input){
        int derOrder = 1;
	return Real(200) * lambda(1, 1, derOrder, Real(100) * input) + Real(1000) * input * lambda(1, 1, derOrder + 1, Real(100) * input);
}

Real swishHundred3rdDer(Real input){
        int derOrder = 2;
	return Real(30000) * lambda(1, 1, derOrder, Real(100) * input) + Real(1000000) * input * lambda(1, 1, derOrder + 1, Real(100) * input);
}

Real getSig4thRemUpperBound(Interval bounds){

	Real Q;
  
        //Region 5
        if(bounds.inf() >= 3.15){
	        Q = sigDer(4, Real(bounds.sup()));
	}

	//Region 4
        else if(bounds.inf() >= 0.85 && bounds.inf() <= 3.15){

	        Q = sigDer(4, Real(bounds.inf()));

		//sup is in Region 5		
	        if(bounds.sup() >= 3.15){
		        Real temp = sigDer(4, Real(bounds.inf()));
			if(temp > Q) Q = temp;
		}
	        
	}

	//Region 3.5
	else if(bounds.inf() >= 0.83 && bounds.inf() <= 0.85){
	        Q = Real(0.1277);
	}

	//Region 3
        else if(bounds.inf() >= -0.83 && bounds.inf() <= 0.83){
	  
	        Q = sigDer(4, Real(bounds.sup()));

		//sup is in Region 4
	        if(bounds.sup() >= 0.83) Q = Real(0.1277);
	}

	//Region 2
	else if(bounds.inf() >= -3.13 && bounds.inf() <= -0.85){

	        Q = sigDer(4, Real(bounds.inf()));

		//sup is in Region 3		
	        if(bounds.sup() >= -0.85 && bounds.sup() <= 0.85){
		        Real temp = sigDer(4, Real(bounds.sup()));
			if(temp > Q) Q = temp;
		}
		//sup is beyond global max
		else if(bounds.sup() >= 0.85){
		        Q = Real(0.1277);
		}
	}

	//Region 1.5
	else if(bounds.inf() >= -3.15 && bounds.inf() <= -3.13){

	        Q = Real(0.01908);

		//sup is in Region 3
	        if(bounds.sup() >= -0.85 && bounds.sup() <= 0.85){
		        Real temp = sigDer(4, Real(bounds.sup()));
			if(temp > Q) Q = temp;
		}

		//sup is beyond global max
		else if(bounds.sup() >= 0.85){
		        Q = Real(0.1277);
		}
	}

	//Region 1
	else if(bounds.inf() <= -3.15){

	        Q = sigDer(4, Real(bounds.sup()));

		//sup is in Region 2
	        if(bounds.sup() >= -3.15 && bounds.sup() <= 0.83){

		        Q = Real(0.01908);

		        Real temp = sigDer(4, Real(bounds.sup()));
			if(temp > Q) Q = temp;
		}

		//sup is beyond global max
		else if(bounds.sup() >= 0.83){
		        Q = Real(0.1277);
		}
	}
		
        return Q; 
}

Real getSig4thRemLowerBound(Interval bounds){

	Real q;
  
        //Region 5
        if(bounds.inf() >= 3.15){
	        q = sigDer(4, Real(bounds.inf()));
	}

	//Region 4.5
	else if(bounds.inf() <= 3.13 && bounds.inf() >= 3.15){
	        q = Real(-0.01908);
	}
	
	//Region 4
        else if(bounds.inf() >= 0.85 && bounds.inf() <= 3.13){

	        q = sigDer(4, Real(bounds.sup()));

		//sup is in Region 5		
	        if(bounds.sup() >= 3.13){
		        q = Real(-0.01908);
		}
	        
	}

	//Region 3
        else if(bounds.inf() >= -0.83 && bounds.inf() <= 0.85){
	  
	        q = sigDer(4, Real(bounds.inf()));

		//sup is in Region 4
	        if(bounds.sup() >= 0.83 && bounds.sup() <= 3.13){
		        Real temp = sigDer(4, Real(bounds.sup()));
			if(q > temp) q = temp;
		}
		//sup is in Region 5
		else if(bounds.sup() >= 3.13){
		        if(q > Real(-0.01908)) q = Real(-0.01908);		  
		}
	}

	//Region 2.5
	else if(bounds.inf() <= -0.85 && bounds.inf() >= -0.83){
	        q = Real(-0.1277);
	}		

	//Region 2
	else if(bounds.inf() >= -3.15 && bounds.inf() <= -0.85){

	        q = sigDer(4, Real(bounds.sup()));

		//sup is in Region 1.5
	        if(bounds.sup() >= -3.15 && bounds.sup() <= -0.85){
		        Real temp = sigDer(4, Real(bounds.inf()));
			if(q > temp) q = temp;
		}
		
		//sup is beyond global min
		else if(bounds.sup() >= 0.85){
		        q = Real(-0.1277);
		}
	}

	//Region 1
	else if(bounds.inf() <= -3.15){

	        q = sigDer(4, Real(bounds.inf()));

		//sup is in Region 2
	        if(bounds.sup() >= -3.15 && bounds.sup() <= -0.85){

		        Real temp = sigDer(4, Real(bounds.sup()));
			if(q > temp) q = temp;
		}

		//sup is beyond global min
		else if(bounds.sup() >= -0.85){
		        q = Real(-0.1277);
		}
	}
		
        return q;  
}

Interval getSig4thDerRemBound(const Interval inputBounds, double apprPoint){
  
        Interval upper = Interval(apprPoint, inputBounds.sup());
	Interval lower = Interval(inputBounds.inf(), apprPoint);
  
        Real Q_u = getSig4thRemUpperBound(upper);
	Real Q_l = getSig4thRemUpperBound(lower);
	Real q_u = getSig4thRemLowerBound(upper);
	Real q_l = getSig4thRemLowerBound(lower);

        Real fact = Real(24);
	Real maxPosDev = Real(inputBounds.sup() - apprPoint);
	Real maxNegDev = Real(inputBounds.inf()- apprPoint);
	maxPosDev.pow_assign(4);
	maxNegDev.pow_assign(4);

	Real u = (maxPosDev * Q_u) / fact;
	Real l = (maxPosDev * q_u) / fact;

	if((maxNegDev * Q_l) / fact > u) u = (maxNegDev * Q_l) / fact;
	if(l > (maxNegDev * q_l) / fact) l = (maxNegDev * q_l) / fact;

	//these checks are necessary because the remainder is always 0 at apprPoint
	if(Real(0) > u) u = Real(0);
	if(l > Real(0)) l = Real(0);
	
        return Interval(l.getValue_RNDD(), u.getValue_RNDU());
    
}

Real getTanh4thRemUpperBound(Interval bounds){

        //first get the extrema of the 5th derivative
        //the hardcoded numbers are from a quadratic equation :)
        Real sqrt105 = Real(105);
	sqrt105.sqrt_assign();

	Real sol1 = (Real(15) + sqrt(105)) / Real(30);
	Real sol2 = (Real(15) - sqrt(105)) / Real(30);

	Real tanh1 = Real(sol1);
	tanh1.sqrt_assign();
	
	Real tanh2 = Real(sol1);
	tanh2.sqrt_assign();
	tanh2 = Real(-1) * tanh2;

	Real tanh3 = Real(sol2);
	tanh3.sqrt_assign();
	
	Real tanh4 = Real(sol2);
	tanh4.sqrt_assign();
	tanh4 = Real(-1) * tanh4;

	//NB: these are ordered in increasing order
	Real ext1 = tanhinv(tanh1);
	Real ext2 = tanhinv(tanh3);
	Real ext3 = tanhinv(tanh4);
	Real ext4 = tanhinv(tanh2);

	Real localMax = tanhDer(4, ext1);
	Real globalMax = tanhDer(4, ext3);

	Real Q;
  
        //Region 5
        if(Real(bounds.inf()) >= ext4){
	        Q = tanhDer(4, Real(bounds.sup()));
	}

	//Region 4
        else if(Real(bounds.inf()) >= ext3 && ext4 >= Real(bounds.inf())){

	        Q = tanhDer(4, Real(bounds.inf()));

		//sup is in Region 5		
	        if(Real(bounds.sup()) >= ext4){
		        Real temp = tanhDer(4, Real(bounds.inf()));
			if(temp > Q) Q = temp;
		}
	        
	}

	//Region 3
        else if(Real(bounds.inf()) >= ext2 && ext3 >= Real(bounds.inf())){
	  
	        Q = tanhDer(4, Real(bounds.sup()));

		//sup is in Region 4
	        if(Real(bounds.sup()) >= ext3) Q = Real(4.0859);
	}

	//Region 2
	else if(Real(bounds.inf()) >= ext1 && ext2 >= Real(bounds.inf())){

	        Q = tanhDer(4, Real(bounds.inf()));

		//sup is in Region 3		
	        if(Real(bounds.sup()) >= ext2 && ext3 >= Real(bounds.sup())){
		        Real temp = tanhDer(4, Real(bounds.sup()));
			if(temp > Q) Q = temp;
		}
		//sup is beyond global max
		else if(Real(bounds.sup()) >= ext3){
		        Q = globalMax;
		}
	}

	//Region 1
	else if(ext1 >= Real(bounds.inf())){

	        Q = tanhDer(4, Real(bounds.sup()));

		//sup is in Region 2
	        if(Real(bounds.sup()) >= ext1 && ext3 >= Real(bounds.sup())){

		        Q = localMax;

		        Real temp = tanhDer(4, Real(bounds.sup()));
			if(temp > Q) Q = temp;
		}

		//sup is beyond global max
		else if(Real(bounds.sup()) >= ext3){
		        Q = globalMax;
		}
	}
		
        return Q;  
}

Real getTanh4thRemLowerBound(Interval bounds){

        //first get the extrema of the 5th derivative
        //the hardcoded numbers are from a quadratic equation :)
        Real sqrt105 = Real(105);
	sqrt105.sqrt_assign();

	Real sol1 = (Real(15) + sqrt(105)) / Real(30);
	Real sol2 = (Real(15) - sqrt(105)) / Real(30);

	Real tanh1 = Real(sol1);
	tanh1.sqrt_assign();
	
	Real tanh2 = Real(sol1);
	tanh2.sqrt_assign();
	tanh2 = Real(-1) * tanh2;

	Real tanh3 = Real(sol2);
	tanh3.sqrt_assign();
	
	Real tanh4 = Real(sol2);
	tanh4.sqrt_assign();
	tanh4 = Real(-1) * tanh4;

	//NB: these are ordered in increasing order
	Real ext1 = tanhinv(tanh1);
	Real ext2 = tanhinv(tanh3);
	Real ext3 = tanhinv(tanh4);
	Real ext4 = tanhinv(tanh2);

	Real localMin = tanhDer(4, ext4);
	Real globalMin = tanhDer(4, ext2);

	Real q;
  
        //Region 5
        if(Real(bounds.inf()) >= ext4){
	        q = tanhDer(4, Real(bounds.inf()));
	}
	
	//Region 4
        else if(Real(bounds.inf()) >= ext3 && ext4 >= Real(bounds.inf())){

	        q = tanhDer(4, Real(bounds.sup()));

		//sup is in Region 5		
	        if(Real(bounds.sup()) >= ext4){
		        q = localMin;
		}
	        
	}

	//Region 3
        else if(Real(bounds.inf()) >= ext2 && ext3 >= Real(bounds.inf())){
	  
	        q = tanhDer(4, Real(bounds.inf()));

		//sup is in Region 4
	        if(Real(bounds.sup()) >= ext3 && ext4 >= Real(bounds.sup())){
		        Real temp = tanhDer(4, Real(bounds.sup()));
			if(q > temp) q = temp;
		}
		
		//sup is in Region 5
		else if(Real(bounds.sup()) >= ext4){
		        if(q > localMin) q = localMin;
		}
	}	

	//Region 2
	else if(Real(bounds.inf()) >= ext1 && ext2 >= Real(bounds.inf())){

	        q = tanhDer(4, Real(bounds.sup()));
		
		//sup is beyond global min
		if(Real(bounds.sup()) >= ext2){
		        q = globalMin;
		}
	}

	//Region 1
	else if(ext1 >= Real(bounds.inf())){

	        q = tanhDer(4, Real(bounds.inf()));

		//sup is in Region 2
	        if(Real(bounds.sup()) >= ext1 && ext2 >= Real(bounds.sup())){
		        Real temp = tanhDer(4, Real(bounds.sup()));
			if(q > temp) q = temp;
		}

		//sup is beyond global min
		else if(Real(bounds.sup()) >= ext2){
		        q = globalMin;
		}
	}
		
        return q;  
}

Interval getTanh4thDerRemBound(const Interval inputBounds, double apprPoint){

        Interval upper = Interval(apprPoint, inputBounds.sup());
	Interval lower = Interval(inputBounds.inf(), apprPoint);
  
        Real Q_u = getTanh4thRemUpperBound(upper);
	Real Q_l = getTanh4thRemUpperBound(lower);
	Real q_u = getTanh4thRemLowerBound(upper);
	Real q_l = getTanh4thRemLowerBound(lower);

        Real fact = Real(24);
	Real maxPosDev = Real(inputBounds.sup() - apprPoint);
	Real maxNegDev = Real(inputBounds.inf()- apprPoint);
	maxPosDev.pow_assign(4);
	maxNegDev.pow_assign(4);

	Real u = (maxPosDev * Q_u) / fact;
	Real l = (maxPosDev * q_u) / fact;

	if((maxNegDev * Q_l) / fact > u) u = (maxNegDev * Q_l) / fact;
	if(l > (maxNegDev * q_l) / fact) l = (maxNegDev * q_l) / fact;

	//these checks are necessary because the remainder is always 0 at apprPoint
	if(Real(0) > u) u = Real(0);
	if(l > Real(0)) l = Real(0);
	
        return Interval(l.getValue_RNDD(), u.getValue_RNDU());
    
}

Real getSwish4thDerBound(Interval bounds){

        if(bounds.inf() <= 2.2 && bounds.sup() >= -2.2){
	        return Real(0.13);
	}

        if(bounds.inf() > 2.2 && bounds.inf() <= 6){
	        return Real(0.04);
	}
	
        if(bounds.sup() >= -6 && bounds.sup() < -2.2){
	        return Real(0.04);
	}

	if(bounds.inf() > 6){
	        return Real(0.013);
	}

	if(bounds.sup() < -6){
	        return Real(0.013);
	}		
  
        return Real(0.13);
}

Real getSwish3rdDerBound(Interval bounds){

        if(bounds.inf() <= 3 && bounds.sup() >= -3){
	        return Real(0.31);
	}

        if(bounds.inf() > 3 && bounds.inf() <= 5){
	        return Real(0.025);
	}

        if(bounds.sup() >= -5 && bounds.sup() < -3){
	        return Real(0.025);
	}			

	if(bounds.inf() > 5){
	        return Real(0.013);
	}

	if(bounds.sup() < -5){
	        return Real(0.013);
	}		
  
        return Real(0.31);
}

Real getSwishTen4thDerBound(Interval bounds){

        if(bounds.inf() <= 0.071 && bounds.sup() >= -0.071){
	        return Real(500);
	}

        if(bounds.inf() > 0.071 && bounds.inf() <= 0.4){
	        return Real(204);
	}

	if(bounds.inf() > 0.4 && bounds.inf() <= 1.2){
	        return Real(10);
	}

        if(bounds.sup() >= -0.4 && bounds.sup() < -0.071){
	        return Real(204);
	}

	if(bounds.sup() >= -1.2 && bounds.sup() < -0.4){
	        return Real(10);
	}	

	if(bounds.inf() > 1.2){
	        return Real(0.05);
	}

	if(bounds.sup() < -1.2){
	        return Real(0.05);
	}		
  
        return Real(500);
}

Real getSwishTen3rdDerBound(Interval bounds){

        if(bounds.inf() <= 0.3 && bounds.sup() >= -0.3){
	        return Real(30.9);
	}

        if(bounds.inf() > 0.3 && bounds.inf() <= 0.91){
	        return Real(2.6);
	}

        if(bounds.sup() >= -0.91 && bounds.sup() < -0.3){
	        return Real(2.6);
	}			

	if(bounds.inf() > 0.91){
	        return Real(0.07);
	}

	if(bounds.sup() < -0.91){
	        return Real(0.07);
	}		
  
        return Real(30.9);
}

Real getSwishHundred4thDerBound(Interval bounds){

        if(bounds.inf() <= 0.007 && bounds.sup() >= -0.007){
	        return Real(500000);
	}

        if(bounds.inf() > 0.007 && bounds.inf() <= 0.045){
	        return Real(200000);
	}

        if(bounds.inf() >= -0.045 && bounds.inf() <= -0.007){
	        return Real(200000);
	}

	if(bounds.inf() > 0.045 && bounds.inf() <= 0.12){
	        return Real(4700);
	}

        if(bounds.sup() >= -0.12 && bounds.sup() < -0.045){
	        return Real(4700);
	}

	if(bounds.sup() >= 0.12 && bounds.sup() < 0.17){
	        return Real(50);
	}

	if(bounds.sup() >= -0.17 && bounds.sup() < -0.12){
	        return Real(50);
	}	

	if(bounds.inf() > 0.17){
	        return Real(0.55);
	}

	if(bounds.sup() < -0.17){
	        return Real(0.55);
	}		
  
        return Real(500000);
}

Real getSwishHundred3rdDerBound(Interval bounds){

        if(bounds.inf() <= 0.2 && bounds.sup() >= -0.2){
	        return Real(7400);
	}

        if(bounds.inf() > 0.2 && bounds.inf() <= 0.5){
	        return Real(8800);
	}

        if(bounds.sup() >= -0.5 && bounds.sup() < -0.2){
	        return Real(8800);
	}

        if(bounds.inf() > 0.5 && bounds.inf() <= 0.8){
	        return Real(3000);
	}

        if(bounds.sup() >= -0.8 && bounds.sup() < -0.5){
	        return Real(3000);
	}

        if(bounds.inf() > 0.8 && bounds.inf() <= 1.2){
	        return Real(260);
	}

        if(bounds.sup() >= -1.2 && bounds.sup() < -0.8){
	        return Real(260);
	}	

	if(bounds.inf() > 1.2){
	        return Real(7.5);
	}

	if(bounds.sup() < -1.2){
	        return Real(7.5);
	}		
  
        return Real(8800);
}

void dnn::sig_reset(TaylorModel &tmReset, const Interval intC, const int varInd, const int numVars){

    double bestPoint = intC.midpoint();
    Interval bestRem = getSig4thDerRemBound(intC, intC.midpoint());

    Real apprPoint = Real(bestPoint);
    Real appr0 = sigmoid(apprPoint);

    Real coef1 = sigDer(1, apprPoint);
    Real coef2 = sigDer(2, apprPoint)/2;
    Real coef3 = sigDer(3, apprPoint)/6;						

    Interval deg0Int = Interval(appr0);

    Interval deg1Int = Interval(coef1);
    Interval deg2Int = Interval(coef2);
    Interval deg3Int = Interval(coef3);

    std::vector<int> deg1(numVars, 0);
    deg1[varInd + 1] = 1;
    std::vector<int> deg2(numVars, 0);
    deg2[varInd + 1] = 2;
    std::vector<int> deg3(numVars, 0);
    deg3[varInd + 1] = 3;
    
    Polynomial deg0Poly = Polynomial(Monomial(deg0Int, numVars));
    
    Polynomial deg1Poly = Polynomial(Monomial(Interval(Real(-1) * coef1 * apprPoint), numVars)) +
      Polynomial(Monomial(deg1Int, deg1));

    Polynomial deg2Poly = Polynomial(Monomial(Interval(coef2 * apprPoint * apprPoint), numVars)) -
      Polynomial(Monomial(Interval(Real(2) * coef2 * apprPoint), deg1)) +
      Polynomial(Monomial(deg2Int, deg2));

    Polynomial deg3Poly = Polynomial(Monomial(Interval(Real(-1) * coef3 * apprPoint * apprPoint * apprPoint),
						  numVars)) +
	  Polynomial(Monomial(Interval(Real(3) * coef3 * apprPoint * apprPoint), deg1)) -
	  Polynomial(Monomial(Interval(Real(3) * coef3 * apprPoint), deg2)) +
	  Polynomial(Monomial(deg3Int, deg3));
    
    Polynomial exp = deg0Poly + deg1Poly + deg2Poly + deg3Poly;

    //printf("upper: %13.10f\n", rem.sup());
    //printf("lower: %13.10f\n", rem.inf());

    tmReset.expansion = exp;
    tmReset.remainder = bestRem;
}

void dnn::swish_reset(TaylorModel &tmReset, const Interval intC, const int varInd, const int numVars){

    Real midPoint = Real(intC.midpoint());  
    Real apprPoint = swish(midPoint);
						
    //First try a 2nd order TS approximation
    
    Real coef1 = swish1stDer(midPoint);
    Real coef2 = swish2ndDer(midPoint)/2;
    Real coef3 = swish3rdDer(midPoint)/6;
						

    Real derBound = getSwish3rdDerBound(intC);
						
    Real maxDev = Real(intC.sup()) - midPoint;
    if (midPoint - Real(intC.inf()) > maxDev){
        maxDev = midPoint - Real(intC.inf());
    }
						
    Real fact = Real(6);					       
    maxDev.pow_assign(3);
    Real remainder = (derBound * maxDev) / fact;

    Interval apprInt = Interval(apprPoint);

    Interval deg1Int = Interval(coef1);
    Interval deg2Int = Interval(coef2);
    Interval deg3Int = Interval(coef3);

    std::vector<int> deg1(numVars, 0);
    deg1[varInd + 1] = 1;
    std::vector<int> deg2(numVars, 0);
    deg2[varInd + 1] = 2;
    std::vector<int> deg3(numVars, 0);
    deg3[varInd + 1] = 3;
						
    Polynomial deg0Poly = Polynomial(Monomial(apprInt, numVars));
						
    Polynomial deg1Poly = Polynomial(Monomial(Interval(Real(-1) * coef1 * midPoint), numVars)) +
      Polynomial(Monomial(deg1Int, deg1));

    Polynomial deg2Poly = Polynomial(Monomial(Interval(coef2 * midPoint * midPoint), numVars)) -
      Polynomial(Monomial(Interval(Real(2) * coef2 * midPoint), deg1)) +
      Polynomial(Monomial(deg2Int, deg2));
						
    Polynomial exp = deg0Poly + deg1Poly + deg2Poly;
    Interval rem;
    remainder.to_sym_int(rem);

    //if uncertainty too large, use a 3rd order approximation
    if (rem.width() > 0.00001){
        fact = 24;
	maxDev = Real(intC.sup()) - midPoint;
	maxDev.pow_assign(4);

	derBound = getSwish4thDerBound(intC);
	
	remainder = (derBound * maxDev) / fact;
	remainder.to_sym_int(rem);

	Polynomial deg3Poly = Polynomial(Monomial(Interval(Real(-1) * coef3 * midPoint * midPoint * midPoint), numVars)) +
	  Polynomial(Monomial(Interval(Real(3) * coef3 * midPoint * midPoint), deg1)) -
	  Polynomial(Monomial(Interval(Real(3) * coef3 * midPoint), deg2)) + Polynomial(Monomial(deg3Int, deg3));
						
	exp += deg3Poly;
    }

    tmReset.expansion = exp;
    tmReset.remainder = rem;
    
}

void dnn::swish10_reset(TaylorModel &tmReset, const Interval intC, const int varInd, const int numVars){

    Real midPoint = Real(intC.midpoint());  
    Real apprPoint = swishTen(midPoint);
						
    //First try a 2nd order TS approximation
    
    Real coef1 = swishTen1stDer(midPoint);
    Real coef2 = swishTen2ndDer(midPoint)/2;
    Real coef3 = swishTen3rdDer(midPoint)/6;
						

    Real derBound = getSwishTen3rdDerBound(intC);
						
    Real maxDev = Real(intC.sup()) - midPoint;
    if (midPoint - Real(intC.inf()) > maxDev){
        maxDev = midPoint - Real(intC.inf());
    }
						
    Real fact = Real(6);					       
    maxDev.pow_assign(3);
    Real remainder = (derBound * maxDev) / fact;

    Interval apprInt = Interval(apprPoint);

    Interval deg1Int = Interval(coef1);
    Interval deg2Int = Interval(coef2);
    Interval deg3Int = Interval(coef3);

    std::vector<int> deg1(numVars, 0);
    deg1[varInd + 1] = 1;
    std::vector<int> deg2(numVars, 0);
    deg2[varInd + 1] = 2;
    std::vector<int> deg3(numVars, 0);
    deg3[varInd + 1] = 3;
						
    Polynomial deg0Poly = Polynomial(Monomial(apprInt, numVars));
						
    Polynomial deg1Poly = Polynomial(Monomial(Interval(Real(-1) * coef1 * midPoint), numVars)) +
      Polynomial(Monomial(deg1Int, deg1));

    Polynomial deg2Poly = Polynomial(Monomial(Interval(coef2 * midPoint * midPoint), numVars)) -
      Polynomial(Monomial(Interval(Real(2) * coef2 * midPoint), deg1)) +
      Polynomial(Monomial(deg2Int, deg2));
						
    Polynomial exp = deg0Poly + deg1Poly + deg2Poly;
    Interval rem;
    remainder.to_sym_int(rem);

    //if uncertainty too large, use a 3rd order approximation
    if (rem.width() > 0.00001){
        fact = 24;
	maxDev = Real(intC.sup()) - midPoint;
	maxDev.pow_assign(4);

	derBound = getSwishTen4thDerBound(intC);
	
	remainder = (derBound * maxDev) / fact;
	remainder.to_sym_int(rem);

	Polynomial deg3Poly = Polynomial(Monomial(Interval(Real(-1) * coef3 * midPoint * midPoint * midPoint), numVars)) +
	  Polynomial(Monomial(Interval(Real(3) * coef3 * midPoint * midPoint), deg1)) -
	  Polynomial(Monomial(Interval(Real(3) * coef3 * midPoint), deg2)) + Polynomial(Monomial(deg3Int, deg3));
						
	exp += deg3Poly;
    }

    tmReset.expansion = exp;
    tmReset.remainder = rem;
    
}

void dnn::tanh_reset(TaylorModel &tmReset, const Interval intC, const int varInd, const int numVars){

    double bestPoint = intC.midpoint();
    Interval bestRem = getTanh4thDerRemBound(intC, intC.midpoint());

    Real apprPoint = Real(bestPoint);
    Real appr0 = tanh(apprPoint);

    Real coef1 = tanhDer(1, apprPoint);
    Real coef2 = tanhDer(2, apprPoint)/2;
    Real coef3 = tanhDer(3, apprPoint)/6;						

    Interval deg0Int = Interval(appr0);

    Interval deg1Int = Interval(coef1);
    Interval deg2Int = Interval(coef2);
    Interval deg3Int = Interval(coef3);

    std::vector<int> deg1(numVars, 0);
    deg1[varInd + 1] = 1;
    std::vector<int> deg2(numVars, 0);
    deg2[varInd + 1] = 2;
    std::vector<int> deg3(numVars, 0);
    deg3[varInd + 1] = 3;
    
    Polynomial deg0Poly = Polynomial(Monomial(deg0Int, numVars));
    
    Polynomial deg1Poly = Polynomial(Monomial(Interval(Real(-1) * coef1 * apprPoint), numVars)) +
      Polynomial(Monomial(deg1Int, deg1));

    Polynomial deg2Poly = Polynomial(Monomial(Interval(coef2 * apprPoint * apprPoint), numVars)) -
      Polynomial(Monomial(Interval(Real(2) * coef2 * apprPoint), deg1)) +
      Polynomial(Monomial(deg2Int, deg2));

    Polynomial deg3Poly = Polynomial(Monomial(Interval(Real(-1) * coef3 * apprPoint * apprPoint * apprPoint),
						  numVars)) +
      Polynomial(Monomial(Interval(Real(3) * coef3 * apprPoint * apprPoint), deg1)) -
      Polynomial(Monomial(Interval(Real(3) * coef3 * apprPoint), deg2)) +
      Polynomial(Monomial(deg3Int, deg3));
    
    Polynomial exp = deg0Poly + deg1Poly + deg2Poly + deg3Poly;

    tmReset.expansion = exp;
    tmReset.remainder = bestRem;
}

void dnn::relu_reset(TaylorModel &tmReset, const Interval intC, const int varInd, const int numVars){


    Polynomial exp;
    Interval rem;

    //ReLU all in 0 area
    if(intC.sup() < 0){

        Interval zeroInt = Interval(0.0, 0.0);
						
	exp = Polynomial(Monomial(zeroInt, numVars));
	rem = zeroInt;
    }
    
    //ReLU all in positive area
    else if(intC.inf() > 0){

        std::vector<int> deg1(numVars, 0);
	deg1[varInd + 1] = 1;

	Polynomial deg1Poly = Polynomial(Monomial(Interval(1.0, 1.0), deg1));

	exp = deg1Poly;
	rem = Interval(0.0, 0.0);

    }					

    //ReLU in both areas: approximate using the Swish function
    else{
        Real midPoint = Real(intC.midpoint());
	
        Interval apprPoint = swishHundred(midPoint);
	
	//NB: This assumes a 2nd order TS approximation

	Real coef1 = swishHundred1stDer(midPoint);
	Real coef2 = swishHundred2ndDer(midPoint)/2;
	Real coef3 = swishHundred3rdDer(midPoint)/6;
						
	Real derBound = getSwishHundred3rdDerBound(intC);

	Real maxDev = Real(intC.sup()) - midPoint;
	if (midPoint - Real(intC.inf()) > maxDev){
	    maxDev = midPoint - Real(intC.inf());
	}
						
	Real fact = Real(6);					       
	maxDev.pow_assign(3);
	Real remainder = (derBound * maxDev) / fact;

	Interval apprInt = Interval(apprPoint);

	Interval deg1Int = Interval(coef1);
	Interval deg2Int = Interval(coef2);
	Interval deg3Int = Interval(coef3);
	
	std::vector<int> deg1(numVars, 0);
	deg1[varInd + 1] = 1;
	std::vector<int> deg2(numVars, 0);
	deg2[varInd + 1] = 2;
	std::vector<int> deg3(numVars, 0);
	deg3[varInd + 1] = 3;
						
	Polynomial deg0Poly = Polynomial(Monomial(apprInt, numVars));
	
	Polynomial deg1Poly = Polynomial(Monomial(Interval(Real(-1) * coef1 * midPoint), numVars)) +
	  Polynomial(Monomial(deg1Int, deg1));
	
	Polynomial deg2Poly = Polynomial(Monomial(Interval(coef2 * midPoint * midPoint), numVars)) -
	  Polynomial(Monomial(Interval(Real(2) * coef2 * midPoint), deg1)) +
	  Polynomial(Monomial(deg2Int, deg2));
	
	exp = deg0Poly + deg1Poly + deg2Poly;
	remainder.to_sym_int(rem);

	//if uncertainty too large, use a 3rd order approximation
	if (rem.width() > 0.00001){
	    fact = 24;
	    maxDev = Real(intC.sup()) - midPoint;
	    maxDev.pow_assign(4);
	    
	    derBound = getSwishHundred4thDerBound(intC);
	    
	    remainder = (derBound * maxDev) / fact;
	    remainder.to_sym_int(rem);
	    
	    Polynomial deg3Poly = Polynomial(Monomial(Interval(Real(-1) * coef3 * midPoint * midPoint * midPoint), numVars)) +
	      Polynomial(Monomial(Interval(Real(3) * coef3 * midPoint * midPoint), deg1)) -
	      Polynomial(Monomial(Interval(Real(3) * coef3 * midPoint), deg2)) +
	      Polynomial(Monomial(deg3Int, deg3));

	    
	    exp += deg3Poly;
	}
	
	//Add approximation error between ReLU and Swish
	// exp += Polynomial(Monomial(Interval(0.0014, 0.0014), numVars));
	// rem += Interval(-0.0014, 0.0014);
	
	//Add approximation error between ReLU and Swish
	
	if (intC.inf() > -0.0127){
	    double maxDer = -swishHundred(Real(intC.inf())).getValue_RNDD();
	  
	    if(intC.sup() - swishHundred(Real(intC.sup())).getValue_RNDD() > maxDer){
	        maxDer = intC.sup() - swishHundred(Real(intC.sup())).getValue_RNDD();
	    }
	    
	    rem += Interval(0.0, maxDer);
	}
	else{						  
	    rem += Interval(0.0, 0.0028);
	}
	printf("input interval: [%f, %f]\n", intC.inf(), intC.sup());
	printf("returned remainder: [%f, %f]\n", rem.inf(), rem.sup());
    }

    tmReset.expansion = exp;
    tmReset.remainder = rem;    
  
}

void dnn::convert_TM_dimension(TaylorModel &tmNew, const TaylorModel &tmOld, int dimNew){

    Polynomial newPoly;

    std::list<Monomial>::const_iterator it = tmOld.expansion.monomials.begin();

    for (; it != tmOld.expansion.monomials.end(); ++it){
	
        std::vector<int> newDegs(dimNew, 0);

	Monomial curM = *it;

	for(int j = 0; j < curM.getDegrees().size(); j++){

	    if(j >= dimNew){
	      
	        if(curM.getDegrees()[j] != 0){
		    printf("Error while adding DNN states to reachability problem. Exiting...\n");
		    exit(-1);
		}
	    }

	    else{
	      
	        newDegs[j] = curM.getDegrees()[j];
		
	    }
	}

	newPoly.monomials.push_back(Monomial(curM.getCoefficient(), newDegs));
      
    }

    tmNew.expansion = Polynomial(newPoly);
    tmNew.remainder = Interval(tmOld.remainder);
  
}

void dnn::load_dnn(std::vector<ResetMap> &resets, std::vector<dnn::activation> &activations,
		  Variables &augmentedStateVars, std::vector<std::string> &augmentedVarNames,
		  const Variables &vars, const std::vector<std::string> &stateVarNames, const std::string filename) {
    Node dnn = LoadFile(filename);

    Node weights = dnn["weights"];
    Node offsets = dnn["offsets"];
    Node acts = dnn["activations"];

    //get number of neural states
    int numNeurStates = offsets[1].size();
    for(int i=0; i<offsets.size(); i++) {
            if(offsets[i+1].size() > numNeurStates) numNeurStates = offsets[i+1].size();
    }

    std::unordered_set <std::string> names_set;
    
    //add all DNN variables
    augmentedVarNames.push_back("local_t");
    augmentedStateVars.declareVar("local_t");
    
    for(int varInd = 0; varInd < stateVarNames.size(); varInd ++){
      
        augmentedVarNames.push_back(stateVarNames[varInd]);
	augmentedStateVars.declareVar(stateVarNames[varInd]);
	names_set.insert(stateVarNames[varInd]);
    }
						
    for(int varInd = 0; varInd < numNeurStates; varInd ++){
        if(names_set.find("_f" + std::to_string(varInd+1)) == names_set.end()){
	    augmentedVarNames.push_back("_f" + std::to_string(varInd+1));
	    augmentedStateVars.declareVar("_f" + std::to_string(varInd+1));
	}
    }
    
    for(int i=0; i<weights.size(); i++) {
  
	int layerId = i+1;
	std::string activationFcn = acts[layerId].as<string>();
	dnn::activation activationAsEnum = LINEAR;

	if(!strncmp(activationFcn.c_str(), "Tanh", strlen("Tanh"))){
	        activationAsEnum = dnn::TANH;
	}
	else if(!strncmp(activationFcn.c_str(), "Sigmoid", strlen("Sigmoid"))){
	        activationAsEnum = dnn::SIGMOID;
	}
	else if(!strncmp(activationFcn.c_str(), "Swish", strlen("Swish"))){
	        activationAsEnum = dnn::SWISH;
	}	
	else if(!strncmp(activationFcn.c_str(), "Relu", strlen("Relu"))){
	        activationAsEnum = dnn::RELU;
	}
	

        map<string,TaylorModel> taylorModels;
        Node layer = weights[layerId];
        int layerSize = layer[0].size();

        int currentNeuron = 1;
        for(const_iterator neuronIt=layer.begin(); neuronIt != layer.end(); ++neuronIt) {
            int currentWeight = 1;
            stringstream buffer;
            for(int i = 0; i < layerSize; i++) {

                buffer << (*neuronIt)[i];
                buffer << " * " << "_f" << currentWeight << + " + ";

                currentWeight++;
            }

            buffer << offsets[layerId][currentNeuron-1];

            taylorModels["_f" + to_string(currentNeuron)] = TaylorModel(buffer.str(), augmentedStateVars);            

            currentNeuron++;
        }

        TaylorModelVec tms;
	std::vector<bool> isIdentity(augmentedStateVars.size() - 1);

        for(int i=0; i < augmentedStateVars.size() - 1; i++ ){
            string varName = augmentedStateVars.varNames[i+1];

            if(taylorModels.find(varName) == taylorModels.end()) {

		if(!strncmp(varName.c_str(), "_f", strlen("_f"))) {
		    tms.tms.push_back(TaylorModel("0", augmentedStateVars));
		    isIdentity[i] = false;
		}
		else{
		    tms.tms.push_back(TaylorModel(varName, augmentedStateVars));
		    isIdentity[i] = true;
		}
            } else {
                tms.tms.push_back(taylorModels[varName]);
                isIdentity[i] = false;
            }
        }

	ResetMap r = ResetMap(tms, isIdentity);
        resets.push_back(ResetMap(tms, isIdentity));
        activations.push_back(activationAsEnum);
	
    }
}
