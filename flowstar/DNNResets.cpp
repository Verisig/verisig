#include "DNNResets.h"
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

	return Real(0.5) * (r2 - r1);
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

bool getTanhDerExtrema(std::vector<Real> &extrema_locations, std::vector<Real> &extrema_magnitudes, const int order){

  
        if(order == 4){
	  
	        // 5th derivative can be factored as: 8 * (1 - y^2) * (15 * y^4 - 15 * y^2 + 2)
	        // Roots are y = (15 +- sqrt(105)) / 30

	        Real sqrt105 = Real(105);
		sqrt105.sqrt_assign();

		Real sol1 = (Real(15) + sqrt105) / Real(30);
		Real sol2 = (Real(15) - sqrt105) / Real(30);

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
		Real root1 = tanhinv(tanh2);
		Real root2 = tanhinv(tanh4);
		Real root3 = tanhinv(tanh3);
		Real root4 = tanhinv(tanh1);

		extrema_locations.push_back(root1);
		extrema_locations.push_back(root2);
		extrema_locations.push_back(root3);
		extrema_locations.push_back(root4);

		extrema_magnitudes.push_back(tanhDer(4, root1));
		extrema_magnitudes.push_back(tanhDer(4, root2));
		extrema_magnitudes.push_back(tanhDer(4, root3));
		extrema_magnitudes.push_back(tanhDer(4, root4));
	}

	if(order == 5){

	        // 6th derivative can be factored as: 16 * y * (1 - y^2) * (-45 * y^4 + 60 * y^2 - 17)
	        // Roots are y = (60 +- sqrt(540)) / 90   # (multiplied top and bottom by -1)

	        Real sqrt540 = Real(540);
		sqrt540.sqrt_assign();

		Real sol1 = (Real(60) + sqrt540) / Real(90);
		Real sol2 = (Real(60) - sqrt540) / Real(90);

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
		Real root1 = tanhinv(tanh2);
		Real root2 = tanhinv(tanh4);
		Real root3 = tanhinv(tanh3);
		Real root4 = tanhinv(tanh1);

		extrema_locations.push_back(root1);
		extrema_locations.push_back(root2);
		extrema_locations.push_back(root3);
		extrema_locations.push_back(root4);

		extrema_magnitudes.push_back(tanhDer(5, root1));
		extrema_magnitudes.push_back(tanhDer(5, root2));
		extrema_magnitudes.push_back(tanhDer(5, root3));
		extrema_magnitudes.push_back(tanhDer(5, root4));		
	  
	}

        return true;
}

bool getSigDerExtrema(std::vector<Real> &extrema_locations, std::vector<Real> &extrema_magnitudes, const int order){

  
        if(order == 4){
	  
	        // 5th derivative can be factored as: (1/8) * (1 - y^2) * (15 * y^4 - 15 * y^2 + 2) # where y = tanh(x/2)
	        // Roots are y = (15 +- sqrt(105)) / 30

	        Real sqrt105 = Real(105);
		sqrt105.sqrt_assign();

		Real sol1 = (Real(15) + sqrt105) / Real(30);
		Real sol2 = (Real(15) - sqrt105) / Real(30);

		Real sig1 = Real(sol1);
		sig1.sqrt_assign();
	
		Real sig2 = Real(sol1);
		sig2.sqrt_assign();
		sig2 = Real(-1) * sig2;

		Real sig3 = Real(sol2);
		sig3.sqrt_assign();
	
		Real sig4 = Real(sol2);
		sig4.sqrt_assign();
		sig4 = Real(-1) * sig4;

		//NB: these are ordered in increasing order
		Real root1 = Real(2) * tanhinv(sig2);
		Real root2 = Real(2) * tanhinv(sig4);
		Real root3 = Real(2) * tanhinv(sig3);
		Real root4 = Real(2) * tanhinv(sig1);

		extrema_locations.push_back(root1);
		extrema_locations.push_back(root2);
		extrema_locations.push_back(root3);
		extrema_locations.push_back(root4);

		extrema_magnitudes.push_back(sigDer(4, root1));
		extrema_magnitudes.push_back(sigDer(4, root2));
		extrema_magnitudes.push_back(sigDer(4, root3));
		extrema_magnitudes.push_back(sigDer(4, root4));
	}

	if(order == 5){

	        // 6th derivative can be factored as: (1/8) * y * (1 - y^2) * (-45 * y^4 + 60 * y^2 - 17) # where y = tanh(x/2)
	        // Roots are y = (60 +- sqrt(540)) / 90   # (multiplied top and bottom by -1)

	        Real sqrt540 = Real(540);
		sqrt540.sqrt_assign();

		Real sol1 = (Real(60) + sqrt540) / Real(90);
		Real sol2 = (Real(60) - sqrt540) / Real(90);

		Real sig1 = Real(sol1);
		sig1.sqrt_assign();
	
		Real sig2 = Real(sol1);
		sig2.sqrt_assign();
		sig2 = Real(-1) * sig2;

		Real sig3 = Real(sol2);
		sig3.sqrt_assign();
	
		Real sig4 = Real(sol2);
		sig4.sqrt_assign();
		sig4 = Real(-1) * sig4;

		//NB: these are ordered in increasing order
		Real root1 = Real(2) * tanhinv(sig2);
		Real root2 = Real(2) * tanhinv(sig4);
		Real root3 = Real(2) * tanhinv(sig3);
		Real root4 = Real(2) * tanhinv(sig1);

		extrema_locations.push_back(root1);
		extrema_locations.push_back(root2);
		extrema_locations.push_back(root3);
		extrema_locations.push_back(root4);

		extrema_magnitudes.push_back(sigDer(5, root1));
		extrema_magnitudes.push_back(sigDer(5, root2));
		extrema_magnitudes.push_back(sigDer(5, root3));
		extrema_magnitudes.push_back(sigDer(5, root4));		
	  
	}

        return true;
}

void getGenericDerBounds(Real &upper, Real &lower, const bool left_segment_increasing, const std::vector<Real> extrema_locations,
			 const std::vector<Real> extrema_magnitudes, const bool tanh_act, const int order, const Interval in_bounds){

        Real cur_segment_low, cur_segment_high;

	bool cur_segment_increasing = left_segment_increasing;

	Real in_lower, in_upper;
	in_bounds.inf(in_lower);
	in_bounds.sup(in_upper);

	Real cur_low, cur_high;

	int num_extrema = extrema_locations.size();

	// there are a total of (num_extrema + 1) segments to consider (hence the <=)
	for(int i = 0; i <= num_extrema; i++){

	        if(i == 0) cur_segment_low = in_lower;
		else cur_segment_low = extrema_locations[i-1];

		if(i == num_extrema) cur_segment_high = in_upper;
		else cur_segment_high = extrema_locations[i];

		// case 0 (input interval is to the right of current segment)
		if(in_lower > cur_segment_high){
		  
		        // toggle the increasing bool
		        cur_segment_increasing = !cur_segment_increasing;
			
			continue;
		}

		// case 1 (entire input interval is contained in current segment)
		if(in_lower >= cur_segment_low && cur_segment_high >= in_upper){

		        if(cur_segment_increasing){
			        if(tanh_act){
				        upper = tanhDer(order, in_upper);
				        lower = tanhDer(order, in_lower);
				}
				else{
				        upper = sigDer(order, in_upper);
				        lower = sigDer(order, in_lower);
				}
			}
			else{
			        if(tanh_act){
				        upper = tanhDer(order, in_lower);
				        lower = tanhDer(order, in_upper);
				}
				else{
				        upper = sigDer(order, in_lower);
				        lower = sigDer(order, in_upper);
				}
			}

			return;
		}

		// case 2 (only lower bound is in current segment)
		if(in_lower >= cur_segment_low && cur_segment_high >= in_lower && in_upper > cur_segment_high){

		        if(cur_segment_increasing){
			        if(tanh_act){
				        cur_high = extrema_magnitudes[i];
				        cur_low = tanhDer(order, in_lower);
				}
				else{
				        cur_high = extrema_magnitudes[i];
				        cur_low = sigDer(order, in_lower);
				}
			}
			else{
			        if(tanh_act){
				        cur_high = tanhDer(order, in_lower);
				        cur_low = extrema_magnitudes[i];
				}
				else{
				        cur_high = sigDer(order, in_lower);
				        cur_low = extrema_magnitudes[i];
				}
			}		  
		}

		// case 3 (current segment is inside input bounds but does not contain either bound)
		if(cur_segment_low > in_lower && in_upper > cur_segment_high){

		        if(cur_segment_increasing){
			        if(extrema_magnitudes[i] > cur_high) cur_high = extrema_magnitudes[i];
			}

			else{
			        if(cur_low > extrema_magnitudes[i]) cur_low = extrema_magnitudes[i];
			}
		}

		// case 4 (current segment contains only upper bound)
		if(cur_segment_low > in_lower && in_upper >= cur_segment_low && cur_segment_high >= in_upper){
		  
		        if(cur_segment_increasing){
			        if(tanh_act){
				        if(tanhDer(order, in_upper) > cur_high) cur_high = tanhDer(order, in_upper);
				}
				else{
				        if(sigDer(order, in_upper) > cur_high) cur_high = sigDer(order, in_upper);
				}
			}

			else{
			        if(tanh_act){
				        if(cur_low > tanhDer(order, in_upper)) cur_low = tanhDer(order, in_upper);
				}
				else{
				        if(cur_low > sigDer(order, in_upper)) cur_low = sigDer(order, in_upper);
				}				
			}

			break;
		}

		// toggle the increasing bool
		cur_segment_increasing = !cur_segment_increasing;
	}

	upper = cur_high;
	lower = cur_low;
}

void getGenericDerRemBound(Interval &remBounds, const bool left_segment_increasing, const std::vector<Real> extrema_locations,
			   const std::vector<Real> extrema_magnitudes, const Interval inputBounds, const double apprPoint,
			   const bool tanh_act, const int order){
  
        Interval upper = Interval(apprPoint, inputBounds.sup());
	Interval lower = Interval(inputBounds.inf(), apprPoint);

	Real Q_u, Q_l, q_u, q_l;

	getGenericDerBounds(Q_u, q_u, left_segment_increasing, extrema_locations, extrema_magnitudes, tanh_act, order, upper);
	getGenericDerBounds(Q_l, q_l, left_segment_increasing, extrema_locations, extrema_magnitudes, tanh_act, order, lower);

	// NB: I am hardcoding the factorial since we only support orders 3 and 4 currently
	// if order == 4
	int factorial = 24;
	
	if(order == 5) factorial = 120;

	if(order > 5){
	        printf("Only orders supported currently are and 3 and 4\n");
		exit(-1);
	}
	
        Real fact = Real(factorial);
	Real maxPosDev = Real(inputBounds.sup() - apprPoint);
	Real maxNegDev = Real(inputBounds.inf()- apprPoint);
	maxPosDev.pow_assign(order);
	maxNegDev.pow_assign(order);

	Real u = (maxPosDev * Q_u) / fact;
	Real l = (maxPosDev * q_u) / fact;

	if((maxNegDev * Q_l) / fact > u) u = (maxNegDev * Q_l) / fact;
	if(l > (maxNegDev * q_l) / fact) l = (maxNegDev * q_l) / fact;

	//these checks are necessary because the remainder is always 0 at apprPoint
	if(Real(0) > u) u = Real(0);
	if(l > Real(0)) l = Real(0);
	
        remBounds = Interval(l.getValue_RNDD(), u.getValue_RNDU());
    
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
	Real ext1 = tanhinv(tanh2);
	Real ext2 = tanhinv(tanh4);
	Real ext3 = tanhinv(tanh3);
	Real ext4 = tanhinv(tanh1);

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

	// sanity checking -- this should not be reached
	else{
	        printf("Unexpected case for tanh reset. Exiting...\n");
		exit(-1);
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
	Real ext1 = tanhinv(tanh2);
	Real ext2 = tanhinv(tanh4);
	Real ext3 = tanhinv(tanh3);
	Real ext4 = tanhinv(tanh1);

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

	// sanity checking -- this should not be reached
	else{
	        printf("Unexpected case for tanh reset. Exiting...\n");
		exit(-1);
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

void dnn::tanh_reset(TaylorModel &tmReset, const Interval intC, const int varInd, const int numVars){

    double bestPoint = intC.midpoint();
    Interval bestRem = getTanh4thDerRemBound(intC, bestPoint);

    int numPointsToTry = 10;

    double range = intC.width() / 10;
    double curPoint = bestPoint;
    double high = bestPoint + range;
    double low = bestPoint - range;

    double der_granularity = intC.width() / 10000;

    Interval upRem, lowRem;
    
    for(int i = 0; i < numPointsToTry; i++){

        upRem = getTanh4thDerRemBound(intC, curPoint);
	lowRem = getTanh4thDerRemBound(intC, curPoint + der_granularity);

	double der = (upRem.width() - lowRem.width()) / der_granularity;

	if (der > 0){
	    low = curPoint;
	}
	else{
	    high = curPoint;
	}

	curPoint = (high + low) / 2;
    }

    bestPoint = curPoint;
    bestRem = upRem;

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

void polynomial_approximation(NNPolynomial &expansion, const std::vector<Real> coefficients, const Real appr_point,
			      const int var_ind, const int num_vars, const int degree){

        std::vector<int> binomial_coefficients, binomial_coefficients_prev;
	binomial_coefficients.push_back(1);
	
        Real cur_factorial = Real(1);

	//create the degree vectors
	std::map<int, std::vector<int>> degrees_map;
	for(int i = 0; i <= degree; i++){

	        std::vector<int> deg(num_vars, 0);
		deg[var_ind + 1] = i;

		degrees_map[i] = deg;
	}

	expansion = NNPolynomial();

        for(int i = 0; i <= degree; i++){

	        NNPolynomial poly_parts = NNPolynomial();

		int num_parts = i + 1;	  

		//compute binomial coefficients
		binomial_coefficients.clear();
		for(int part = 0; part < num_parts; part++){
		  
		        if(part == 0 || part == num_parts - 1){
		                binomial_coefficients.push_back(1);
				continue;
			}
		      
			binomial_coefficients.push_back(binomial_coefficients_prev[part-1] + binomial_coefficients_prev[part]);
		}

		//update factorial
		if(i > 0) cur_factorial = cur_factorial * Real(i);

		for(int part = 0; part < num_parts; part++){

		       int sign = 1;
		       if(part % 2 == 1) sign = -1;

		       int var_degree = i - part;
		       int const_degree = part;

		       Real const_coef = appr_point;
		       const_coef.pow_assign(const_degree);

		       Interval cur_coef = Interval( (coefficients[i] * const_coef *
						     Real(sign) * Real(binomial_coefficients[part]) ) / cur_factorial);

		       shared_ptr<NNMonomial> newMono(new NNMonomial(cur_coef, degrees_map[var_degree]));
		       NNPolynomial cur_poly = NNPolynomial(newMono);

		       poly_parts.add_assign(cur_poly);
		}

		//update polynomial
		expansion.add_assign(poly_parts);

		//store binomial coefficients
		binomial_coefficients_prev = binomial_coefficients;
	}
}

void dnn::act_reset(NNTaylorModel &tmReset, const Interval intC, const int varInd, const int numVars, const bool tanh_act, const int order){

    if(order > 4 || order < 3){
            printf("Verisig currently only supports 3rd and 4th order polynomial approximations of the neural network.\n");
	    exit(1);
    }

    double bestPoint = intC.midpoint();
    Interval bestRem;

    std::vector<Real> extrema_locations;
    std::vector<Real> extrema_magnitudes;

    bool left_segment_increasing;

    if(tanh_act){
            left_segment_increasing = getTanhDerExtrema(extrema_locations, extrema_magnitudes, order+1);
    }
    else{
            left_segment_increasing = getSigDerExtrema(extrema_locations, extrema_magnitudes, order+1);
    }

    //getGenericDerRemBound(bestRem, left_segment_increasing, extrema_locations,
    //			  extrema_magnitudes, intC, bestPoint, tanh_act, order+1);
    
    int numPointsToTry = 10;

    double range = intC.width() / 10;
    double curPoint = bestPoint;
    double high = bestPoint + range;
    double low = bestPoint - range;

    double der_granularity = intC.width() / 10000;

    Interval upRem, lowRem;
    
    for(int i = 0; i < numPointsToTry; i++){

	getGenericDerRemBound(upRem, left_segment_increasing, extrema_locations,
			      extrema_magnitudes, intC, curPoint, tanh_act, order+1);

	getGenericDerRemBound(lowRem, left_segment_increasing, extrema_locations,
			      extrema_magnitudes, intC, curPoint + der_granularity, tanh_act, order+1);

	double der = (upRem.width() - lowRem.width()) / der_granularity;

	if (der > 0){
	    low = curPoint;
	}
	else{
	    high = curPoint;
	}

	curPoint = (high + low) / 2;
    }

    bestPoint = curPoint;
    bestRem = upRem;

    Real apprPoint = Real(bestPoint);

    std::vector<Real> coefficients;

    if(tanh_act){

        coefficients.push_back(tanh(apprPoint));

	for(int i = 1; i < order+1; i++){
	    coefficients.push_back(tanhDer(i, apprPoint));
	}
    }
    else{
        coefficients.push_back(sigmoid(apprPoint));

	for(int i = 1; i < order+1; i++){
	    coefficients.push_back(sigDer(i, apprPoint));
	}
    }

    NNPolynomial exp;

    polynomial_approximation(exp, coefficients, apprPoint, varInd, numVars, order);

    tmReset.expansion = exp;
    
    tmReset.remainder = bestRem;


}

void dnn::convert_TM_dimension(NNTaylorModel &tmNew, const NNTaylorModel &tmOld, const int dimNew,
			       const int varInd, const std::vector<std::string> & varNames, const std::map<int, int> &indexMap){

    NNPolynomial newPoly;

    for (auto it = tmOld.expansion.monomials_map.begin(); it != tmOld.expansion.monomials_map.end(); ++it){
	
        std::vector<int> newDegs(dimNew, 0);

	shared_ptr<NNMonomial> curM = it->second;

	for(int j = 0; j < curM->getNumVars(); j++){

	    if(curM->degree() == 0) break;

	    if(indexMap.size() == 0){
	        newDegs[varInd+1] = 1;

	    }

	    // skip time
	    else if(j > 0 && indexMap.find(j-1) != indexMap.end()){		  	  
	        newDegs[indexMap.at(j-1)+1] = curM->getDegree(j);

	    }
	}

	shared_ptr<NNMonomial> newMono(new NNMonomial(curM->getCoefficient(), newDegs));

	//newPoly.monomials.push_back(newMono);
	newPoly.add_assign(newMono, varNames);
      
    }

    tmNew.expansion = NNPolynomial(newPoly);
    tmNew.remainder = Interval(tmOld.remainder);
  
}

void dnn::get_state_to_f_map(std::map<int, int> &stateToF, std::map<int, int> &fToState,
			const std::vector<std::string> &stateVarNames, const TaylorModelVec &tmv){

    std::vector<bool> init_states(stateVarNames.size());

    for(int i = 0; i < stateVarNames.size(); i++){

	if(!strncmp(stateVarNames[i].c_str(), "_f", strlen("_f"))){

	    std::list<Monomial>::const_iterator it;
	    for (it = tmv.tms[i].expansion.monomials.begin(); it != tmv.tms[i].expansion.monomials.end(); ++it){

	        for(int varInd = 1; varInd < it->getNumVars(); varInd++){ //the 1 is for time

		    if(it->getDegree(varInd) != 0){

		        init_states[varInd-1] = true;

		    }
		}
	    }
	}
    }

    int curFInd = 0;

    for(int i = 0; i < stateVarNames.size(); i++){

	if(init_states[i]){

	    fToState[curFInd] = i;
	    stateToF[i] = curFInd;

	    curFInd++;

	}
    }
}

void dnn::load_dnn(std::vector<iMatrix> &reset_weights, std::vector<iMatrix> &reset_biases, std::vector<dnn::activation> &activations,
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
	

        map<string, NNTaylorModel> taylorModels;
        Node layer = weights[layerId];
        int layerSize = layer[0].size();

	int num_vars = augmentedStateVars.size() - 1;
	
	iMatrix cur_weights(num_vars, num_vars);
	iMatrix cur_biases(layer.size(), 1);

        int currentNeuron = 1;
        for(const_iterator neuronIt=layer.begin(); neuronIt != layer.end(); ++neuronIt) {

	    //int currentWeight = 1;

	    //stringstream buffer;
            for(int i = 0; i < layerSize; i++) {

	      //buffer << (*neuronIt)[i];
	      //buffer << " * " << "_f" << currentWeight << + " + ";

	      //currentWeight++;

		Real cur_weight((double) (*neuronIt)[i].as<double>());

		cur_weights.setDataAt(currentNeuron - 1, i, cur_weight);
            }

	    Real cur_bias((double) offsets[layerId][currentNeuron-1].as<double>());

	    cur_biases.setDataAt(currentNeuron - 1, 0, cur_bias);

            //buffer << offsets[layerId][currentNeuron-1];

            //taylorModels["_f" + to_string(currentNeuron)] = NNTaylorModel(buffer.str(), augmentedStateVars);            

            currentNeuron++;
        }

        // NNTaylorModelVec tms;

        // for(int i=0; i < augmentedStateVars.size() - 1; i++ ){
        //     string varName = augmentedStateVars.varNames[i+1];

        //     if(taylorModels.find(varName) == taylorModels.end()) {

	// 	if(!strncmp(varName.c_str(), "_f", strlen("_f"))) {
	// 	    tms.tms.push_back(NNTaylorModel("0", augmentedStateVars));
	// 	}
	// 	else{
	// 	    tms.tms.push_back(NNTaylorModel(varName, augmentedStateVars));
	// 	}
        //     } else {
        //         tms.tms.push_back(taylorModels[varName]);
        //     }
        // }

        reset_weights.push_back(cur_weights);
	reset_biases.push_back(cur_biases);
        activations.push_back(activationAsEnum);
	
    }
}
