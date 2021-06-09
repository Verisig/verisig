#include "AddedResets.h"

Real sec(Real input){
        Real out = Real(input);
	
	out.cos_assign();

	Real one = Real(1);
	
        return one / out;
}

Real divide(Real input){
  
        return Real(1) / input;
}

Real arctan(double input){

        return Real(atan(input));
}

Real arctan(const Real input){

        Real result = input;
	result.atan_assign();

        return result;
}

Real tan(Real input){
        Real denom = Real(input);
	Real num = Real(input);

	num.sin_assign();
	denom.cos_assign();
  
        return num/denom;
}

Real cosine(Real input){

        Real out = Real(input);

	out.cos_assign();
  
        return out;
}

Real sine(Real input){

        Real out = Real(input);

	out.sin_assign();
  
        return out;
}

Real sqrt(Real input){
        Real output = Real(input);

	output.sqrt_assign();
  
        return output;
}

Real sec1stDer(Real input){
        return sec(input) * tan(input);
}

Real sec2ndDer(Real input){
        //return sec(input) * pow(tan(input), 2) + pow(sec(input), 3);
        return sec(input) * tan(input) * tan(input) + sec(input) * sec(input) * sec(input);
}

Real sec3rdDer(Real input){
        //return sec(input) * pow(tan(input), 3) + 5 * pow(sec(input), 3) * tan(input);
        return sec(input) * tan(input) * tan(input) * tan(input) + Real(5) * sec(input) * sec(input) * sec(input) * tan(input);
}

Real sec4thDer(Real input){
        //return sec(input) * pow(tan(input), 4) + 18 * pow(sec(input), 3) * pow(tan(input), 2) + 5 * pow(sec(input), 5);
	return sec(input) * tan(input) * tan(input) * tan(input) * tan(input) +
	  Real(18) * sec(input) * sec(input) * sec(input) * tan(input) * tan(input) +
	  Real(5) * sec(input) * sec(input) * sec(input) * sec(input) * sec(input);
}

Real sec5thDer(Real input){
	return sec(input) * tan(input) * tan(input) * tan(input) * tan(input) * tan(input) +
	  Real(58) * sec(input) * sec(input) * sec(input) * tan(input) * tan(input) * tan(input) +
	  Real(61) * sec(input) * sec(input) * sec(input) * sec(input) * sec(input) * tan(input);
}

Real div1stDer(Real input){
  
        return Real(-1) / (input * input);
}

Real div2ndDer(Real input){
  
        return Real(2) / (input * input * input);
}

Real div3rdDer(Real input){
  
        return Real(-6) / (input * input * input * input);
}

Real div4thDer(Real input){
  
        return Real(24) / (input * input * input * input * input);
}

Real div5thDer(Real input){
  
        return Real(-120) / (input * input * input * input * input * input);
}

Real arctanCoef(int order){

        if (order % 2 == 0) return Real(0);

	if((order + 1) % 4 == 2) return Real(1.0/order);

	else return Real(-1.0/order);
}

// Real arctanDer(int order, Real input){

// 	Real sum = 0;

// 	bool posSign = true;

// 	for(int i = 1; i <= order; i++){
// 	        if (i % 2 != 0){
// 		        Real nextEl = Real(input);

// 			nextEl.pow_assign(i);

// 			if(posSign) sum += nextEl/Real(i);

// 			else sum -= nextEl/Real(i);

// 			posSign = !posSign;
// 		}

// 	}

// 	return sum;
// }

Real arctanDer(const int order, const Real input){

        if(order == 1){

	        Real result = input;
		result.pow_assign(2);

		return Real(1) / (Real(1) + result);
	}

	if(order == 2){

	        Real num = input;
		num = Real(-2) * num;

		Real denom = input;
		denom.pow_assign(2);
		denom = Real(1) + denom;
		denom.pow_assign(2);

		return num / denom;
	}

	if(order == 3){

	        Real num = input;
		num.pow_assign(2);
		num = Real(6) * num + Real(-2);

		Real denom = input;
		denom.pow_assign(2);
		denom = Real(1) + denom;
		denom.pow_assign(3);

		return num / denom;		
	}

	if(order == 4){

	        Real num = input;
		Real num2 = input;
		num2.pow_assign(2);

		num = Real(-24) * num * (num2 - Real(1));

		Real denom = input;
		denom.pow_assign(2);
		denom = Real(1) + denom;
		denom.pow_assign(4);

		return num / denom;				
	}

	if(order == 5){

		Real num2 = input;
		num2.pow_assign(2);
		Real num4 = input;
		num4.pow_assign(4);
		Real num = Real(24) * (Real(5) * num4 + Real(-10) * num2 + Real(1));

		Real denom = input;
		denom.pow_assign(2);
		denom = Real(1) + denom;
		denom.pow_assign(5);

		return num / denom;
	}
  
}

Real tan1stDer(Real input){
        return Real(1) + tan(input) * tan(input);
}

Real tan2ndDer(Real input){
        return Real(2) * tan(input) + Real(2) * tan(input) * tan(input) * tan(input);
}

Real tan3rdDer(Real input){
        return Real(2) + Real(8) * tan(input) * tan(input) + Real(6) * tan(input) * tan(input) * tan(input) * tan(input);
}

Real tan4thDer(Real input){
        return Real(16) * tan(input) + Real(40) * tan(input) * tan(input) * tan(input) +
	  Real(24) * tan(input) * tan(input) * tan(input) * tan(input) * tan(input);
}

Real cos1stDer(Real input){

        Real out = Real(input);

	out.sin_assign();
  
        return Real(-1.0) * out;
}

Real cos2ndDer(Real input){

        Real out = Real(input);

	out.cos_assign();
  
        return Real(-1.0) * out;
}

Real cos3rdDer(Real input){

        Real out = Real(input);

	out.sin_assign();
  
        return out;
}

Real cos4thDer(Real input){

        Real out = Real(input);

	out.cos_assign();
  
        return out;
}

Real cosDer(Real input, int order){

        Real out = Real(input);

        if(order % 4 == 1){
		out.sin_assign();
		out = Real(-1.0) * out;
	}

        if(order % 4 == 2){
		out.cos_assign();
		out = Real(-1.0) * out;
	}

        if(order % 4 == 3){
		out.sin_assign();
	}

	if(order % 4 == 0){
		out.cos_assign();
	}

	return out;
}

Real sinDer(Real input, int order){

        Real out = Real(input);

        if(order % 4 == 0){
		out.sin_assign();
	}

	if(order % 4 == 1){
		out.cos_assign();
	}

        if(order % 4 == 2){
		out.sin_assign();
		out = Real(-1.0) * out;
	}

        if(order % 4 == 3){
		out.cos_assign();
		out = Real(-1.0) * out;
	}	

	return out;
}

Real sin1stDer(Real input){

        Real out = Real(input);

	out.cos_assign();
  
        return out;
}

Real sin2ndDer(Real input){

        Real out = Real(input);

	out.sin_assign();
  
        return Real(-1.0) * out;
}

Real sin3rdDer(Real input){

        Real out = Real(input);

	out.cos_assign();
  
        return Real(-1.0) * out;
}

Real sin4thDer(Real input){

        Real out = Real(input);

	out.sin_assign();
  
        return out;
}

Real sqrt1stDer(Real input){

        Real inSqrt = Real(input);

	inSqrt.sqrt_assign();
  
        return Real(0.5) / inSqrt;
}

Real sqrt2ndDer(Real input){
  
        Real inSqrt = Real(input);

	inSqrt.sqrt_assign();
  
        return Real(-0.25) / (inSqrt * Real(input));
}

Real sqrt3rdDer(Real input){
        Real inSqrt = Real(input);

	inSqrt.sqrt_assign();

	Real out = Real( 3.0 / 8.0 ) / (inSqrt * Real(input) * Real(input));
  
        return out;
}

Real sqrt4thDer(Real input){
        Real inSqrt = Real(input);

	inSqrt.sqrt_assign();
  
        return Real(- 15.0 / 16.0 ) / (inSqrt * Real(input) * Real(input) * Real(input));
}

Real getSecDerBound(Interval intC, int order, bool upper){

        Real derBound;
	if (intC.inf() > -M_PI/2 && intC.sup() < M_PI/2){

	        if (order == 5){

		        if(upper) derBound = sec5thDer(intC.sup());

			else derBound = sec5thDer(intC.inf());
			
		}
	  
	        if (order == 4){
		        if(upper){
			        derBound = sec4thDer(intC.sup());

				if(sec4thDer(intC.inf()) > derBound) derBound = sec4thDer(intC.inf());
			}

			else{
			        derBound = sec4thDer(intC.inf());

				if(derBound > sec4thDer(intC.sup())) derBound = intC.sup();

				if(Real(0) >= Real(intC.inf()) && Real(intC.sup()) >= Real(0)) derBound = Real(5);
			}
		}

	        if (order == 3){
		        if(upper) derBound = sec3rdDer(intC.sup());

			else derBound = sec3rdDer(intC.inf());
		}	
	}
	else if (intC.inf() > M_PI/2 && intC.sup() <= M_PI){
	        if (order == 5){
		        if(upper) derBound = sec5thDer(intC.inf());

			else derBound = sec5thDer(intC.sup());
		}
	  
	        if (order == 4){
		        if(upper) derBound = sec4thDer(intC.sup());

			else derBound = sec4thDer(intC.inf());
		}

		if (order == 3){
		        if(upper) derBound = sec3rdDer(intC.inf());

			else derBound = sec3rdDer(intC.sup());
		}
	}
	else if (intC.inf() >= -M_PI && intC.sup() < -M_PI/2){
	        if (order == 5){
		        if(upper) derBound = sec5thDer(intC.inf());

			else derBound = sec5thDer(intC.sup());
		}
	  
	        if (order == 4){
		        if(upper) derBound = sec4thDer(intC.inf());

			else derBound = sec4thDer(intC.sup());
		}

		if (order == 3){
		        if(upper) derBound = sec3rdDer(intC.inf());

			else derBound = sec3rdDer(intC.sup());
		}
	}
	else{
	        printf("Uncertainty too large. Please try decreasing the initial set size.\n");
		exit(-1);
	}

	return derBound;
}

// Interval getSecDerRemBound(const Interval inputBounds, double apprPoint, int order){

//         Interval upper = Interval(apprPoint, inputBounds.sup());
// 	Interval lower = Interval(inputBounds.inf(), apprPoint);
  
//         Real Q_u = getSecDerBound(upper, order, true);
// 	Real Q_l = getSecDerBound(upper, order, false);
// 	Real q_u = getSecDerBound(lower, order, true);
// 	Real q_l = getSecDerBound(lower, order, false);

//         Real fact = Real(24);
// 	Real maxPosDev = Real(inputBounds.sup() - apprPoint);
// 	Real maxNegDev = Real(inputBounds.inf()- apprPoint);

// 	if(order == 4){
// 	        maxPosDev.pow_assign(4);
// 		maxNegDev.pow_assign(4);
// 	}

// 	if(order == 3){
// 	        fact = Real(6);
// 		maxPosDev.pow_assign(3);
// 		maxNegDev.pow_assign(3);
// 	}

// 	Real u = (maxPosDev * Q_u) / fact;
// 	Real l = (maxPosDev * q_u) / fact;

// 	if((maxNegDev * Q_l) / fact > u) u = (maxNegDev * Q_l) / fact;
// 	if(l > (maxNegDev * q_l) / fact) l = (maxNegDev * q_l) / fact;

// 	//these checks are necessary because the remainder is always 0 at apprPoint
// 	if(Real(0) > u) u = Real(0);
// 	if(l > Real(0)) l = Real(0);
	
//         return Interval(l.getValue_RNDD(), u.getValue_RNDU());
    
// }

Real getSecRemUpperBound(Interval intC, int order){

        Real derBound;
	if (intC.inf() > -M_PI/2 && intC.sup() < M_PI/2){

	        if (order == 5){
		        derBound = sec5thDer(intC.sup());
		}
	  
	        if (order == 4){
		        derBound = sec4thDer(intC.sup());
			
			if (sec4thDer(intC.inf()) > derBound)
			        derBound = sec4thDer(intC.inf());
		}

	        if (order == 3){
		        derBound = sec3rdDer(intC.sup());
		}	
	}
	else if (intC.inf() > M_PI/2 && intC.sup() <= M_PI){
	        if (order == 5)
		        derBound = sec5thDer(intC.inf());
	  
	        if (order == 4)
		        derBound = sec4thDer(intC.sup());

		if (order == 3)
		        derBound = sec3rdDer(intC.inf());
	}
	else if (intC.inf() >= -M_PI && intC.sup() < -M_PI/2){
	        if (order == 5)
		        derBound = sec5thDer(intC.inf());
	  
	        if (order == 4)
		        derBound = sec4thDer(intC.inf());

		if (order == 3)		  
		        derBound = sec3rdDer(intC.inf());
	}
	else{
	        printf("Uncertainty too large. Please try decreasing the initial set size.\n");
		exit(-1);
	}

	return derBound;
}

Real getSecRemLowerBound(Interval intC, int order){

        Real derBound;
	if (intC.inf() > -M_PI/2 && intC.sup() < M_PI/2){

	        if (order == 5){
		        derBound = sec5thDer(intC.inf());
		}
	  
	        if (order == 4){
		        derBound = sec4thDer(intC.sup());
			
			if (derBound > sec4thDer(intC.inf()))
			        derBound = sec4thDer(intC.inf());
		}

	        if (order == 3){
		        derBound = sec3rdDer(intC.inf());
		}	
	}
	else if (intC.inf() > M_PI/2 && intC.sup() <= M_PI){
	        if (order == 5)
		        derBound = sec5thDer(intC.sup());
	  
	        if (order == 4)
		        derBound = sec4thDer(intC.inf());

		if (order == 3)
		        derBound = sec3rdDer(intC.sup());
	}
	else if (intC.inf() >= -M_PI && intC.sup() < -M_PI/2){
	        if (order == 5)
		        derBound = sec5thDer(intC.sup());
	  
	        if (order == 4)
		        derBound = sec4thDer(intC.sup());

		if (order == 3)		  
		        derBound = sec3rdDer(intC.sup());
	}
	else{
	        printf("Uncertainty too large. Please try decreasing the initial set size.\n");
		exit(-1);
	}

	return derBound;
}

Interval getSecDerRemBound(const Interval inputBounds, const double apprPoint, const int order){
  
        Interval upper = Interval(apprPoint, inputBounds.sup());
	Interval lower = Interval(inputBounds.inf(), apprPoint);
  
        Real Q_u = getSecRemUpperBound(upper, order);
	Real Q_l = getSecRemUpperBound(lower, order);
	Real q_u = getSecRemLowerBound(upper, order);
	Real q_l = getSecRemLowerBound(lower, order);

	// default for order 4
        Real fact = Real(24);
	
	if(order == 5)
	  fact = Real(120);

	if(order == 3)
	  fact = Real(6);
	
	Real maxPosDev = Real(inputBounds.sup() - apprPoint);
	Real maxNegDev = Real(inputBounds.inf() - apprPoint);
	maxPosDev.pow_assign(order);
	maxNegDev.pow_assign(order);

	Real u = (maxPosDev * Q_u) / fact;
	Real l = (maxPosDev * q_u) / fact;

	if((maxNegDev * Q_l) / fact > u) u = (maxNegDev * Q_l) / fact;
	if(l > (maxNegDev * q_l) / fact) l = (maxNegDev * q_l) / fact;

	//these checks are necessary because the remainder is always 0 at apprPoint
	if(Real(0) > u) u = Real(0);
	if(l > Real(0)) l = Real(0);
	
        return Interval(l.getValue_RNDD(), u.getValue_RNDU());
    
}

Real getDivRemUpperBound(const Interval intC, const int order){

        Real bound = Real(0);

        if(order == 3){
	        //3rd derivative is negative and decreasing for negative numbers
	        if (intC.sup() < 0){
		        bound = div3rdDer(Real(intC.inf()));
		}
		//3rd derivative is negative and increasing for positive numbers
		else if (intC.inf() > 0){
		        bound = div3rdDer(Real(intC.sup()));
		}
	}

	else if(order == 4){
	        //4th derivative is negative and decreasing for negative numbers
	        if (intC.sup() < 0){
		        bound = div3rdDer(Real(intC.inf()));
		}
		//4th derivative is positive and decreasing for positive numbers
		else if (intC.inf() > 0){
		        bound = div3rdDer(Real(intC.inf()));
		}
	}

	else if(order == 5){
	        //5th derivative is negative and decreasing for negative numbers
	        if (intC.sup() < 0){
		        bound = div3rdDer(Real(intC.inf()));
		}
		//5th derivative is negative and increasing for positive numbers
		else if (intC.inf() > 0){
		        double maxVal = fabs(div5thDer(Real(intC.inf())).getValue_RNDD());

		        bound = div3rdDer(Real(intC.sup()));
		}
	}

	return bound;
}

Real getDivRemLowerBound(const Interval intC, const int order){

        Real bound = Real(0);

        if(order == 3){
	        //3rd derivative is negative and decreasing for negative numbers
	        if (intC.sup() < 0){
		        bound = div3rdDer(Real(intC.sup()));
		}
		//3rd derivative is negative and increasing for positive numbers
		else if (intC.inf() > 0){
		        bound = div3rdDer(Real(intC.inf()));
		}
	}

	else if(order == 4){
	        //4th derivative is negative and decreasing for negative numbers
	        if (intC.sup() < 0){
		        bound = div3rdDer(Real(intC.sup()));
		}
		//4th derivative is positive and decreasing for positive numbers
		else if (intC.inf() > 0){
		        bound = div3rdDer(Real(intC.sup()));
		}
	}

	else if(order == 5){
	        //5th derivative is negative and decreasing for negative numbers
	        if (intC.sup() < 0){
		        bound = div3rdDer(Real(intC.sup()));
		}
		//5th derivative is negative and increasing for positive numbers
		else if (intC.inf() > 0){
		        bound = div3rdDer(Real(intC.inf()));
		}
	}

	return bound;
}

Interval getDivDerRemBound(const Interval inputBounds, const double apprPoint, const int order){
  
        Interval upper = Interval(apprPoint, inputBounds.sup());
	Interval lower = Interval(inputBounds.inf(), apprPoint);
  
        Real Q_u = getDivRemUpperBound(upper, order);
	Real Q_l = getDivRemUpperBound(lower, order);
	Real q_u = getDivRemLowerBound(upper, order);
	Real q_l = getDivRemLowerBound(lower, order);

        Real fact = Real(24);
	if(order == 5)
	  fact = Real(120);
	
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
	
        return Interval(l.getValue_RNDD(), u.getValue_RNDU());
    
}

Real getDivDerBound(int order, Interval intC){

        Real bound = Real(0);

        if(order == 3){
	        //3rd derivative is negative and decreasing for negative numbers
	        if (intC.sup() < 0){
		        double maxVal = fabs(div3rdDer(Real(intC.sup())).getValue_RNDD());

			bound = Real(maxVal);
		}
		//3rd derivative is negative and increasing for positive numbers
		else if (intC.inf() > 0){
		        double maxVal = fabs(div3rdDer(Real(intC.inf())).getValue_RNDD());

			bound = Real(maxVal);
		}
	}

	else if(order == 4){
	        //4th derivative is negative and decreasing for negative numbers
	        if (intC.sup() < 0){
		        double maxVal = fabs(div4thDer(Real(intC.sup())).getValue_RNDD());

			bound = Real(maxVal);
		}
		//4th derivative is positive and decreasing for positive numbers
		else if (intC.inf() > 0){
		        double maxVal = fabs(div4thDer(Real(intC.inf())).getValue_RNDD());

			bound = Real(maxVal);
		}
	}

	else if(order == 5){
	        //5th derivative is negative and decreasing for negative numbers
	        if (intC.sup() < 0){
		        double maxVal = fabs(div5thDer(Real(intC.sup())).getValue_RNDD());

			bound = Real(maxVal);
		}
		//5th derivative is negative and increasing for positive numbers
		else if (intC.inf() > 0){
		        double maxVal = fabs(div5thDer(Real(intC.inf())).getValue_RNDD());

			bound = Real(maxVal);
		}
	}

	return bound;
}

void getGenericDerBounds(Real &upper, Real &lower, const bool left_segment_increasing, const std::vector<Real> extrema_locations,
			 const std::vector<Real> extrema_magnitudes, const int reset_type, const int order, const Interval in_bounds){

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
			        if(reset_type == ARC){
				        upper = arctanDer(order, in_upper);
				        lower = arctanDer(order, in_lower);
				}
			}
			else{
			        if(reset_type == ARC){
				        upper = arctanDer(order, in_lower);
				        lower = arctanDer(order, in_upper);
				}
			}

			return;
		}

		// case 2 (only lower bound is in current segment)
		if(in_lower >= cur_segment_low && cur_segment_high >= in_lower && in_upper > cur_segment_high){

		        if(cur_segment_increasing){
			        if(reset_type == ARC){
				        cur_high = extrema_magnitudes[i];
				        cur_low = arctanDer(order, in_lower);
				}
			}
			else{
			        if(reset_type == ARC){
				        cur_high = arctanDer(order, in_lower);
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
			        if(reset_type == ARC){
				        if(arctanDer(order, in_upper) > cur_high) cur_high = arctanDer(order, in_upper);
				}
			}

			else{
			        if(reset_type == ARC){
				        if(cur_low > arctanDer(order, in_upper)) cur_low = arctanDer(order, in_upper);
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
			   const int reset_type, const int order){
  
        Interval upper = Interval(apprPoint, inputBounds.sup());
	Interval lower = Interval(inputBounds.inf(), apprPoint);

	Real Q_u, Q_l, q_u, q_l;

	getGenericDerBounds(Q_u, q_u, left_segment_increasing, extrema_locations, extrema_magnitudes, reset_type, order, upper);
	getGenericDerBounds(Q_l, q_l, left_segment_increasing, extrema_locations, extrema_magnitudes, reset_type, order, lower);

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

bool getArctanDerExtrema(std::vector<Real> &extrema_locations, std::vector<Real> &extrema_magnitudes, const int order){

  
        if(order == 4){
	  
	        // 5th derivative numerator can be factored as: 24 * (5 * x^4 - 10 * x^2 + 1)
	        // Roots are x^2 = (10 +- sqrt(20)) / 10

	        Real sqrt20 = Real(20);
		sqrt20.sqrt_assign();

		Real sol1 = (Real(10) + sqrt20) / Real(20);
		Real sol2 = (Real(10) - sqrt20) / Real(20);

		Real root1 = Real(sol1);
		root1.sqrt_assign();
	
		Real root2 = Real(sol1);
		root2.sqrt_assign();
		root2 = Real(-1) * root2;

		Real root3 = Real(sol2);
		root3.sqrt_assign();
	
		Real root4 = Real(sol2);
		root4.sqrt_assign();
		root4 = Real(-1) * root4;

		//NB: these are ordered in increasing order
		extrema_locations.push_back(root2);
		extrema_locations.push_back(root4);
		extrema_locations.push_back(root3);
		extrema_locations.push_back(root1);

		extrema_magnitudes.push_back(arctanDer(4, root2));
		extrema_magnitudes.push_back(arctanDer(4, root4));
		extrema_magnitudes.push_back(arctanDer(4, root3));
		extrema_magnitudes.push_back(arctanDer(4, root1));
	}

	if(order == 5){
	  
	        // 6th derivative numerator can be factored as: -240 * x * (3 * x^4 - 10 * x^2 + 3)
	        // Roots are x = 0, x^2 = (10 +- sqrt(64)) / 6

	        Real sqrt64 = Real(64);
		sqrt64.sqrt_assign();

		Real sol1 = (Real(10) + sqrt64) / Real(6);
		Real sol2 = (Real(10) - sqrt64) / Real(6);

		Real root1 = Real(sol1);
		root1.sqrt_assign();
	
		Real root2 = Real(sol1);
		root2.sqrt_assign();
		root2 = Real(-1) * root2;

		Real root3 = Real(sol2);
		root3.sqrt_assign();
	
		Real root4 = Real(sol2);
		root4.sqrt_assign();
		root4 = Real(-1) * root4;

		Real root5 = Real(0);

		//NB: these are ordered in increasing order
		extrema_locations.push_back(root2);
		extrema_locations.push_back(root4);
		extrema_locations.push_back(root5);
		extrema_locations.push_back(root3);
		extrema_locations.push_back(root1);

		extrema_magnitudes.push_back(arctanDer(4, root2));
		extrema_magnitudes.push_back(arctanDer(4, root4));
		extrema_magnitudes.push_back(arctanDer(4, root5));
		extrema_magnitudes.push_back(arctanDer(4, root3));
		extrema_magnitudes.push_back(arctanDer(4, root1));
	}	

        return true;
}

//NB: this assumes intC \subset [-1, 1]
Real getArcTanDerBound(double dev, int order){

	return Real( (pow(dev, 2 * order + 1)) / (2 * order + 1));
}

Real getTanDerBound(int order, Interval intC){

	if(order == 4){
	        double low = fabs(tan4thDer(Real(intC.inf())).getValue_RNDD());
		double high = fabs(tan4thDer(Real(intC.sup())).getValue_RNDU());

		if (low > high) return Real(low);
		else return Real(high);
		  
	}

	else{ //if order == 3
	        double low = fabs(tan3rdDer(Real(intC.inf())).getValue_RNDD());
		double high = fabs(tan3rdDer(Real(intC.sup())).getValue_RNDU());

		if (low > high) return Real(low);
		else return Real(high);
		  
	}
}

Real getCosDerBound(int order, Interval intC){

        Real bound = Real(1);

        if(order == 3){ //derivative is sin(x)

	        if(intC.sup() - intC.inf() < M_PI){
		        int pio2Mult = ceil((intC.inf() + M_PI/2) / M_PI);

			double rem = pio2Mult * M_PI - intC.inf() - M_PI/2;

			//printf("rem: %f\n", rem);

			if (intC.sup() - intC.inf() < rem){
			        double maxVal = fabs(cos3rdDer(Real(intC.inf())).getValue_RNDD());

				if (fabs(cos3rdDer(Real(intC.sup())).getValue_RNDD()) > maxVal){

				        maxVal = fabs(cos3rdDer(Real(intC.sup())).getValue_RNDD());
					
				}

				bound = Real(maxVal);
			}
		}

	}

	else if(order == 4){ //derivative is cos(x)
	        if(intC.sup() - intC.inf() < M_PI){
		        int pio2Mult = ceil(intC.inf() / M_PI);

			double rem = pio2Mult * M_PI - intC.inf();

			if (intC.sup() - intC.inf() < rem){
			        double maxVal = fabs(cos4thDer(Real(intC.inf())).getValue_RNDD());

				if (fabs(cos4thDer(Real(intC.sup())).getValue_RNDD()) > maxVal){

				        maxVal = fabs(cos4thDer(Real(intC.sup())).getValue_RNDD());
					
				}

				bound = Real(maxVal);
			}
		}
	}

	return bound;

}

Real getSinRemUpperBound(Interval intC, int order){

        Real bound = Real(1);

	if(intC.sup() - intC.inf() < 2 * M_PI){
	  
	        int next_peak_index = 0;
		double dist_to_peak = 0;

		if(order % 4 == 0){ //derivative is sin(x)
		        next_peak_index = floor((intC.inf() + (3 * M_PI) / 2) / (2 * M_PI));
			dist_to_peak = (M_PI) / 2 + next_peak_index * 2 * M_PI - intC.inf();
		}		

		else if(order % 4 == 1){ //derivative is cos(x)
		        next_peak_index = ceil(intC.inf() / (2 * M_PI));
			dist_to_peak = next_peak_index * 2 * M_PI - intC.inf();
		}

		else if(order % 4 == 2){ //derivative is -sin(x)
		        next_peak_index = floor((intC.inf() + (M_PI) / 2) / (2 * M_PI));
			dist_to_peak = (3 * M_PI) / 2 + next_peak_index * 2 * M_PI - intC.inf();
		}

		else if(order % 4 == 3){ //derivative is -cos(x)
		        next_peak_index = floor(intC.inf() + M_PI / (2 * M_PI));
			dist_to_peak = M_PI + next_peak_index * 2 * M_PI - intC.inf();
		}		

		if (intC.sup() - intC.inf() < dist_to_peak){

			Real maxVal = cosDer(Real(intC.inf()), order);

			if (cosDer(Real(intC.sup()), order) > maxVal){

				maxVal = cosDer(Real(intC.sup()), order);

			}

			bound = maxVal;
		}
	}

	return bound;

}

Real getSinRemLowerBound(Interval intC, int order){

        Real bound = Real(-1);

	if(intC.sup() - intC.inf() < 2 * M_PI){
	  
	        int next_valley_index = 0;
		double dist_to_valley = 0;

		if(order % 4 == 0){ //derivative is sin(x)
		        next_valley_index = floor((intC.inf() + (M_PI) / 2) / (2 * M_PI));
			dist_to_valley = (3 * M_PI) / 2 + next_valley_index * 2 * M_PI - intC.inf();
		}		

		else if(order % 4 == 1){ //derivative is cos(x)
		        next_valley_index = floor(intC.inf() + M_PI / (2 * M_PI));
			dist_to_valley = M_PI + next_valley_index * 2 * M_PI - intC.inf();
		}

		else if(order % 4 == 2){ //derivative is -sin(x)
		        next_valley_index = floor((intC.inf() + (3 * M_PI) / 2) / (2 * M_PI));
			dist_to_valley = (M_PI) / 2 + next_valley_index * 2 * M_PI - intC.inf();
		}

		else if(order % 4 == 3){ //derivative is -cos(x)
		        next_valley_index = ceil((intC.inf()) / (2 * M_PI));
			dist_to_valley = next_valley_index * 2 * M_PI - intC.inf();
		}		

		if (intC.sup() - intC.inf() < dist_to_valley){

			Real minVal = cosDer(Real(intC.inf()), order);

			if (minVal > cosDer(Real(intC.sup()), order)){

				minVal = cosDer(Real(intC.sup()), order);

			}

			bound = minVal;
		}
	}

	return bound;

}

Real getCosRemUpperBound(Interval intC, int order){

        Real bound = Real(1);

	if(intC.sup() - intC.inf() < 2 * M_PI){
	  
	        int next_peak_index = 0;
		double dist_to_peak = 0;

		if(order % 4 == 0){ //derivative is cos(x)
		        next_peak_index = ceil(intC.inf() / (2 * M_PI));
			dist_to_peak = next_peak_index * 2 * M_PI - intC.inf();
		}

		else if(order % 4 == 1){ //derivative is -sin(x)
		        next_peak_index = floor((intC.inf() + (M_PI) / 2) / (2 * M_PI));
			dist_to_peak = (3 * M_PI) / 2 + next_peak_index * 2 * M_PI - intC.inf();
		}

		if(order % 4 == 2){ //derivative is cos(x)
		        next_peak_index = floor(intC.inf() + M_PI / (2 * M_PI));
			dist_to_peak = M_PI + next_peak_index * 2 * M_PI - intC.inf();
		}		

		else if(order % 4 == 3){ //derivative is sin(x)
		        next_peak_index = floor((intC.inf() + (3 * M_PI) / 2) / (2 * M_PI));
			dist_to_peak = (M_PI) / 2 + next_peak_index * 2 * M_PI - intC.inf();
		}

		if (intC.sup() - intC.inf() < dist_to_peak){

			Real maxVal = cosDer(Real(intC.inf()), order);

			if (cosDer(Real(intC.sup()), order) > maxVal){

				maxVal = cosDer(Real(intC.sup()), order);

			}

			bound = maxVal;
		}
	}

	return bound;

}

Real getCosRemLowerBound(Interval intC, int order){

        Real bound = Real(-1);

	if(intC.sup() - intC.inf() < 2 * M_PI){
	  
	        int next_valley_index = 0;
		double dist_to_valley = 0;

		if(order % 4 == 0){ //derivative is cos(x)
		        next_valley_index = floor(intC.inf() + M_PI / (2 * M_PI));
			dist_to_valley = M_PI + next_valley_index * 2 * M_PI - intC.inf();
		}

		else if(order % 4 == 1){ //derivative is -sin(x)
		        next_valley_index = floor((intC.inf() + (3 * M_PI) / 2) / (2 * M_PI));
			dist_to_valley = (M_PI) / 2 + next_valley_index * 2 * M_PI - intC.inf();
		}

		else if(order % 4 == 2){ //derivative is -cos(x)
		        next_valley_index = ceil((intC.inf()) / (2 * M_PI));
			dist_to_valley = next_valley_index * 2 * M_PI - intC.inf();
		}		

		else if(order % 4 == 3){ //derivative is sin(x)
		        next_valley_index = floor((intC.inf() + (M_PI) / 2) / (2 * M_PI));
			dist_to_valley = (3 * M_PI) / 2 + next_valley_index * 2 * M_PI - intC.inf();
		}

		if (intC.sup() - intC.inf() < dist_to_valley){

			Real minVal = cosDer(Real(intC.inf()), order);

			if (minVal > cosDer(Real(intC.sup()), order)){

				minVal = cosDer(Real(intC.sup()), order);

			}

			bound = minVal;
		}
	}

	return bound;

}

Interval getCosDerRemBound(const Interval inputBounds, const double apprPoint, const int order){
  
        Interval upper = Interval(apprPoint, inputBounds.sup());
	Interval lower = Interval(inputBounds.inf(), apprPoint);
  
        Real Q_u = getCosRemUpperBound(upper, order);
	Real Q_l = getCosRemUpperBound(lower, order);
	Real q_u = getCosRemLowerBound(upper, order);
	Real q_l = getCosRemLowerBound(lower, order);

        Real fact = Real(24);
	if(order == 5)
	  fact = Real(120);
	if(order == 3)
	  fact = Real(6);
	
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
	
        return Interval(l.getValue_RNDD(), u.getValue_RNDU());
    
}

Interval getSinDerRemBound(const Interval inputBounds, const double apprPoint, const int order){
  
        Interval upper = Interval(apprPoint, inputBounds.sup());
	Interval lower = Interval(inputBounds.inf(), apprPoint);
  
        Real Q_u = getSinRemUpperBound(upper, order);
	Real Q_l = getSinRemUpperBound(lower, order);
	Real q_u = getSinRemLowerBound(upper, order);
	Real q_l = getSinRemLowerBound(lower, order);

        Real fact = Real(24);
	if(order == 5)
	  fact = Real(120);
	if(order == 3)
	  fact = Real(6);
	
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
	
        return Interval(l.getValue_RNDD(), u.getValue_RNDU());
}

Real getSinDerBound(int order, Interval intC){

        Real bound = Real(1);

        if(order == 4){ //derivative is sin(x)

	        if(intC.sup() - intC.inf() < M_PI){
		        int pio2Mult = ceil((intC.inf() + M_PI/2) / M_PI);

			double rem = pio2Mult * M_PI - intC.inf() - M_PI/2;

			//printf("rem: %f\n", rem);

			if (intC.sup() - intC.inf() < rem){
			        double maxVal = fabs(cos3rdDer(Real(intC.inf())).getValue_RNDD());

				if (fabs(cos3rdDer(Real(intC.sup())).getValue_RNDD()) > maxVal){

				        maxVal = fabs(cos3rdDer(Real(intC.sup())).getValue_RNDD());
					
				}

				bound = Real(maxVal);
			}
		}

	}

	else if(order == 3){ //derivative is -cos(x)
	        if(intC.sup() - intC.inf() < M_PI){
		        int pio2Mult = ceil(intC.inf() / M_PI);

			double rem = pio2Mult * M_PI - intC.inf();

			if (intC.sup() - intC.inf() < rem){
			        double maxVal = fabs(cos4thDer(Real(intC.inf())).getValue_RNDD());

				if (fabs(cos4thDer(Real(intC.sup())).getValue_RNDD()) > maxVal){

				        maxVal = fabs(cos4thDer(Real(intC.sup())).getValue_RNDD());
					
				}

				bound = Real(maxVal);
			}
		}
	}

	return bound;

}

Real getSqrtDerBound(int order, Interval intC){

        Real bound = Real(0);

        if(order == 3){ //3rd derivative is positive and decreasing so the max is the absolute value of intC.inf()

	        double maxVal = fabs(sqrt3rdDer(Real(intC.inf())).getValue_RNDD());


		bound = Real(maxVal);
	}

	else if(order == 4){ //4th derivative is (negative and) increasing so the max is the absolute value of intC.inf()
	        double maxVal = fabs(sqrt4thDer(Real(intC.inf())).getValue_RNDD());

		bound = Real(maxVal);

	}

	return bound;

}

void sqrt_reset(TaylorModel &tmReset, const Interval intC, const int varStoreInd, const int varInputInd, const int numVars){

    Polynomial exp;
    Interval rem;

					
    Real midPoint = Real(intC.midpoint());
    
    Real apprPoint = sqrt(intC.midpoint());
    
    //NB: This assumes a 2nd order TS approximation
    Real coef1 = sqrt1stDer(midPoint);
    Real coef2 = sqrt2ndDer(midPoint)/2;
    Real coef3 = sqrt3rdDer(midPoint)/6;
    
    Real derBound = getSqrtDerBound(3, intC);

    Real maxDev = Real(intC.sup()) - midPoint;
    if (midPoint - Real(intC.inf()) > maxDev){
        maxDev = midPoint - Real(intC.inf());
    }

    Real fact = 6;
    maxDev.pow_assign(3);
    
    Real remainder = (derBound * maxDev) / fact;
    
    Interval apprInt = Interval(apprPoint);
    
    Interval deg1Int = Interval(coef1);
    Interval deg2Int = Interval(coef2);
    Interval deg3Int = Interval(coef3);

    std::vector<int> deg1(numVars, 0);
    deg1[varInputInd + 1] = 1;
    std::vector<int> deg2(numVars, 0);
    deg2[varInputInd + 1] = 2;
    std::vector<int> deg3(numVars, 0);
    deg3[varInputInd + 1] = 3;

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
	derBound = getSqrtDerBound(4, intC);
						
	remainder = (derBound * maxDev) / fact;
	remainder.to_sym_int(rem);

	Polynomial deg3Poly = Polynomial(Monomial(Interval(Real(-1) * coef3 * midPoint * midPoint * midPoint), numVars)) +
	  Polynomial(Monomial(Interval(Real(3) * coef3 * midPoint * midPoint), deg1)) -
	  Polynomial(Monomial(Interval(Real(3) * coef3 * midPoint), deg2)) +
	  Polynomial(Monomial(deg3Int, deg3));
						
	exp += deg3Poly;
    }					

    tmReset.expansion = exp;
    tmReset.remainder = rem;    
}

void sin_reset(TaylorModel &tmReset, const Interval intC, const int varStoreInd, const int varInputInd, const int numVars){

    Polynomial exp;
    Interval rem;

    Real midPoint = Real(intC.midpoint());

    Real apprPoint = sine(intC.midpoint());

    //NB: This assumes a 2nd order TS approximation
    Real coef1 = sinDer(midPoint, 1);
    Real coef2 = sinDer(midPoint, 2)/2;
    Real coef3 = sinDer(midPoint, 3)/6;

    rem = getSinDerRemBound(intC, intC.midpoint(), 4);
    
    Interval apprInt = Interval(apprPoint);
    
    Interval deg1Int = Interval(coef1);
    Interval deg2Int = Interval(coef2);
    Interval deg3Int = Interval(coef3);
    
    std::vector<int> deg1(numVars, 0);
    deg1[varInputInd + 1] = 1;
    std::vector<int> deg2(numVars, 0);
    deg2[varInputInd + 1] = 2;
    std::vector<int> deg3(numVars, 0);
    deg3[varInputInd + 1] = 3;
					
    Polynomial deg0Poly = Polynomial(Monomial(apprInt, numVars));

    Polynomial deg1Poly = Polynomial(Monomial(Interval(Real(-1) * coef1 * midPoint), numVars)) +
      Polynomial(Monomial(deg1Int, deg1));

    Polynomial deg2Poly = Polynomial(Monomial(Interval(coef2 * midPoint * midPoint), numVars)) -
      Polynomial(Monomial(Interval(Real(2) * coef2 * midPoint), deg1)) +
      Polynomial(Monomial(deg2Int, deg2));

    Polynomial deg3Poly = Polynomial(Monomial(Interval(Real(-1) * coef3 * midPoint * midPoint * midPoint), numVars)) +
      Polynomial(Monomial(Interval(Real(3) * coef3 * midPoint * midPoint), deg1)) -
      Polynomial(Monomial(Interval(Real(3) * coef3 * midPoint), deg2)) +
      Polynomial(Monomial(deg3Int, deg3));    
					
    exp = deg0Poly + deg1Poly + deg2Poly + deg3Poly;
		    
    tmReset.expansion = exp;
    tmReset.remainder = rem;    
}

void cos_reset(TaylorModel &tmReset, const Interval intC, const int varStoreInd, const int varInputInd, const int numVars){

    Polynomial exp;
    Interval rem;

    Real midPoint = Real(intC.midpoint());

    Real apprPoint = cosine(intC.midpoint());

    //NB: This assumes a 3rd order TS approximation
    Real coef1 = cosDer(midPoint, 1);
    Real coef2 = cosDer(midPoint, 2)/2;
    Real coef3 = cosDer(midPoint, 3)/6;
    
    rem = getCosDerRemBound(intC, intC.midpoint(), 4);
    
    Interval apprInt = Interval(apprPoint);
    
    Interval deg1Int = Interval(coef1);
    Interval deg2Int = Interval(coef2);
    Interval deg3Int = Interval(coef3);
    
    std::vector<int> deg1(numVars, 0);
    deg1[varInputInd + 1] = 1;
    std::vector<int> deg2(numVars, 0);
    deg2[varInputInd + 1] = 2;
    std::vector<int> deg3(numVars, 0);
    deg3[varInputInd + 1] = 3;

    Polynomial deg0Poly = Polynomial(Monomial(apprInt, numVars));

    Polynomial deg1Poly = Polynomial(Monomial(Interval(Real(-1) * coef1 * midPoint), numVars)) +
      Polynomial(Monomial(deg1Int, deg1));

    Polynomial deg2Poly = Polynomial(Monomial(Interval(coef2 * midPoint * midPoint), numVars)) -
      Polynomial(Monomial(Interval(Real(2) * coef2 * midPoint), deg1)) +
      Polynomial(Monomial(deg2Int, deg2));

    Polynomial deg3Poly = Polynomial(Monomial(Interval(Real(-1) * coef3 * midPoint * midPoint * midPoint), numVars)) +
      Polynomial(Monomial(Interval(Real(3) * coef3 * midPoint * midPoint), deg1)) -
      Polynomial(Monomial(Interval(Real(3) * coef3 * midPoint), deg2)) +
      Polynomial(Monomial(deg3Int, deg3));    
					
    exp = deg0Poly + deg1Poly + deg2Poly + deg3Poly;

    tmReset.expansion = exp;
    tmReset.remainder = rem;    
}

void tan_reset(TaylorModel &tmReset, const Interval intC, const int varStoreInd, const int varInputInd, const int numVars){

    Polynomial exp;
    Interval rem;

    Real midPoint = Real(intC.midpoint());

    Real apprPoint = tan(intC.midpoint());

    //NB: This assumes a 2nd order TS approximation
    Real coef1 = tan1stDer(midPoint);
    Real coef2 = tan2ndDer(midPoint)/2;
    Real coef3 = tan3rdDer(midPoint)/6;
    
    Real derBound = getTanDerBound(3, intC);
					
    Real maxDev = Real(intC.sup()) - midPoint;
    if (midPoint - Real(intC.inf()) > maxDev){
        maxDev = midPoint - Real(intC.inf());
    }

    Real fact = 6;
    maxDev.pow_assign(3);
    
    Real remainder = (derBound * maxDev) / fact;

    Interval apprInt = Interval(apprPoint);
    
    Interval deg1Int = Interval(coef1);
    Interval deg2Int = Interval(coef2);
    Interval deg3Int = Interval(coef3);

    std::vector<int> deg1(numVars, 0);
    deg1[varInputInd + 1] = 1;
    std::vector<int> deg2(numVars, 0);
    deg2[varInputInd + 1] = 2;
    std::vector<int> deg3(numVars, 0);
    deg3[varInputInd + 1] = 3;
    
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
	derBound = getTanDerBound(4, intC);
						
	remainder = (derBound * maxDev) / fact;
	remainder.to_sym_int(rem);

	Polynomial deg3Poly = Polynomial(Monomial(Interval(Real(-1) * coef3 * midPoint * midPoint * midPoint), numVars)) +
	  Polynomial(Monomial(Interval(Real(3) * coef3 * midPoint * midPoint), deg1)) -
	  Polynomial(Monomial(Interval(Real(3) * coef3 * midPoint), deg2)) +
	  Polynomial(Monomial(deg3Int, deg3));
						
	exp += deg3Poly;
    }

    tmReset.expansion = exp;
    tmReset.remainder = rem;    
    
}

void arc_reset(TaylorModel &tmReset, const Interval intC, const int varStoreInd, const int varInputInd, const int numVars){

    Polynomial exp;
    Interval rem;

    //NB: This performs a 4th order TS approximation
    int order = 4;
    int reset_type = ARC;

    Real midPoint = Real(intC.midpoint());
    
    Real apprPoint = arctan(midPoint);

    std::vector<Real> extrema_locations;
    std::vector<Real> extrema_magnitudes;

    bool left_segment_increasing = getArctanDerExtrema(extrema_locations, extrema_magnitudes, order+1);

    getGenericDerRemBound(rem, left_segment_increasing, extrema_locations,
			  extrema_magnitudes, intC, intC.midpoint(), reset_type, order+1);

    Real coef1 = arctanDer(1, midPoint);
    Real coef2 = arctanDer(2, midPoint);
    Real coef3 = arctanDer(3, midPoint);
    Real coef4 = arctanDer(4, midPoint);
    
    Interval apprInt = Interval(apprPoint);
    
    Interval deg1Int = Interval(coef1);
    Interval deg2Int = Interval(coef2);
    Interval deg3Int = Interval(coef3);
    Interval deg4Int = Interval(coef4);

    std::vector<int> deg1(numVars, 0);
    deg1[varInputInd + 1] = 1;
    std::vector<int> deg2(numVars, 0);
    deg2[varInputInd + 1] = 2;
    std::vector<int> deg3(numVars, 0);
    deg3[varInputInd + 1] = 3;
    std::vector<int> deg4(numVars, 0);
    deg4[varInputInd + 1] = 4;

    Polynomial deg0Poly = Polynomial(Monomial(apprInt, numVars));

    Polynomial deg1Poly = Polynomial(Monomial(Interval(Real(-1) * coef1 * midPoint), numVars)) +
      Polynomial(Monomial(deg1Int, deg1));

    Polynomial deg2Poly = Polynomial(Monomial(Interval(coef2 * midPoint * midPoint), numVars)) -
      Polynomial(Monomial(Interval(Real(2) * coef2 * midPoint), deg1)) +
      Polynomial(Monomial(deg2Int, deg2));

    Polynomial deg3Poly = Polynomial(Monomial(Interval(Real(-1) * coef3 * midPoint * midPoint * midPoint), numVars)) +
      Polynomial(Monomial(Interval(Real(3) * coef3 * midPoint * midPoint), deg1)) -
      Polynomial(Monomial(Interval(Real(3) * coef3 * midPoint), deg2)) +
      Polynomial(Monomial(deg3Int, deg3));

    Polynomial deg4Poly =
      Polynomial(Monomial(Interval(coef4 * midPoint * midPoint * midPoint * midPoint), numVars)) +
      Polynomial(Monomial(Interval(Real(-4) * coef4 * midPoint * midPoint * midPoint), deg1)) +
      Polynomial(Monomial(Interval(Real(6) * coef4 * midPoint * midPoint), deg2)) +
      Polynomial(Monomial(Interval(Real(-4) * coef4 * midPoint), deg3)) +
      Polynomial(Monomial(deg4Int, deg4));    
					
    exp = deg0Poly + deg1Poly + deg2Poly + deg3Poly + deg4Poly;

    tmReset.expansion = exp;
    tmReset.remainder = rem;
    
}

void arc_reset_original(TaylorModel &tmReset, const Interval intC, const int varStoreInd, const int varInputInd, const int numVars){

    Polynomial exp;
    Interval rem;

    Real midPoint = Real(intC.midpoint());
    
    Real apprPoint = arctan(intC.midpoint());


    //NB: This assumes a 2nd order TS approximation
    Real coef1 = arctanCoef(1);
    Real coef2 = arctanCoef(2);
    Real coef3 = arctanCoef(3);

    Real maxDev = Real(intC.sup()) - midPoint;
    if (midPoint - Real(intC.inf()) > maxDev){
        maxDev = midPoint - Real(intC.inf());
    }

    Real remainder = getArcTanDerBound(maxDev.getValue_RNDD(), 3);
    
    Interval apprInt = Interval(apprPoint);
    
    Interval deg1Int = Interval(coef1);
    Interval deg2Int = Interval(coef2);
    Interval deg3Int = Interval(coef3);

    std::vector<int> deg1(numVars, 0);
    deg1[varInputInd + 1] = 1;
    std::vector<int> deg2(numVars, 0);
    deg2[varInputInd + 1] = 2;
    std::vector<int> deg3(numVars, 0);
    deg3[varInputInd + 1] = 3;

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

        remainder = getArcTanDerBound(maxDev.getValue_RNDD(), 4);
	remainder.to_sym_int(rem);

	Polynomial deg3Poly = Polynomial(Monomial(Interval(Real(-1) * coef3 * midPoint * midPoint * midPoint), numVars)) +
	  Polynomial(Monomial(Interval(Real(3) * coef3 * midPoint * midPoint), deg1)) -
	  Polynomial(Monomial(Interval(Real(3) * coef3 * midPoint), deg2)) +
	  Polynomial(Monomial(deg3Int, deg3));
						
	exp += deg3Poly;
    }					

    if(rem.width() > 1){
        printf("Uncertainty too large. Please increase Taylor Model order.\n");
	exit(-1);
    }					
    
    tmReset.expansion = exp;
    tmReset.remainder = rem;
    
}

void sec_reset(TaylorModel &tmReset, const Interval intC, const int varStoreInd, const int varInputInd, const int numVars,
	       const TaylorModelVec tmvAggregation, const std::vector<Interval> doAggregation){

    Polynomial exp;
    Interval rem;

    //Real midPoint = (Real(intC.sup()) + Real(intC.inf()))/2;
    Real midPoint = Real(intC.midpoint());

    Real apprPoint = sec(midPoint);
					
    //NB: This performs a 2nd order TS approximation since higher
    //order polynomials have very large coefficients that don't work
    //well with interval analysis
    Real coef1 = sec1stDer(midPoint);
    Real coef2 = sec2ndDer(midPoint)/2;

    // Real coef4 = sec4thDer(midPoint)/24;

    //Real derBound = getSecDerBound(intC, 5);

    rem = getSecDerRemBound(intC, intC.midpoint(), 3);
    
    Interval apprInt = Interval(apprPoint);
    
    Interval deg1Int = Interval(coef1);
    Interval deg2Int = Interval(coef2);
    // Interval deg3Int = Interval(coef3);
    // Interval deg4Int = Interval(coef4);
    
    std::vector<int> deg1(numVars, 0);
    deg1[varInputInd + 1] = 1;
    std::vector<int> deg2(numVars, 0);
    deg2[varInputInd + 1] = 2;
    // std::vector<int> deg4(numVars, 0);
    // deg4[varInputInd + 1] = 4;

    /*
      Poly approx. = apprPoint + coef1 * (x - midPoint) 
      + coef2 * (x^2 - 2 * x * midPoint + midPoint^2)
      + coef3 * (x^3 - 3 * x^2 * midpoint + 3 * x * midPoint^2 - midPoint^3)
      + coef4 * (x^4 - 4 * x^3 * midpoint + 6 * x^2 * midPoint^2 - 4 * x * midPoint^3 + midPoint^4)
    */

    Polynomial deg0Poly = Polynomial(Monomial(apprInt, numVars));

    Polynomial deg1Poly = Polynomial(Monomial(Interval(Real(-1) * coef1 * midPoint), numVars)) +
      Polynomial(Monomial(deg1Int, deg1));

    Polynomial deg2Poly = Polynomial(Monomial(Interval(coef2 * midPoint * midPoint), numVars)) -
      Polynomial(Monomial(Interval(Real(2) * coef2 * midPoint), deg1)) +
      Polynomial(Monomial(deg2Int, deg2));

    exp = deg0Poly + deg1Poly + deg2Poly;

    //if uncertainty too large, use a 3rd order approximation
    if (rem.width() > 0.00001){

      printf("third order\n");

            std::vector<int> deg3(numVars, 0);
	    deg3[varInputInd + 1] = 3;

	    Real coef3 = sec3rdDer(midPoint)/6;
	    Interval deg3Int = Interval(coef3);
	    
	    rem = getSecDerRemBound(intC, intC.midpoint(), 4);
	    
	    Polynomial deg3Poly = Polynomial(Monomial(Interval(Real(-1) * coef3 * midPoint * midPoint * midPoint), numVars)) +
	      Polynomial(Monomial(Interval(Real(3) * coef3 * midPoint * midPoint), deg1)) -
	      Polynomial(Monomial(Interval(Real(3) * coef3 * midPoint), deg2)) +
	      Polynomial(Monomial(deg3Int, deg3));
	    
	    exp += deg3Poly;

    }				        

    // Polynomial deg4Poly =
    //   Polynomial(Monomial(Interval(coef4 * midPoint * midPoint * midPoint * midPoint), numVars)) +
    //   Polynomial(Monomial(Interval(Real(-4) * coef4 * midPoint * midPoint * midPoint), deg1)) +
    //   Polynomial(Monomial(Interval(Real(6) * coef4 * midPoint * midPoint), deg2)) +
    //   Polynomial(Monomial(Interval(Real(-4) * coef4 * midPoint), deg3)) +
    //   Polynomial(Monomial(deg4Int, deg4));

    printf("reset remainder: [%13.10f, %13.10f]\n", rem.inf(), rem.sup());
    
    tmReset.expansion = exp;
    tmReset.remainder = rem;

}

void div_reset(TaylorModel &tmReset, const Interval intC, const int varStoreInd, const int varDenInd, const int numVars){

    Polynomial exp;
    Interval rem;
				  
    Real midPoint = Real(intC.midpoint());

    Real apprPoint = divide(midPoint);
					
    //NB: This performs a 3rd order TS approximation
    Real coef1 = div1stDer(midPoint);
    Real coef2 = div2ndDer(midPoint)/2;
    Real coef3 = div3rdDer(midPoint)/6;
    //Real coef4 = div4thDer(midPoint)/24;
    
    rem = getDivDerRemBound(intC, intC.midpoint(), 4);
    
    Interval apprInt = Interval(apprPoint);
    
    Interval deg1Int = Interval(coef1);
    Interval deg2Int = Interval(coef2);
    Interval deg3Int = Interval(coef3);
    //Interval deg4Int = Interval(coef4);
    
    std::vector<int> deg1(numVars, 0);
    deg1[varDenInd + 1] = 1;
    std::vector<int> deg2(numVars, 0);
    deg2[varDenInd + 1] = 2;
    std::vector<int> deg3(numVars, 0);
    deg3[varDenInd + 1] = 3;
    //std::vector<int> deg4(numVars, 0);
    //deg4[varDenInd + 1] = 4;
					
    /*
      Poly approx. = apprPoint + coef1 * (x - midPoint) 
      + coef2 * (x^2 - 2 * x * midPoint + midPoint^2)
      + coef3 * (x^3 - 3 * x^2 * midpoint + 3 * x * midPoint^2 - midPoint^3)
      + coef4 * (x^4 - 4 * x^3 * midpoint + 6 * x^2 * midPoint^2 - 4 * x * midPoint^3 + midPoint^4)
    */
    
    Polynomial deg0Poly = Polynomial(Monomial(apprInt, numVars));
    
    Polynomial deg1Poly = Polynomial(Monomial(Interval(Real(-1) * coef1 * midPoint), numVars)) +
      Polynomial(Monomial(deg1Int, deg1));
    
    Polynomial deg2Poly = Polynomial(Monomial(Interval(coef2 * midPoint * midPoint), numVars)) -
      Polynomial(Monomial(Interval(Real(2) * coef2 * midPoint), deg1)) +
      Polynomial(Monomial(deg2Int, deg2));

    Polynomial deg3Poly =
      Polynomial(Monomial(Interval(Real(-1) * coef3 * midPoint * midPoint * midPoint), numVars)) +
      Polynomial(Monomial(Interval(Real(3) * coef3 * midPoint * midPoint), deg1)) -
      Polynomial(Monomial(Interval(Real(3) * coef3 * midPoint), deg2)) +
      Polynomial(Monomial(deg3Int, deg3));

    // Polynomial deg4Poly =
    //   Polynomial(Monomial(Interval(coef4 * midPoint * midPoint * midPoint * midPoint), numVars)) +
    //   Polynomial(Monomial(Interval(Real(-4) * coef4 * midPoint * midPoint * midPoint), deg1)) +
    //   Polynomial(Monomial(Interval(Real(6) * coef4 * midPoint * midPoint), deg2)) +
    //   Polynomial(Monomial(Interval(Real(-4) * coef4 * midPoint), deg3)) +
    //   Polynomial(Monomial(deg4Int, deg4));
    
					
    exp = deg0Poly + deg1Poly + deg2Poly + deg3Poly;// + deg4Poly;

    tmReset.expansion = exp;
    tmReset.remainder = rem;
    
}
