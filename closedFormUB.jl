### Closed-Form UB for the MAD Case

#Assumes distribution already standardized:  mu = 1, M = 1
delta_l(S) = -log(S) / lambertw(-exp(-1) / S)
delta_m(S) = log(S) / (1 + log(S))
delta_h(S) = (S - 1) / S 

#Assumes Standardized:  mu = 1, S = 1
function _vopp_ub_MAD_cf(S, D)
	if 0 <= D <= delta_l(S)
		return -lambertw((D - 1) / exp(1)) / (1 - D)
	elseif D <= delta_m(S)
		return log(S) / D
	elseif D <= delta_h(S)
		return -lambertw( -exp(-1) / S / (1 - D) )
	else
		throw("D $D < $(delta_h(S)) is too large for Scale $S.  Infeasible")
	end
	return -1.
end

function vopp_ub_MAD(S, M, D; method = :ClosedForm)
	if method == :ClosedForm
		Sc = vopp.comp_Sc(S, M)
		return _vopp_ub_MAD_cf(Sc, D / M)
	elseif method == :Optimization
		return vopp_ub_MAD_opt(S, M, D)[1]
	else
		throw("Method must be one of :ClosedForm or :Optimization : $method")
	end
	return -1.
end

#assumes things already standardized
function tight_dist_ub_MAD(x, S, D)
	alpha = 1/vopp_ub_MAD(S, 1., D)
	#match which regime 
	if 0 <= D <= delta_l(S)  #low heterogeneity
		if x <= alpha 
			return 1.
		elseif x <= 1
			return alpha / x
		elseif x <= S
			return D / log(S) / x
		else
			return 0.
		end
	elseif D <= delta_m(S)  #medium
		if x == 0 
			return 1
		elseif x <= exp(1) * S^(1 - 1/D)
			return alpha / exp(1) * S^(1/D - 1)
		elseif x <= S
			return alpha / x
		else
			return 0.
		end
	elseif D <= delta_h(S) #high
		if x == 0 
			return 1.
		elseif x <= alpha/(1-D)
			return 1 - D
		elseif x <= S
			return alpha / x
		else
			return 0.
		end	
	else
		throw("D not valid; $D not in [0, $(delta_h(S))")
	end
end

