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
		throw("D $D is too large for Scale $S.  Infeasible")
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