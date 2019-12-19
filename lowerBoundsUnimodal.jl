#lowerBoundsUnimodal

#Based on Theorem 3
function _lower_bound_unimodal(D, M)
    if D/M <= 1/3
        return 1 / (1 - D/M)
    end

    return 8D / M / (1 + D/M )^2
end

####
## Separators for the Unimodal Case
####
#min_{t in [l, u]} a * H(t) + bt
#assumes h is convex, H(t) = 1/(mode-t) Int_t^mode h
function _sep_LB_Unimodal_1(h, H, a, b, mode, l, u)
	#check endpoints.
	val_star, t_star = a * H(l) + b * l, l

	if a * H(u) + b * u < val_star
		val_star = a * H(u) + b * u
		t_star = u
	end

	if a > 0 #convex case
		#dirty hack to get around t == mode
		function deriv(t) 
			step = 1e-8
			if t == mode  #analytically = h_prime(mode)/ 2
				return a/4/step * (h(mode + step) - h(mode - step)) + b
			end
			a/(mode - t) * (H(t) - h(t)) + b
		end
		val_crit, t_crit = Inf, l
		if  deriv(l) * deriv(u) < 0 #critical point in interval
			t_crit = find_zero(deriv, (l,u), Bisection())
			val_crit = a * H(t_crit) + b * t_crit
		end
		if val_crit < val_star
			t_star = t_crit
			val_star = val_crit
		end
	end

	return val_star, t_star
end

#min_{t in [l, u]} a t^2 + b t + c int_[t, mode] h
#specifically for MAD:  h(t) = M/2 * abs(t-1) - D
#Notice this is specialized ot MAD bc otherwise optimization might be tricky
function _sep_LB_Unimodal_MAD_2(M, D, a, b, c, mode, l, u)
	h(t) = M/2 * abs(t-1) - D
	H_un(t) = quadgk(h, t, mode)[1]

	#check endpoints first
	t_star = l
	val_star = a * t_star^2 + b * t_star + c * H_un(t_star) 

	val_u = a * u^2 + b * u + c * H_un(u) 
	if val_u < val_star
		val_star = val_u
		t_star = u
	end

	#the <=1 critical point
	t_less = (-2b - 2c * D + c * M) / (4a + c * M)
	val_less = +Inf
	if t_less <=1 && l <= t_less <= u
		val_less = a * t_less^2 + b * t_less + c * H_un(t_less) 
	end
	if val_less < val_star
		val_star = val_less
		t_star = t_less
	end

	#the >= 1 critcical point
	t_greater = (-2b  - 2c * D - c * M) / (4a - c * M)
	val_greater = Inf
	if t_greater >= 1 && l <= t_greater <= u
		val_greater = a * t_greater^2 + b * t_greater + c * H_un(t_greater) 
	end
	if val_greater < val_star
		val_star = val_greater
		t_star = t_greater
	end

	return val_star, t_star
end


#min_{t in [l, u]} a t^2 + b t + c int_[t, mode] h
#specifically for CV:  h(t) = M^2 (t-1)^2 - C^2
#Notice this is specialized to CV bc otherwise optimization might be tricky
function _sep_LB_Unimodal_CV_2(M, C, a, b, c, mode, l, u)
	#define the functions 
	h(t) = M^2 * (t-1)^2 - C^2
	H_unc(t) = quadgk(h, t, mode)[1]
	f(t) = a * t^2 + b*t + c * H_unc(t)

	#check end points
	t_star = l
	val_star = f(t_star)

	if f(u) < val_star
		t_star = u
		val_star = f(u)
	end

	#compute critical points and check feasibility
	disc = a^2 + c * (-2a + b - C^2) * M^2
	if disc >= 0
		t1 = a - c*M^2 - sqrt(disc)
		t1 /= c * M^2
		if l <= t1 <= u && f(t1) < val_star
			val_star = f(t1)
			t_star = t1
		end

		t2 = a - c*M^2 + sqrt(disc)
		t2 /= c*M^2
		if l <= t2 <= u && f(t2) < val_star
			val_star = f(t2)
			t_star = t2
		end
	end

	return val_star, t_star
end


#####################
#computes the maximum over the geometric grid
#sep_1(a, b, l, u) solves  min_{t in [l, u]}  a H(t) + b t
#sep_2(a, b, c, l, u) solves min_{t in [l, u]} a t^2 + b t + c int_{t,m} h 
#user must ensure sep_1, sep_2 and h are consistent.
function _vopp_lb_unimodal_(S, mode, h, delta, sep_1, sep_2, numCuts, TOL, print_trace, ps)
	#Integrated functions
	function H_un(t)
		if min(mode, t) <= 1 <= max(mode, t)
			return quadgk(h, t, 1, mode)[1]
		else
			return quadgk(h, t, mode)[1]
		end
	end
	H(t) = t == mode ? h(mode) : H_un(t) / (mode - t)

	rstar = -Inf
	pstar = 0. 
	rj = 0.

	for pj in ps
		#dispatch to appropriate case
		if pj < mode
			rj = _rev_p_lb_unimodal_case1(S, mode, h, sep_1, sep_2, H, H_un, pj, numCuts, print_trace, TOL)
			print_trace && println("\n Price Pj:\t", pj, "\t", "Rj:\t", rj, "\n")
		elseif pj == mode 
			rj = _rev_p_lb_unimodal_case2(S, mode, h, sep_1, sep_2, H, H_un, pj, numCuts, print_trace, TOL)
			print_trace && println("\n Price Pj:\t", pj, "\t", "Rj:\t", rj, "\n")
		else
			rj = _rev_p_lb_unimodal_case3(S, mode, h, sep_1, sep_2, H, H_un, pj, numCuts, print_trace, TOL)
			print_trace && println("\n Price Pj:\t", pj, "\t", "Rj:\t", rj, "\n")
		end
		if rj > rstar
			rstar = rj 
			pstar = pj
		end
	end
	return 1/rstar, pstar
end

#computes the maximum over the geometric grid
#sep_1(a, b, l, u) solves  min_{t in [l, u]}  a H(t) + b t
#sep_2(a, b, c, l, u) solves min_{t in [l, u]} a t^2 + b t + c int_{t,m} h 
#user must ensure sep_1, sep_2 and h are consistent.
#H_un(t) = int_[t, m] h.  H(t) = H_un(t)/(m-t)
function _rev_p_lb_unimodal_case1(S, mode, h, sep_1, sep_2, H, H_un, pj, numCuts, print_trace, TOL)
	@assert pj < mode "Wrong Case called for LB Unimodal: price point: $pj mode: $mode"

	m = Model(solver=GurobiSolver(OutputFlag=false))
	@variable(m, theta)
	@variable(m, lam1)
	@variable(m, lam2)
	@objective(m, Min, theta + (2 - mode) * lam2)

	#add endpoints to get started
	#t = 0
	@constraint(m, theta * mode  + lam1 * H_un(0) >= pj * (mode - pj))
	#t = pj
	@constraint(m, theta + lam1 * H(pj) + lam2 * pj >= pj)
	#t = S
	@constraint(m, theta + lam1 * H(S) + lam2 * S >= pj)

	for iter = 1:numCuts
		status = solve(m)
		feasible = true

		print_trace && println("Iteration: \t", iter)
		#1st constraint corresponds to second separator
		a = -getvalue(lam2)
		b = getvalue(lam2) * mode - getvalue(theta)
		c = getvalue(lam1)
		rhs = pj * (mode - pj) - getvalue(theta) * mode

		val, tstar = sep_2(a, b, c, 0, pj)
		print_trace && println("Iter: $iter Val $val RHS $rhs  : tstar $tstar")
		if val < rhs - TOL
			@constraint(m, theta * (mode - tstar) + 
							lam1 * H_un(tstar) + 
							lam2 * tstar * (mode - tstar) >= pj * (mode - pj)
							)
			feasible = false
		end

		#2nd constraint corresponds to first separator
		val, tstar = sep_1(getvalue(lam1), getvalue(lam2), pj, S)
		rhs = pj - getvalue(theta)
		print_trace && println("Iter: $iter Val $val RHS $rhs : tstar $tstar")
		if val < rhs - TOL 
			cnst = @constraint(m, theta + lam1 * H(tstar) + lam2 * tstar >= pj)
			feasible = false
		end

		if feasible
			if status != :Optimal
				println("Problem likely infeasible.  Dual Status:\t", status)
				return -1.0
			else #correctly solved
				return getobjectivevalue(m)
			end
		end
	end
	##Exhausted iteration count
	throw("Max Iteration Count Reached:  h(1): $(h(1)) ")
	return -1.
end
								
function _rev_p_lb_unimodal_case2(S, mode, h, sep_1, sep_2, H, H_un, pj, numCuts, print_trace, TOL)
	@assert abs(pj - mode)<= 1e-10 "Wrong Case called for LB Unimodal: price point: $pj mode: $mode"

	m = Model(solver=GurobiSolver(OutputFlag=false))
	@variable(m, theta)
	@variable(m, lam1)
	@variable(m, lam2)
	@objective(m, Min, theta + (2 - mode) * lam2)

	#add endpoints to get started
	#t = 0
	@constraint(m, theta + lam1 * H(0) >= 0)
	#t = mode
	@constraint(m, theta + lam1 * H(mode) + lam2 * mode >= mode)
	#t = S
	@constraint(m, theta + lam1 * H(S) + lam2 * S >= S)

	for iter = 1:numCuts
		status = solve(m)
		feasible = true

		print_trace && println("Iteration: \t", iter)
		#1st constraint corresponds to first opt type
		val, tstar = sep_1(getvalue(lam1), getvalue(lam2), 0, mode)
		rhs = 0
		if val < rhs - TOL 
			@constraint(m, theta + lam1 * H(tstar) + lam2 * tstar >= 0)
			feasible = false
		end

		#2nd constraint corresponds ot second type
		rhs = mode

		val, tstar = sep_1(getvalue(lam1), getvalue(lam2), mode, S)
		if val < rhs - TOL
			@constraint(m, theta + lam1 * H(tstar) + lam2 * tstar >= mode)
			feasible = false
		end

		if feasible
			if status != :Optimal
				println("Problem likely infeasible.  Dual Status:\t", status)
				return -1.0
			else #correctly solved
				return getobjectivevalue(m)
			end
		end
	end
	##Exhausted iteration count
	throw("Max Iteration Count Reached")
	return -1.
end
							
function _rev_p_lb_unimodal_case3(S, mode, h, sep_1, sep_2, H, H_un, pj, numCuts, print_trace, TOL)
	@assert pj > mode "Wrong Case called for LB Unimodal: price point: $pj mode: $mode"

	m = Model(solver=GurobiSolver(OutputFlag=false))
	@variable(m, theta >= -300)
	@variable(m, lam1 >= -400)
	@variable(m, lam2 >= -500)
	@objective(m, Min, theta + (2 - mode) * lam2)

	#add endpoints to get started
	#t = 0
	@constraint(m, theta  + lam1 * H(0) >= 0)
	#t = pj
	@constraint(m, theta + lam1 * H(pj) + lam2 * pj >= 0)
	#t = S
	@constraint(m, theta * (S - mode) - lam1 * H_un(S) + lam2 * S * (S - mode) >= pj * (S - pj))

	for iter = 1:numCuts
		status = solve(m)
		feasible = true

		print_trace && println("Iteration: \t", iter, "\t Status: ", status)
		#1st constraint corresponds to first separator
		val, tstar = sep_1(getvalue(lam1), getvalue(lam2), 0, pj)
		rhs = -getvalue(theta)
		print_trace && println("Iter: val ", val, " rhs ", rhs, " : ", tstar)

		if val < rhs - TOL 
			cnst = @constraint(m, theta + lam1 * H(tstar) + lam2 * tstar >= 0)
			feasible = false
		end

		#2nd constraint corresponds to second separator
		a = getvalue(lam2)
		b = getvalue(theta) - getvalue(lam2) * mode - pj
		c = -getvalue(lam1)
		rhs = getvalue(theta) * mode - pj^2

		val, tstar = sep_2(a, b, c, pj, S)
		print_trace && println("Iter: val ", val, " rhs ", rhs, " : ", tstar)

		if val < rhs - TOL
			@constraint(m, theta * (tstar - mode) - 
							lam1 * H_un(tstar) + 
							lam2 * tstar * (tstar - mode) >= pj * (tstar - pj)
							)
			feasible = false
		end

		if feasible
			if status != :Optimal
				println("Problem likely infeasible.  Dual Status:\t", status)
				return -1.0
			else #correctly solved
				return getobjectivevalue(m)
			end
		end
	end
	##Exhausted iteration count
	throw("Max Iteration Count Reached")
	return -1.
end

#safe_fail forces a return of -1 if D seems infeasible
function vopp_lb_unimodal_MAD(mu, S, M, D, mode; 
			method=:Opt, delta = S/100, numCuts=100, print_trace=false, TOL=1e-6, safe_fail=false, pj = -1)
	#standardize
	c = mu * (1 - M)
	Sc = vopp.comp_Sc(S, M)
	mode_c = (mode - c) / (mu - c)

	if D/M > .25 
		println("Deviation $D > .25, problem likely infeasible")
		if safe_fail 
			println("Safe_fail set to true. Aborting Calc")
			return -1.0
		end
	end

	if method == :Formula
		return _lower_bound_unimodal(D, M)
	end
	if method != :Opt
		throw("Method must be one of :Opt or :Formula")
	end

	#create the separators
	h(t) = M/2 * abs(t-1) - D
	H(t) = mode_c == t ? h(mode_c) : quadgk(h, t, mode_c)[1] / (mode_c - t)
	sep_1(a, b, l, u) = _sep_LB_Unimodal_1(h, H, a, b, mode_c, l, u)
	sep_2(a, b, c, l, u) = _sep_LB_Unimodal_MAD_2(M, D, a, b, c, mode_c, l, u)

	if pj > 0
		ps = [pj]
	else
		#search over the a geometric price ladder
		ps = geom_price_ladder(S, delta)
	end

	#kick it to the workhorse
	_vopp_lb_unimodal_(Sc, mode_c, h, delta, sep_1, sep_2, numCuts, TOL, print_trace, ps)
end

function vopp_lb_unimodal_CV(mu, S, M, C, mode; delta = S/100, numCuts=100, print_trace=false, TOL=1e-6, pj = -1)
	#standardize
	c = mu * (1 - M)
	Sc = vopp.comp_Sc(S, M)
	mode_c = (mode - c) / (mu - c)

	if C > max_cv_unimodal_guess(mu, S, M, mode)
		println("Warning:  Desired Coeff of Variation might exceed max-possible for uniform distributions.")
	end

	if pj > 0
		ps = [pj]
	else
		#search over the a geometric price ladder
		ps = geom_price_ladder(S, delta)
	end

	#create the separators
	h(t) = M^2 * (t - 1)^2 - C^2
	H(t) = mode_c == t ? h(mode_c) : quadgk(h, t, mode_c)[1] / (mode_c - t)
	sep_1(a, b, l, u) = _sep_LB_Unimodal_1(h, H, a, b, mode_c, l, u)
	sep_2(a, b, c, l, u) = _sep_LB_Unimodal_CV_2(M, C, a, b, c, mode_c, l, u)

	#kick it to the workhorse
	_vopp_lb_unimodal_(Sc, mode_c, h, delta, sep_1, sep_2, numCuts, TOL, print_trace, ps)
end


