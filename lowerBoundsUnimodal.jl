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



#min_{t in [l, u]} a t^2 + b t + c int_[t, mode] h
#specifically for GM:  h(t) = -log(M ( t-1) + 1) + log(Bbar)
#Notice this is specialized to GM bc otherwise optimization might be tricky
#Bbar  = B/mu
function _sep_LB_Unimodal_GM_2(M, Bbar, a, b, c, mode, l, u)
	#define the functions 
	h(t) = -log(M * (t - 1) + 1 ) + log(Bbar)
	H_unc(t) = quadgk(h, t, mode)[1]
	f(t) = a * t^2 + b * t + c * H_unc(t)

	#check end points
	t_star, val_star = l, f(l)

	if f(u) < val_star
		t_star, val_star = u, f(u)
	end

	#compute critical points and check feasibility
	#Degenerate case
	if c == 0.
		tp = -b /2/a
		if l <= tp <= u && f(tp) < val_star
			t_star, val_star = tp, f(tp)
		end
		return val_star, t_star
	end

	#else real work to do 
	inner = 2a/c/M * Bbar * exp(2a/c * (1/M - 1) - b/c)

	#no interior solution
	if inner < -1/MathConstants.e
		return val_star, t_star
	end

	if inner < 0
		#check the -1 branch
		tp = 1 - 1/M + c/2/a * lambertw(inner)
		if l <= tp <= u && f(tp) < val_star
			val_star, t_star = f(tp), tp
		end
	end
	#always check the principal branch
	tp = 1 - 1/M + c/2/a * lambertw0(inner)
	if l <= tp <= u && f(tp) < val_star
		val_star, t_star = f(tp), tp
	end

	return val_star, t_star
end





#####################
#Convenience call that integrates h and calls workhorse
#INtegration routines can introduce numerical error. Better ot use other method when integrals known analytically
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
	return _vopp_lb_unimodal_(S, mode, h, H, H_un, delta, sep_1, sep_2, numCuts, TOL, print_trace, ps)
end


#computes the maximum over the geometric grid
#sep_1(a, b, l, u) solves  min_{t in [l, u]}  a H(t) + b t
#sep_2(a, b, c, l, u) solves min_{t in [l, u]} a t^2 + b t + c int_{t,m} h 
#user must ensure sep_1, sep_2 and h are consistent.
#both return val, tstar
#Note that H is normalized in sep 1 and the integral is UNnormalized in sep2
function _vopp_lb_unimodal_(S, mode, h, H, H_un, delta, sep_1, sep_2, numCuts, TOL, print_trace, ps)
	rstar = -Inf
	pstar = 0. 
	rj = 0.

	for pj in ps
		#dispatch to appropriate case
		print_trace && println("\n Price Pj:\t", pj)
		if pj < mode
			rj = _rev_p_lb_unimodal_case1(S, mode, h, sep_1, sep_2, H, H_un, pj, numCuts, print_trace, TOL)
			print_trace && println("\t \t", "Rj:\t", rj, "\n")
		elseif pj == mode 
			rj = _rev_p_lb_unimodal_case2(S, mode, h, sep_1, sep_2, H, H_un, pj, numCuts, print_trace, TOL)
			print_trace && println("\t \t", "Rj:\t", rj, "\n")
		else
			rj = _rev_p_lb_unimodal_case3(S, mode, h, sep_1, sep_2, H, H_un, pj, numCuts, print_trace, TOL)
			print_trace && println("\t \t", "Rj:\t", rj, "\n")
		end
		if rj > rstar
			rstar = rj 
			pstar = pj
		end
	end
	return 1/rstar, pstar
end


#called by _vopp_lb_unimodal_
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


	#Add a silly constraint to prevent unboundedness
	@constraint(m, theta + (2 - mode) * lam2 >= -100)

	for iter = 1:numCuts
		status = solve(m)
		feasible = true

		print_trace && println("\n   Iteration: \t", iter, "\t Status:\t", status)

		#1st constraint corresponds to second separator
		a = -getvalue(lam2)
		b = getvalue(lam2) * mode - getvalue(theta)
		c = getvalue(lam1)
		rhs = pj * (mode - pj) - getvalue(theta) * mode

		val, tstar = sep_2(a, b, c, 0, pj)
		print_trace && println("\t Iter: $iter  \t Sep2: LHS-RHS >=? 0:  $(val-rhs) : tstar $tstar")

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
		print_trace && println("\t Iter: $iter  \t Sep1: LHS-RHS >=? 0:  $(val-rhs) : tstar $tstar")

		if val < rhs - TOL 
			###DEBUG
			###VG Confirm that the constraint is actually violated right now?
			temp_lhs = getvalue(lam1) * H(tstar) + getvalue(lam2) * tstar 
			temp_rhs = pj - getvalue(theta) 
			print_trace && println("Violation check: tstar $tstar \n", 
					"LHS >=? RHS: \t", temp_lhs, "\t", temp_rhs, "\t", temp_lhs - temp_rhs)

			print_trace && println("Violation Check 2\n", 
				"val:\t $val", "temp_lhs:\t", temp_lhs)
			### END DEBUG


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
	print_trace && println("Max Iteration Count Reached. \n", 
							"theta:\t", getvalue(theta), 
							"\t lam1:\t", getvalue(lam1), 
							"\t lam2:\t", getvalue(lam2))
	throw("Max Iteration Count Reached")
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
	@constraint(m, theta + lam1 * H(S) + lam2 * S >= mode)

	#A trivial (bad) constraint for unboundedness?
	@constraint(m, theta + (2-mode) * lam2 >= -100)

	for iter = 1:numCuts
		status = solve(m)
		feasible = true

		print_trace && println("\n   Iteration: \t", iter, "Status:\t", status)
		#1st constraint corresponds to first opt type
		val, tstar = sep_1(getvalue(lam1), getvalue(lam2), 0, mode)
		rhs = -1 * getvalue(theta)

		print_trace && println("\t Iter: $iter  \t Sep1: LHS -RHS >=? 0:  $(val-rhs) :\t tstar $tstar")

		if val < rhs - TOL 
			@constraint(m, theta + lam1 * H(tstar) + lam2 * tstar >= 0)
			feasible = false
		end

		rhs = mode - getvalue(theta)
		val, tstar = sep_1(getvalue(lam1), getvalue(lam2), mode, S)

		print_trace && println("\t Iter: $iter \t Sep1: LHS - RHS ?>= 0:  $(val-rhs) :\t tstar $tstar")

		if val < rhs - TOL
			@constraint(m, theta + lam1 * H(tstar) + lam2 * tstar >= mode)
			feasible = false
		end

		if feasible
			if status != :Optimal
				println("Problem likely infeasible.  Dual Status:\t", status)
				return -1.0
			else #correctly solved

				###DEBUG 
				if getobjectivevalue(m) < 0
					print_trace && println("Weird Objective\n", 
						"theta:\t", getvalue(theta), 
						"\t lam1:\t", getvalue(lam1), 
						"\t lam2:\t", getvalue(lam2))
				end

				return getobjectivevalue(m)
			end
		end
	end
	##Exhausted iteration count
	print_trace && println("Max Iteration Count Reached. \n", 
							"theta:\t", getvalue(theta), 
							"\t lam1:\t", getvalue(lam1), 
							"\t lam2:\t", getvalue(lam2))
	throw("Max Iteration Count Reached")
	return -1.
end
							
function _rev_p_lb_unimodal_case3(S, mode, h, sep_1, sep_2, H, H_un, pj, numCuts, print_trace, TOL)
	@assert pj > mode "Wrong Case called for LB Unimodal: price point: $pj mode: $mode"

	m = Model(solver=GurobiSolver(OutputFlag=false))
	@variable(m, theta)
	@variable(m, lam1)
	@variable(m, lam2)

	@objective(m, Min, theta + (2 - mode) * lam2)

	#add endpoints to get started
	#t = 0
	@constraint(m, theta  + lam1 * H(0) >= 0)
	#t = pj
	@constraint(m, theta + lam1 * H(pj) + lam2 * pj >= 0)
	#t = S
	@constraint(m, theta * (S - mode) - lam1 * H_un(S) + lam2 * S * (S - mode) >= pj * (S - pj))

	#VG Try to prevent unbounded errors.
	@constraint(m, theta + (2 - mode) * lam2 >= -10.)

	for iter = 1:numCuts
		status = solve(m)
		feasible = true

		print_trace && println("\nIteration: \t", iter, "\t Status: ", status)
		#1st constraint corresponds to first separator
		val, tstar = sep_1(getvalue(lam1), getvalue(lam2), 0, pj)
		rhs = -getvalue(theta)
		print_trace && println("\t Iter: $iter Sep1 val -rhs ?>= 0: ", val-rhs, " : ", tstar)

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
		print_trace && println("\t Iter: $iter Sep 2 val-rhs ?>= 0: ", val-rhs, " : ", tstar)

		if val < rhs - TOL
			# print_trace && println("\t VG Extra Detailed Debug")
			# print_trace && println("\t", getvalue(theta), "\t", getvalue(lam1), "\t", getvalue(lam2) )
			# print_trace && println("\t tstar:\t", tstar, "\t H_un", H_un(tstar), "\t mode\t", mode)

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
		println("Deviation $(D/M) > .25, problem likely infeasible")
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
		ps = [ (pj - c)/(mu - c) ]
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

	if C > max_cv(S, M)
		println("Warning:  Desired Coeff of Variation might exceed max-possible for uniform distributions.")
	end

	if pj > 0
		ps = [ (pj - c)/(mu - c) ]
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


function vopp_lb_unimodal_GM(mu, S, M, B, mode; delta = S/100, numCuts=100, print_trace=false, TOL=1e-6, pj = -1)
	#standardize
	c = mu * (1 - M)
	Sc = vopp.comp_Sc(S, M)
	mode_c = (mode - c) / (mu - c)
	Bbar = B/mu

	if Bbar > 1 || B < min_gm(S, M, mu) 
		println("Warning:  Desired Geometric Mean might be infeasible.  $B not in [$(min_gm(S, M, mu)), 1]")
	end

	if pj > 0
		ps = [ (pj - c)/(mu - c) ]
	else
		#search over the a geometric price ladder
		ps = geom_price_ladder(S, delta)
	end

	#create the separators
	h(t) = -log(M * (t - 1) + 1) + log(Bbar)
	#notice scaled
	H(t) = mode_c == t ? h(mode_c) : quadgk(h, t, mode_c)[1] / (mode_c - t)

	#First separator works because of convexity
	sep_1(a, b, l, u) = _sep_LB_Unimodal_1(h, H, a, b, mode_c, l, u)
	sep_2(a, b, c, l, u) = _sep_LB_Unimodal_GM_2(M, Bbar, a, b, c, mode_c, l, u)

	#kick it to the workhorse
	_vopp_lb_unimodal_(Sc, mode_c, h, delta, sep_1, sep_2, numCuts, TOL, print_trace, ps)
end


function vopp_lb_unimodal_IC(mu, S, M, phat, q, mode; delta = S/100, numCuts=100, print_trace=false, TOL=1e-6, pj = -1)
	#standardize
	c = mu * (1 - M)
	Sc = vopp.comp_Sc(S, M)
	mode_c = (mode - c) / (mu - c)
	v0 = (phat - 1)/ M + 1

	@assert c <= phat <= S "Price outside valuation bounds:\t $phat not in [$c , $S ]"
	@assert 0 <= q <= 1 "Market fraction not $q not in [0, 1]"

	if pj > 0
		ps = [(pj - c)/(mu - c)]
	else
		#search over the a geometric price ladder
		ps = geom_price_ladder(S, delta)
		ps = [ps; range(max(0, phat/2), stop=min(2 * phat, 1), length=100)]
		push!(ps, v0)
		push!(ps, mode_c)
		sort!(ps)
	end

	#create the separators
	h(t) = t >= v0 ? 1 - q : -q
	#notice scaled
	function H(t)  #G(v_0, mode_c, t) - q
    	if v0 > max(mode_c, t)  
        	return -q
    	elseif v0 < min(mode_c, t)  
        	return 1 - q
    	elseif mode_c == t
        	return mode_c >= v0 ? 1 - q : -q
    	end
    	#non-trivial answer only if v0 in interval [mode_c, t]
    	return (max(mode_c, t) - v0) / abs(mode_c - t) - q
	end

	#sep_1(a, b, l, u) solves  min_{t in [l, u]}  a H(t) + b t
	function sep_1(a, b, l, u)
		f(t) = a * H(t) + b * t
		#check endpoints
		val, tstar = f(l), l
		if f(u) < val
			val, tstar = f(u), u
		end

		#check mode if its feasible
		if l <= mode_c <= u && f(mode_c) < val
			val, tstar = f(mode_c), mode_c
		end

		#check v0 if its feasible
		if l <= v0 <= u && f(v0) < val
			val, tstar = f(v0), v0
		end

		if a * b < 0
			if v0 > mode_c  
				t1 = mode_c + sqrt( a/b * (mode_c - v0) )
				if l <= t1 <= u && f(t1) < val
					val, tstar = f(t1), t1
				end
			end
			if v0 < mode_c 
				t4 = mode_c - sqrt( a/b * (v0 - mode_c) )
				if l <= t4 <= u && f(t4) < val
					val, tstar = f(t4), t4
				end
			end
		end  #end check for internal optimizers
		return val, tstar
	end  #ends def of separator

	#sep_2(a, b, c, l, u) solves min_{t in [l, u]} a t^2 + b t + c int_{t,m} h 
	#Notice this reverses rolls of a, b from document
	
	H_un(t) = H(t) * (mode_c - t)
	function sep_2(a, b, c, l, u)
		obj(t) = a * t^2 + b * t + c * H_un(t)

		#check endpoints
		tstar, val = l, obj(l)
		if obj(u) < val
			tstar, val = u, obj(u)
		end

		#handle degenerate case where c == 0
		if c == 0.
			tp = -b/2/a		
			if l <= tp <= u && obj(tp) < val
					val, tstar = obj(tp), tp
			end
			return val, tstar
		end

		#else real work to be done
		#check the two possible critical points
		if a != 0
			t1 = c/2/a * (1 - q - b/c)
			if t1 >= v0 && l <= t1 < u && obj(t1) < val
					val, tstar = obj(t1), t1
			end

			t0 = -c/2/a * (q + b/c)
			if t0 < v0 && l <= t0 < u && obj(t0) < val
					val, tstar = obj(t0), t0
			end
		end
		return val, tstar
	end  #end def of separator

	#kick it to the workhorse
	_vopp_lb_unimodal_(Sc, mode_c, h, H, H_un, delta, sep_1, sep_2, numCuts, TOL, print_trace, ps)
end





