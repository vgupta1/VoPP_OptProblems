#lowerBoundsSymmetric

#computes the maximum over the geometric grid
#sep_1(a, b, l, u) solves  min_{t in [l, u]}  at + 2bt H(t)
#sep_2(a, l, u) solves min_{t in [l, u]} a H(t)
#user must ensure sep_1, sep_2 and h are consistent.
#H(t) = 1/2t int_[-t, t] h(1+s)ds
function _rev_p_lb_symmetric(S, h, sep_1, sep_2, H, H_un, pj, numCuts, print_trace, TOL)
	m = Model(solver=GurobiSolver(OutputFlag=false))
	@variable(m, theta)
	@variable(m, lam)
	@objective(m, Min, theta)

	#add endpoints to get things started
	#t = 0
	ind = pj <= 1 ? pj : 0
	@constraint(m, lam * H(0) >= ind - theta)

	#t = abs(pj -1)
	tstar = abs(pj -1)
	@constraint(m, lam * H(tstar) >= ind - theta)
	@constraint(m, (2theta - pj) * tstar + lam * H_un(tstar) >= pj * (1 - pj))

	#t = S-1
	@constraint(m, (2theta - pj) * (S - 1) + lam * H_un(S - 1) >= pj * (1 - pj))

	for iter = 1:numCuts
		status = solve(m)
		feasible = true

		print_trace && println("Iteration: \t", iter)
		#First constraint [0, abs(pj -1)]
		rhs =  ind - getvalue(theta)
		val, tstar = sep_2(getvalue(lam), 0, abs(pj -1))
		print_trace && println("Iter: $iter Val $val RHS $rhs  : tstar $tstar")
		if val < rhs - TOL
			@constraint(m, lam * H(tstar) >= ind - theta)
			feasible = false
		end

		#2nd constraint [abs(pj-1) S-1]
		rhs = pj * (1 - pj)
		val, tstar = sep_1(2 * getvalue(theta) - pj, getvalue(lam), abs(pj-1), S-1)

		print_trace && println("Iter: $iter Val $val RHS $rhs : tstar $tstar")
		if val < rhs - TOL 
			cnst = @constraint(m, (2theta - pj) * tstar + lam * H_un(tstar) >= rhs)
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
	throw("Max Iteration Count Reached:  h(1): $(h(1))")
	return -1.
end

function vopp_lb_symmetric_MAD(mu, S, M, D; delta = S/100, numCuts=100, print_trace=false, TOL=1e-6, pj=-1, safe_fail=false)
	#standardize
	c = mu * (1 - M)
	Sc = vopp.comp_Sc(S, M)

	#round down Sc if it's too large and throw otherwise
	Sc = min(Sc, 2)
	@assert Sc > 1 "Scale $S is too small for any distribution"

	if D/M > .25 
		println("Deviation $D > .25, problem likely infeasible")
		if safe_fail 
			println("Safe_fail set to true. Aborting Calc")
			return -1.0
		end
	end

	#create the functions/separators closed-form
	h(t) = M/2 * abs(t-1) - D
	H(t) = M * t / 4 - D
	H_un(t) = 2t * H(t)

	#sep_1(a, b, l, u) solves  min_{t in [l, u]}  at + b * H_un(t)
	#min_{t in [l, u]} (a - 2bD) t + bM/2 * t^2
	function sep_1(a, b, l, u)
		f(t) = (a - 2b * D) * t + b * M/2 * t^2
		valstar, tstar = f(l), l
		if f(u) < valstar
			valstar, tstar = f(u), u
		end
		if b > 0 
			tmid = (2b * D - a) / b / M
			valmid = Inf
			if l <= tmid <= u 
				valmid = f(tmid)
			end
			if valmid < valstar
				valstar, tstar = valmid, tmid
			end
		end
		valstar, tstar
	end
	
	#solves min_{t in [l, u]} a H(t)
	function sep_2(a, l, u)
		if a * H(u) < a * H(l)
			return a * H(u), u
		end
		return a * H(l), l
	end

	#search over the geometric price ladder
	if pj > 0
		ps = [pj]
	else
		ps = geom_price_ladder(Sc, delta)
	end
	rstar = -Inf
	pstar = 0. 
	rj = 0.

	for pj in ps
		#kick i tto the workhorse
		rj = _rev_p_lb_symmetric(Sc, h, sep_1, sep_2, H, H_un, pj, numCuts, print_trace, TOL)
		print_trace && println("\n Price Pj:\t", pj, "\t", "Rj:\t", rj, "\n")

		if rj > rstar
			rstar = rj 
			pstar = pj
		end
	end

	return 1/rstar, pstar	
end




