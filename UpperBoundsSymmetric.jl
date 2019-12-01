#UpperBoundsSymmetric.jl
#symmetric and unimodal assumptions!!!


#generic separation formulation
#E[h(V)] = 0
#H(t) = int_[-t, t] h(1 + s) ds
#sep_fun(a, lam, l, u) optimizes
#   #max a + lam H(t) 
#   #s.t. t in [l, u]
#   returns val, tstar
function _vopp_sym_unimodal(S, h, H, sep_fun, N, numCuts, TOL, print_trace)
	#create the symmetric grid
	@assert 1 < S <= 2 "Standardized S must be between 1 and 2 by symmetry"
	ps = collect(range(1, stop=S, length=N + 1))
	ps = vcat(2 .- ps[2:end], ps)  
	sort!(ps)  #access by ix + N + 1  so that ix spans from -N:N

	m = Model(solver=GurobiSolver(OutputFlag = false))
	@variable(m, theta)
	@variable(m, lam)
	@variable(m, Q[-N:N] >= 0)
	@constraint(m, sum(Q) == 1)

	#t = 0 case
	@constraint(m, theta + lam * h(1) <= sum(ps[j + N + 1] * Q[j] for j = -N:0))
	@objective(m, Max, theta)

	for iter = 1:numCuts
		status = solve(m)
		feasible = true

		print_trace && println("Iteration: \t", iter)
		for k = 0:N-1
			a = 2 * getvalue(theta)
			a -= 2 * sum(ps[j + N + 1] * getvalue(Q[j]) for j = -N:-(k+1))
			a -= sum(ps[j + N + 1] * getvalue(Q[j]) for j = -k:k)
			rhs = sum(ps[j + N + 1] * getvalue(Q[j]) * (1 - ps[j + N + 1]) for j = -k:k)

			val, tstar = sep_fun(a, getvalue(lam), ps[k + N + 1] - 1, 
												   ps[k + 1  + N + 1] - 1)

			if val > rhs + TOL
				feasible = false
				print_trace && println("iter $iter : \t $val \t $rhs")

				@constraint(m, 2tstar * theta + lam * H(tstar) <= 
					2tstar * sum(ps[j + N + 1] * Q[j] for j = -N:-(k+1)) + 
					sum(ps[j + N + 1] * Q[j] * (1 + tstar - ps[j + N + 1]) for j=-k:k)
					)
			end
		end

		if feasible
			if status != :Optimal
				println("Problem likely infeasible.  Dual Status:\t", status)
				return -1.0
			else #correctly solved
				return 1 / getobjectivevalue(m)
			end
		end
	end
	##Exhausted iteration count
	throw("Max Iteration Count Reached")
	return -1.
end


#exposed function for MAD
function vopp_sym_unimodal_MAD(mu, S, M, D, mode; numCuts=100, N=100, print_trace=false, TOL=1e-6)
	#standardize
	c = mu * (1 - M)
	Sc = vopp.comp_Sc(S, M)
	mode_c = (mode - c) / (mu - c)

	@assert Sc >= 1 "S is too small for given mean.  Problem infeasible."
	if Sc > 2
		#just silently round down for now
		Sc = 2
	end

	#create subfunctions
	h(t) = M * abs(t - 1)/2 - D

	#H(t) = int_[-t, t] h(1 + t)
	function H(t)
		if t < 1e-10
			return 0.
		end
		quadgk(s->h(1+s), -t, 0., t)[1]
	end

	#separator
	#max a + lam H(t) 
	#s.t. t in [l, u]
	#returns val, tstar
	function sep_fun(a, lam, l, u)
		val = a * u + lam * H(u)
		tstar = u

		val_l = a * l  + lam * H(l)
		if val_l > val 
			val = val_l
			tstar = l
		end

		#compute the middle one if its reasonable.
		if lam != 0 
			tmid = 2D / M - a / lam / M 
		else
			tmid = -1.
		end

		if l <= tmid <= u 
			val_mid = a * tmid + lam * H(tmid)
		else
			val_mid = -Inf
		end

		if val_mid > val
			val = val_mid
			tstar = tmid
		end

		return val, tstar
	end

	#kick it off
	_vopp_sym_unimodal(Sc, h, H, sep_fun, N, numCuts, TOL, print_trace)
end


function vopp_sym_unimodal_CV(mu, S, M, C, mode; TOL=1e-6, N=100, numCuts=500, print_trace=false)
	#standardize
	c = mu * (1 - M)
	Sc = vopp.comp_Sc(S, M)
	mode_c = (mode - c) / (mu - c)

	#create subfunctions
	h(t) = M^2 * (t - 1)^2 - C^2

	#H(t) = int_[-t, t] h(1 + t)
	function H(t)
		if t < 1e-10
			return 0.
		end

		quadgk(s->h(1+s), -t, 0, t)[1]
	end

	#separator
	#max a + lam H(t) 
	#s.t. t in [l, u]
	#returns val, tstar
	function sep_fun(a, lam, l, u)
		val = a * u + lam * H(u)
		tstar = u

		val_l = a * l + lam * H(l)
		if val_l > val 
			val = val_l
			tstar = l
		end

		if lam != 0 
			tmid = (2lam * C^2 - a) / (2lam * M^2)
		else
			tmid = -1.
		end

		if tmid > 0
			tmid = sqrt(tmid)
		end
		if l <= tmid <= u 
			val_mid = a * tmid + lam * H(tmid)
		else 
			val_mid = -Inf
		end

		if val_mid > val
			val = val_mid
			tstar = tmid
		end

		return val, tstar
	end

	#kick it off
	_vopp_sym_unimodal(Sc, h, H, sep_fun, N, numCuts, TOL, print_trace)
end


