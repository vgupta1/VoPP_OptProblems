#UpperBoundsUnimodal.jl

#####
#Function assumes that we have already translated to V_c space
###
#h 1d function that represents  E[h(V)] = 0
#H is a 1d function that represents int_[t, mode] h  (notice ordered limits)
#sep_fun(a, b, c, mode, pk, p_{k+1}) returns maximizers of 
#	max_{t in [pk, p_{k+1}}  at + b t^2 + c H(t)
#returns val, tstar
function _vopp_unimodal(S, mode, h, H, sep_fun; N = 100, numCuts=100, TOL=1e-6, print_trace=false)
	m = Model(solver=GurobiSolver(OutputFlag = false))
	@variable(m, theta)
	@variable(m, lam[1:2])

	#build a grid containing mode
	ps = collect(range(0, stop=S, length=N+1)) #doc indexes from 0 to N
	if !(mode in ps)
		push!(ps, mode)
		sort!(ps)
	end
	N = length(ps) - 1
	jstar = findall(mode .== ps)[1]

	@variable(m, Q[1:N+1] >= 0)  #doc indexes from 0 to N
	@constraint(m, sum(Q) == 1.)

	#the constraint corresponding to t = m
	@constraint(m, theta + mode * lam[1] + h(mode) * lam[2] <= 
				sum(ps[j] * Q[j] for j = 1:jstar)
				)

	@objective(m, Max, theta + (2 - mode) * lam[1])


	for iter = 1:numCuts
		status = solve(m)
		feasible = true

		print_trace && println("Starting Iteration: \t", iter, "\t", feasible)
		#t < mode constraints
		for k = 1:jstar-1
			a = sum(ps[j] * getvalue(Q[j]) for j = 1:k)
			a = a - getvalue(theta) + mode * getvalue(lam[1])
			b = -getvalue(lam[1])
			c =  getvalue(lam[2])

			val, tstar = sep_fun(a, b, c, mode, ps[k], ps[k+1])
			rhs = mode * sum(ps[j] * getvalue(Q[j]) for j = 1:k)
			rhs += sum(ps[j] * getvalue(Q[j]) * (mode - ps[j]) for j = (k+1):jstar)
			rhs -= mode * getvalue(theta)

			if val > rhs + TOL
				feasible = false
				print_trace && println("< iter $iter : \t $val \t $rhs")

				@constraint(m, (mode - tstar) * theta + (mode - tstar) * tstar * lam[1] + 
								H(tstar) * lam[2] <= sum(ps[j] * Q[j] * (mode - tstar) for j = 1:k) + 
												     sum(ps[j] * Q[j] * (mode - ps[j]) for j = (k+1):jstar)

					)
			end
		end
		print_trace && feasible && println("Iter $iter \t All t < m cnsts satisfied")

		#t > mode constraints
		for k = (jstar+1):N
			a = getvalue(theta) - mode * getvalue(lam[1])
			a -= sum(ps[j] * getvalue(Q[j]) for j = 1:k)
			b = getvalue(lam[1])
			c = -getvalue(lam[2])

			val, tstar = sep_fun(a, b, c, mode, ps[k], ps[k+1])
			rhs = mode * getvalue(theta)
			rhs -= mode * sum(ps[j] * getvalue(Q[j]) for j = 1:jstar)
			rhs -= sum(ps[j]^2 * getvalue(Q[j]) for j = (jstar + 1):k)

			if val > rhs + TOL
				feasible = false
				print_trace && println("> iter $iter : \t $val \t $rhs")

				@constraint(m, (tstar - mode) * theta + (tstar- mode) * tstar * lam[1] - 
								H(tstar) * lam[2] <= sum(ps[j] * Q[j] * (tstar - mode) for j = 1:jstar) + 
												     sum(ps[j] * Q[j] * (tstar - ps[j]) for j = (jstar+1):k)

					)
			end
		end
		print_trace && feasible && println("Iter $iter \t All t > m cnsts satisfied")

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

function vopp_ub_unimodal_MAD(mu, S, M, D, mode; numCuts=100, N=100, print_trace=false)
	#standardize
	c = mu * (1 - M)
	Sc = vopp.comp_Sc(S, M)
	mode_c = (mode - c) / (mu - c)

	#create subfunctions
	h(t) = M * abs(t - 1)/2 - D

	#H(t) = int_[t, mode_c] h
	function H(t)
		quadgk(h, t, 1, mode_c)[1]
	end

	#return maximizers of 
	#max_{[pk, pk_p]}  at + bt^2 + c H(t)
	function sep_fun(a, b, c, mode, pk, pk_p)
		obj(t) = a * t + b * t^2 + c * H(t)
		val_k = obj(pk)
		val_k_p = obj(pk_p)

		tstar, val = 0., 0.
		if val_k > val_k_p
			tstar, val = pk, val_k
		else
			tstar, val = pk_p, val_k_p
		end

		t_less = (-2a - 2D + c * M) / (4b + c * M)
		if t_less <= 1 && pk <= t_less <= pk_p
			val_less = obj(t_less)
			if val_less >= val
				tstar, val = t_less, val_less
			end
		end

		t_more = (-2a - 2D - c * M) / (4b - c * M)
		if t_more >= 1 && pk <= t_more <= pk_p
			val_more = obj(t_more)
			if val_more >= val
				tstar, val = t_more, val_more
			end
		end
		return val, tstar
	end

	#Now kick it off!	   
	return _vopp_unimodal(Sc, mode_c, h, H, sep_fun, N=N, numCuts=numCuts, print_trace=print_trace)
end


function vopp_ub_unimodal_CV(mu, S, M, C, mode; N=100, numCuts=500, print_trace=false)
	#standardize
	c = mu * (1 - M)
	Sc = vopp.comp_Sc(S, M)
	mode_c = (mode - c) / (mu - c)

	#create subfunctions
	h(t) = M^2 * (t - 1)^2 - C^2

	#H(t) = int_[t, mode_c] h
	function H(t)
		quadgk(h, t, 1, mode_c)[1]
	end

	#return maximizers of 
	#max_{[pk, pk_p]}  at + bt^2 + c H(t)
	function sep_fun(a, b, c, mode, pk, pk_p)
		obj(t) = a * t + b * t^2 + c * H(t)
		val_k = obj(pk)
		val_k_p = obj(pk_p)

		tstar, val = 0., 0.
		if val_k > val_k_p
			tstar, val = pk, val_k
		else
			tstar, val = pk_p, val_k_p
		end

		tp = b + c*M^2 + sqrt(b^2 + c *(a + 2b + c*C^2) * M^2)
		tp /= c*M^2

		if pk <= tp <= pk_p
			valp = obj(tp)
			if valp >= val
				tstar, val = tp, valp
			end
		end

		tm = b + c*M^2 - sqrt(b^2 + c *(a + 2b + c*C^2) * M^2)
		tm /= c*M^2

		if pk <= tm <= pk_p
			valm = obj(tm)
			if valm >= val
				tstar, val = tm, valm
			end
		end
		return val, tstar
	end

	#Now kick it off!	   
	return _vopp_unimodal(Sc, mode_c, h, H, sep_fun, N=N, numCuts=numCuts, print_trace=print_trace)
end

function vopp_ub_unimodal_GM(mu, S, M, B, mode; N=100, numCuts=500, print_trace=false)
	#standardize
	c = mu * (1 - M)
	Sc = vopp.comp_Sc(S, M)
	mode_c = (mode - c) / (mu - c)

	#create subfunctions
	h(t) = -log(M * t + 1 - M) + log(B/mu)

	#H(t) = int_[t, mode_c] h
	function H(t)
		quadgk(h, t, 1, mode_c)[1]
	end

	#return maximizers of 
	#max_{[pk, pk_p]}  at + bt^2 + c H(t)
	function sep_fun(a, b, c, mode, pk, pk_p)
		obj(t) = a * t + b * t^2 + c * H(t)
		val_k = obj(pk)
		val_k_p = obj(pk_p)

		tstar, val = 0., 0.
		if val_k > val_k_p
			tstar, val = pk, val_k
		else
			tstar, val = pk_p, val_k_p
		end

		#handle degenerate case where c == 0
		if c == 0.
			tp = -a/2/b		
			if pk <= tp <= pk_p
				valp = obj(tp)
				if valp >= val
					tstar, val = tp, valp
				end
			end
			return val, tstar
		end

		#else real work to be done
		inner = 2b/c * B/mu / M
    	inner *= exp( 2b/c * (1/M - 1) - a/c)

    	if inner < -1/MathConstants.e
    		tp = pk - 1
    	else
    		if inner < 0. ##two solutions.  Need to check -1 branch
				tp = 1- 1/M + c/2b * vopp.lambertw(inner)     		
				if pk <= tp <= pk_p
					valp = obj(tp)
					if valp >= val
						tstar, val = tp, valp
					end
				end
			end
			#always check principal branch
			tp= 1- 1/M + c/2b * vopp.lambertw0(inner) 
		end

		if pk <= tp <= pk_p
			valp = obj(tp)
			if valp >= val
				tstar, val = tp, valp
			end
		end

		return val, tstar
	end

	#Now kick it off!	   
	return _vopp_unimodal(Sc, mode_c, h, H, sep_fun, N=N, numCuts=numCuts, print_trace=print_trace)
end


function vopp_ub_unimodal_IC(mu, S, M, phat, q, mode; N=100, numCuts=500, print_trace=false)
	#standardize
	c = mu * (1 - M)
	Sc = vopp.comp_Sc(S, M)
	mode_c = (mode - c) / (mu - c)

	#create subfunctions
	v0 = (phat - 1)/M + 1
	h(t) = t >= v0 ? 1. - q : -q

	#H(t) = int_[t, mode_c] h
	#lazy implementation.
	function H(t)
		quadgk(h, t, 1, mode_c)[1]
	end

	#return maximizers of 
	#max_{[pk, pk_p]}  at + bt^2 + c H(t)
	function sep_fun(a, b, c, mode, pk, pk_p)
		obj(t) = a * t + b * t^2 + c * H(t)

		#check endpoints
		tstar, val = pk, obj(pk)
		if obj(pk_p) > val
			tstar = pk_p
			val = obj(pk_p)
		end

		#handle degenerate case where c == 0
		if c == 0.
			tp = -a/2/b		
			if pk <= tp <= pk_p
				if obj(tp) >= val
					tstar, val = tp, obj(tp)
				end
			end
			return val, tstar
		end

		#else real work to be done
		#check the two possible critical points
		if b != 0
			t1 = c/2/b * (1 - q - a/c)
			if t1 >= v0 && (pk <= t1 < pk_p)
				if obj(t1) >= val
					val = obj(t1)
					tstar = t1
				end
			end

			t0 = -c/2/b * (q + a/c)
			if t0 < v0 && (pk <= t0 < pk_p)
				if obj(t0) >= val
					val = obj(t0)
					tstar = t0
				end
			end
		end
		return val, tstar
	end  #end definition of separator

	#Now kick it off!	   
	return _vopp_unimodal(Sc, mode_c, h, H, sep_fun, N=N, numCuts=numCuts, print_trace=print_trace)
end



