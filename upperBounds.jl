###  
# Optimization problems for VoPP Upper Bounds
###

####Separation oracles
#hbar(t) = M * abs(t-1)/2 - D
#notice no dependence on mu.  
#maximum occurs at an end point or 1
function sep_MAD(a, b, lam1, lam2, M, D)
	lhs = lam1 * a + lam2/2 * abs(a-1) * M - lam2 * D
	rhs = lam1 * b + lam2/2 * abs(b-1) * M - lam2 * D
	mid = a <= 1. <= b ? lam1 - lam2 * D : -Inf
	val = max(lhs, rhs, mid)

	if lhs >= val
		return val, a
	elseif rhs >= val
		return val, b
	elseif mid >= val
		return val, 1.
	end
	throw("Unreachable")
	return -1., -1.
end

#hbar(t) = M^2 (t - 1)^2 - C^2
#Again, no dependence on mu
#maximum occurs at an end point or at unique minimizer
function sep_CV(a, b, lam1, lam2, M, C)
	lhs = lam1 * a + lam2 * M^2 * (a - 1)^2 - lam2 * C^2
	rhs = lam1 * b + lam2 * M^2 * (b - 1)^2 - lam2 * C^2
	vstar = 1 - lam1 / 2 / lam2 / M^2
	mid = a <= vstar <= b ? lam1 * vstar + lam2 * M^2 * (vstar - 1)^2 - lam2 * C^2 : -Inf
	val = max(lhs, rhs, mid)

	if lhs >= val
		return val, a
	elseif rhs >= val
		return val, b
	elseif mid >= val
		return val, vstar
	end
	throw("Unreachable")
	return -1., -1.	
end

#h(t) = -log(t/mu) + log (B/mu)   equivalent to Exp[log[V]]= log[B]
#hbar(t)  = -log(M(t -1 ) + 1) + log (B/mu)
function sep_GM(a, b, lam1, lam2, M, mu, B)
	lhs = lam1 * a - lam2 * log(M * (a - 1) + 1 ) + lam2 * log(B/mu)
	rhs = lam1 * b - lam2 * log(M * (b - 1) + 1 ) + lam2 * log(B/mu)

	vstar = lam2 <= 0 ? 1 - 1/M + lam2 / lam1 : a - 1e-6
	mid = a <= vstar <= b ? lam1 * vstar - lam2 * log(M * (vstar - 1) + 1) + lam2 * log(B/mu) : -Inf
	val = max(lhs, rhs, mid)
	
	if lhs >= val
		return val, a
	elseif rhs >= val
		return val, b
	elseif mid >= val
		return val, vstar
	end
	throw("Unreachable")
	return -1., -1.	
end

#### Functions based on a generic separation argument 
#h 1d function that represents moment constraint (includes moment):  E[h(V)] = 0
#sep_fun(a, b, lam1, lam2, M, mu) maximizes of lam1 * v  + lam2 h(mu * (Mv + 1-M) ).  returns optval, vstar
function _vopp_moment_dual(S, M, mu, h, sep_fun; N, numCuts=100, TOL =1e-6, print_trace=false)
	#translate everything to the Vc space
	hbar(t) = h( mu * (M * t + (1 - M)) )
	Sc = vopp.comp_Sc(S, M)

	@assert isfinite(Sc) "Currenty only supports finite Sc"

	p_grid = range(0, stop=Sc, length=N+1) #indexes from 1, doc indexes from 0
	m = Model(solver=GurobiSolver(OutputFlag = false))
	@variable(m, theta)
	@variable(m, lam[1:2])
	@variable(m, Qs[1:N+1] >= 0)  #indexes form 1, doc indexes from 0

	@constraint(m, sum(Qs) == 1)
	@constraint(m, theta + lam[1] * Sc + lam[2] * hbar(Sc) <= sum([p_grid[i] * Qs[i] for i = 1:N+1]))
	@objective(m, Max, theta + lam[1])
	
	for iter = 1:numCuts
		iter == numCuts && throw("Max Iteration Count Reached")
		status = solve(m)  #check status?
		feasible = true
		#check if any cut is violated
		for k = 2:N+1
			val, vstar = sep_fun(p_grid[k-1], p_grid[k], getvalue(lam[1]), getvalue(lam[2]), M, mu)
			lhs = val + getvalue(theta)
			rhs = sum([p_grid[j] * getvalue(Qs[j]) for j = 1:k-1])

			if lhs >= rhs + TOL
				@constraint(m, theta + vstar * lam[1] + hbar(vstar) * lam[2] <= sum( [p_grid[j] * Qs[j] for j = 1:k-1]) )
				feasible = false
				print_trace && println("iter $(iter): \t $lhs \t $rhs")
				continue
			end
		end
		if feasible
			@assert status == :Optimal "Cutting Procedure failed: $(status)"
			return 1/(getobjectivevalue(m) - TOL)
		end
	end
	##unreachable
	return -1.
end

######################
####
# Specialized function for MAD
# Not using separation.
# Algorithm relies on fact that minimizers are always at endpt or 1.    
####
function vopp_ub_MAD_opt(S, M, D; N = 500)
	#Solve the two cases for lambda2
	val_neg = _vopp_ub_MAD(S, M, D, false, N)
	val_pos = _vopp_ub_MAD(S, M, D, true, N)

	#take minimum bc we already took reciprocals
	min(val_neg, val_pos), val_neg, val_pos
end

function _vopp_ub_MAD(S, M, D, is_lam2_pos, N)
	Sc = comp_Sc(S, M)
	@assert isfinite(Sc) "Currenty only supports finite Sc"

	p_grid = range(0, stop=Sc, length=N+1) #indexes from 1, doc indexes from 0
	m = Model(solver=GurobiSolver(OutputFlag=false))
	@variable(m, theta)
	@variable(m, lam[1:2])
	@variable(m, Qs[1:N+1] >= 0)  #indexes form 1, doc indexes from 0

	if is_lam2_pos
		@constraint(m, lam[2] >= 0)
	else
		@constraint(m, lam[2] <= 0)
	end

	hbar(t) = M * abs(t-1)/2 - D

	@constraint(m, sum(Qs) == 1)
	@constraint(m, theta + lam[1] * Sc + lam[2] * hbar(Sc) <= sum([p_grid[i] * Qs[i] for i = 1:N+1]))
	@objective(m, Max, theta + lam[1])
	
	for k = 2:N+1
		@constraint(m, theta + p_grid[k-1] * lam[1] + hbar(p_grid[k-1]) * lam[2] <= sum( [p_grid[j] * Qs[j] for j = 1:k-1]) )
		@constraint(m, theta + p_grid[k]   * lam[1] + hbar(p_grid[k])   * lam[2] <= sum( [p_grid[j] * Qs[j] for j = 1:k-1]) )

		if !is_lam2_pos && p_grid[k-1] <= 1 <= p_grid[k]
			@constraint(m, theta + lam[1] + hbar(1.) * lam[2] <= sum( [p_grid[j] * Qs[j] for j = 1:k-1]) )			
		end
	end
	status = solve(m)
	@assert status == :Optimal "Problem Failed $status"
	return 1/(getobjectivevalue(m))
end

### Specialized function for CV
# useSep determines use separator or else Reformulation 
# Notice that the minimizer is not always at 1 because of linear offset
#Weirdly slow for some reason.  
function vopp_ub_CV(S, M, mu, C; method=:ConstraintGeneration, N=500, numCuts=100, TOL=1e-6, print_trace=false)
	if method == :Slemma
		return _vopp_CV_SLemma(S, M, mu, C, N)
	elseif method == :ConstraintGeneration
		h(v) = (v - mu)^2 / mu^2 - C^2
		sep_fun(a, b, lam1, lam2, M, mu) = sep_CV(a, b, lam1, lam2, M, C)
		return _vopp_moment_dual(S, M, mu, h, sep_fun, N=N, numCuts=numCuts, TOL=TOL, print_trace=print_trace)
	elseif method == :Reformulation
		return min( _vopp_ub_CV(S, M, mu, C, true, N), _vopp_ub_CV(S, M, mu, C, false, N) )
	else
		throw("Method must be one of :Slemma , :ConstraintGeneration , :Reformulation")
	end
end

function _vopp_ub_CV(S, M, mu, C, is_lam2_pos, N=500)
	#translate everything to the Vc space
	hbar(t) = M^2 * (t - 1)^2 - C^2  ##This functional form implicitly coded into dual constraints
	Sc = vopp.comp_Sc(S, M)
	@assert isfinite(Sc) "Currenty only supports finite Sc"

	p_grid = range(0, stop=Sc, length=N+1) #indexes from 1, doc indexes from 0
	m = Model(solver=GurobiSolver(OutputFlag = false))

	@variable(m, theta)
	@variable(m, lam[1:2])

	if is_lam2_pos
		@constraint(m, lam[2] >= 0)  #means solutions only occur at endpoints
	else
		@constraint(m, lam[2] <= 0)
	end

	@variable(m, Qs[1:N+1] >= 0)  #indexes form 1, doc indexes from 0
	@constraint(m, sum(Qs) == 1)

	@constraint(m, theta + lam[1] * Sc + lam[2] * hbar(Sc) <= sum([p_grid[i] * Qs[i] for i = 1:N+1]))

	#Add the auxiliary variables for each constraint in lam_2_neg
	if !is_lam2_pos
		@variable(m, alphas[2:N+1])
		@variable(m, betas[2:N+1])
		@variable(m, ts[2:N+1] >= 0)
	end

	for k = 2:N+1
		#always need the end points
		@constraint(m, theta + p_grid[k-1] * lam[1] + hbar(p_grid[k-1]) * lam[2] <= sum( [p_grid[j] * Qs[j] for j = 1:k-1]) )		
		@constraint(m, theta + p_grid[k] * lam[1] + hbar(p_grid[k]) * lam[2] <= sum( [p_grid[j] * Qs[j] for j = 1:k-1]) )		

		if !is_lam2_pos #concave case
			#use auxiliary variables to specify dual constraints
			@constraint(m, alphas[k] - lam[1] == betas[k])
			@constraint(m, betas[k]^2 <= -4 * M^2 * lam[2] * ts[k])
			@constraint(m, ts[k] - alphas[k] + lam[1] + alphas[k] * p_grid[k-1] <= 
								C^2 * lam[2] - theta + sum([p_grid[i] * Qs[i] for i = 1:k-1]))

			@constraint(m, ts[k] - alphas[k] + lam[1] + alphas[k] * p_grid[k] <= 
								C^2 * lam[2] - theta + sum([p_grid[i] * Qs[i] for i = 1:k-1]))
		end
	end
	@objective(m, Max, theta + lam[1])

	##Solve and extract solution 
	status = solve(m) 
	@assert status == :Optimal "Reformulation failed: $(status)"

	return 1/getobjectivevalue(m)
end


function _vopp_CV_SLemma(S, M, mu, C, N)
	#translate everything to the Vc space
	hbar(t) = M^2 * (t - 1)^2 - C^2  ##This functional form implicitly coded into dual constraints
	Sc = vopp.comp_Sc(S, M)
	@assert isfinite(Sc) "Currenty only supports finite Sc"

	ps = range(0, stop=Sc, length=N+1) #indexes from 1, doc indexes from 0
	m = Model(solver=GurobiSolver(OutputFlag = false))

	@variable(m, theta)
	@variable(m, lam[1:2])
	@variable(m, Qs[1:N+1] >= 0)  #doc indexes from 0 to N
	@constraint(m, sum(Qs) == 1)

	@constraint(m, theta + lam[1] * Sc + lam[2] * hbar(Sc) <= sum([ps[i] * Qs[i] for i = 1:N+1]))

	for k = 2:N+1
		#add auxiliary variables and specify modified S-Lemma
		x = @variable(m)
		y = @variable(m)
		z = @variable(m)

		@constraint(m, x >= 0)
		@constraint(m, y >= 0)
		@constraint(m, y >= -lam[2] * M^2)

		@constraint(m, ps[k-1] * ps[k] * (y + lam[2]*M^2) - x 
						>= 
						theta - sum(ps[j] * Qs[j] for j = 1:k) + 
						lam[2] * M^2 - lam[2] * C^2
					)

		@constraint(m, 4 * y * x >= z^2)
		@constraint(m, z == 2 * lam[2] * M^2 - lam[1] - 
							(ps[k-1] + ps[k]) * (y + lam[2] * M^2)
					)
	end
	@objective(m, Max, theta + lam[1])

	##Solve and extract solution 
	status = solve(m) 
	@assert status == :Optimal "Reformulation failed: $(status)"

	return 1/getobjectivevalue(m)
end






### Specialized function for BM
# useSep determines use separator or else Reformulation 
function vopp_ub_GM(S, M, mu, B; useSep=false, N=500, numCuts=100, TOL=1e-6, print_trace=false)
	if useSep 
		h(v) = -log(v/mu) + log(B/mu)
		sep_fun(a, b, lam1, lam2, M, mu) = sep_GM(a, b, lam1, lam2, M, mu, B)
		return _vopp_moment_dual(S, M, mu, h, sep_fun, N=N, numCuts=numCuts, TOL=TOL, print_trace=print_trace)
	end
	return min( _vopp_ub_GM(S, M, mu, B, true, N), _vopp_ub_GM(S, M, mu, B, false, N) )
end

### Specialized function for GM using convex conjugates
#uses a combination of IPOPT and prayer. 
#h(t) = -log(t/mu) + log(B/mu)
function _vopp_ub_GM(S, M, mu, B, is_lam2_pos, N)
	#translate everything to the Vc space
	hbar(t) = -log(M*t + 1- M) + log(B/mu)  ##This functional form implicitly coded into dual constraints

	Sc = vopp.comp_Sc(S, M)
	@assert isfinite(Sc) "Currenty only supports finite Sc"

	p_grid = range(0, stop=Sc, length=N+1) #indexes from 1, doc indexes from 0
	m = Model(solver=IpoptSolver())

	@variable(m, theta)
	@variable(m, lam[1:2])

	if is_lam2_pos
		@constraint(m, lam[2] >= 0)  #means solutions only occur at endpoints
	else
		@constraint(m, lam[2] <= 0)
	end

	@variable(m, Qs[1:N+1] >= 0)  #indexes form 1, doc indexes from 0
	@constraint(m, sum(Qs) == 1)

	@constraint(m, theta + lam[1] * Sc + lam[2] * hbar(Sc) <= sum([p_grid[i] * Qs[i] for i = 1:N+1]))

	#Add the auxiliary variables for each constraint in lam_2_neg
	if !is_lam2_pos
		@variable(m, rhos[2:N+1])
		@variable(m, ts[2:N+1])
	end

	for k = 2:N+1
		#always need the end points
		@constraint(m, theta + p_grid[k-1] * lam[1] + hbar(p_grid[k-1]) * lam[2] <= sum( [p_grid[j] * Qs[j] for j = 1:k-1]) )		
		@constraint(m, theta + p_grid[k]   * lam[1] + hbar(p_grid[k])   * lam[2] <= sum( [p_grid[j] * Qs[j] for j = 1:k-1]) )		

		if !is_lam2_pos #concave case
			#use auxiliary variables to specify dual constraints
			@constraint(m, rhos[k] >= lam[1])
			@NLconstraint(m, -lam[2] * log(-M * lam[2] / (rhos[k] - lam[1])) <= ts[k])
			@constraint(m, (1-M)/M * (rhos[k] - lam[1]) + lam[2] + lam[2]*log(B/mu) + ts[k] + 
						rhos[k] * p_grid[k-1] <= sum([p_grid[i] * Qs[i] for i = 1:k-1]) - theta)
			@constraint(m, (1-M)/M * (rhos[k] - lam[1]) + lam[2] + lam[2]*log(B/mu) + ts[k] + 
						rhos[k] * p_grid[k] <= sum([p_grid[i] * Qs[i] for i = 1:k-1]) - theta)

			# @constraint(m, -lam[2] * log( -rhos[k] / M /lam[2] ) <= ts[k])
			# @constraint(m, lam[1]* p_grid[k-1] - rhos[k] * p_grid[k-1] - lams[2] * log(B/mu) - lams[2] + 
			# 				(1 - 1/M) * rhos[k] + ts[k] <= sum([p_grid[i] * Qs[i] for i = 1:k-1]) - theta )
			# @constraint(m, lam[1]* p_grid[k] - rhos[k] * p_grid[k] - lams[2] * log(B/mu) - lams[2] + 
			# 				(1 - 1/M) * rhos[k] + ts[k] <= sum([p_grid[i] * Qs[i] for i = 1:k-1]) - theta )
		end
	end
	@objective(m, Max, theta + lam[1])

	##Solve and extract solution 
	status = solve(m) 
	@assert status == :Optimal "Reformulation failed: $(status)"

	return 1/getobjectivevalue(m)
end
