####
#Debugging
####
#The following functions are primarily meant for debugging.  
#Not used in computations.  

#Uses a discretization in the primal space
#h is a 1D function defining moment constraint for the original problem 
#S, M, mu are the original parameters
function _vopp_ub_primal(S, M, mu, h; N=1000)
	#translate everything to the Vc space
	c = mu * (1 - M)
	hbar(t) = h((mu - c) * t + c)
	Sc = vopp.comp_Sc(S, M)

	#discretize the space
	@assert isfinite(Sc) "Currently only supports finite S"

	grid = range(0, stop=Sc, length=N)  #grid includes 0 and Sc
	m = Model(solver=GurobiSolver(OutputFlag=0))
	@variable(m, probs[1:N] >= 0)
	@constraint(m, sum(probs) == 1)  #normalization
	@constraint(m, sum([probs[i] * grid[i] for i = 1:N]) == 1)  #mean

	@constraint(m, sum( [hbar(grid[i]) * probs[i] for i = 1:N] ) == 0)  #moment
	@variable(m, y)

	for i = 1:N
		@constraint(m, y >= grid[i] * sum(probs[i:N]))  #pricing inequality
	end
	@objective(m, Min, y)
	solve(m)
	1/getvalue(y)
end

#Uses a discretization in the primal space of rectangular distributions
#S, M, mu, mode are the original parameters of NON-standardized dist
function _vopp_lb_primal_unimodal_MAD(mu, S, M, D, mode; N=1000, pj=-1)
	#convert to standardized distribution
	c = mu * (1 - M)
	Sc = vopp.comp_Sc(S, M)
	mode_c = (mode - c) / (mu - c)

	h(t) = M * abs(t - 1)/2 - D

	#H(t) = 1/abs(mode_c - t) * int_<t, mode_c> h if t != mode_c 
	function H(t)
		TOL = 1e-8
	    #Dirac case
	    if abs(t - mode_c) <= TOL
	        return h(mode_c)
	    end
	    
	    if 1 >= min(t, mode_c) && 1 <= max(t, mode_c)
	    	return quadgk(h, t, 1, mode_c)[1] / (mode_c - t)
	    end
	    
	    return quadgk(h, t, mode_c)[1]/ (mode_c - t)
	end

	function G(p, m, t)
		if abs(m - t) <= 1e-8
			return m >= p ? 1. : 0.
		elseif p > max(m, t)
			return 0.
		elseif p < min(m, t)
			return 1.
		end
		return (max(m, t) - p) / abs(m - t)
	end

	#ensure the dirac is in mixture
	t_grid = collect(range(0, stop=Sc, length=N))
	if !(mode_c in t_grid)
		push!(t_grid, mode_c)
		sort!(t_grid)
	end
	N = length(t_grid)
	
	#maximize pricing at p
	r_mp = -Inf
	best_price = -Inf
	dist = zero(t_grid)

	if pj < 0
		p_grid = t_grid[:] 
	else
		p_grid = [pj]
	end

	for p in p_grid
		m = Model(solver=GurobiSolver(OutputFlag=0))
		@variable(m, probMs[1:N] >= 0)
		@constraint(m, sum(probMs) == 1)  #normalization

		#mean
		@constraint(m, dot(t_grid .+ mode_c, probMs) == 2)
	
		#moment 
		@constraint(m, dot(map(H, t_grid), probMs) == 0)

		@objective(m, Max, p * dot(map(t-> G(p, mode_c, t), t_grid), 
									probMs)
					)

		status = solve(m)
		if status != :Optimal
			println("Status: $status for price $p  mode $mode_c")
		elseif getobjectivevalue(m) > r_mp 
			r_mp = getobjectivevalue(m)
			best_price = p
			dist[:] = getvalue(probMs)
		end
	end #end loop over prices
	return 1/r_mp, best_price, t_grid, dist
end

#NOTE:  Worst-Case distribution is a dirac at 0 (which is mode) and a uniform <0, S>
#As S goes to infinity, you can achieve a D of 2. 
#returns maximum coef of dev, mode and also the mixture distribution
#recall all mixtures written as <t, m> uniforms
#Inputs should be normalized
##VG Possibly Deprecate?
function max_D_unimodal(S; N=100)
	t_grid = range(0, stop=S, length=N)

	h(t) = abs(t-1)
	function H(t, m)
		TOL = 1e-8
		if abs(t-m) <= TOL
			return h(m)
		end
		if 1 <= min(t, m) || 1 >= max(t, m)
			return quadgk(s-> h(s), t, m)[1] / (m-t)
		else
			return quadgk(s-> h(s), t, 1, m)[1] / (m-t)
		end
		#unreachable
	end


	max_D = -Inf
	max_mode = -Inf
	dist = zero(t_grid)
	for mode in range(0, stop=2, length=N)
		m = Model(solver=GurobiSolver(OutputFlag=0))
		@variable(m, probMs[1:N] >= 0)
		@constraint(m, sum(probMs) == 1)  #normalization

		#mean
		@constraint(m, dot(t_grid .+ mode, probMs) == 2)
	
		#max coef of dev
		@objective(m, Max, .5 * dot(map(t->H(t, mode), t_grid), probMs))
		status = solve(m)
		if status != :Optimal
			println("Status: $status for mode $mode")
		elseif getobjectivevalue(m) > max_D 
			max_D = getobjectivevalue(m)
			max_mode = mode
			dist[:] = getvalue(probMs)
		end
	end #loop over mode

	return max_D, max_mode, dist, t_grid
end


#Uses a discretization in the primal space of rectangular distributions
#S, M, mu, mode are the original parameters of NON-standardized dist
function _vopp_lb_primal_unimodal_CV(mu, S, M, C, mode; N=1000, price=-1)
	#convert to standardized distribution
	c = mu * (1 - M)
	Sc = vopp.comp_Sc(S, M)
	mode_c = (mode - c) / (mu - c)

	h(t) = M^2 * (t-1)^2 - C^2

	#H(t) = 1/abs(mode_c - t) * int_<t, mode_c> h if t != mode_c 
	H(t) = abs(mode_c - t) <= 1e-8 ? h(mode_c) : quadgk(h, t, mode_c)[1]/ (mode_c - t)

	function G(p, m, t)
		if abs(m - t) <= 1e-8
			return m >= p ? 1. : 0.
		elseif p > max(m, t)
			return 0.
		elseif p < min(m, t)
			return 1.
		end
		return (max(m, t) - p) / abs(m - t)
	end

	#ensure the dirac is in mixture
	t_grid = collect(range(0, stop=Sc, length=N))
	if !(mode_c in t_grid)
		push!(t_grid, mode_c)
		sort!(t_grid)
	end
	N = length(t_grid)
	
	#maximize pricing at p
	r_mp = -Inf
	best_price = -Inf
	dist = zero(t_grid)

	if price > 0 
		p_grid = [price]
	else 
		p_grid = t_grid
	end

	for p in p_grid
		m = Model(solver=GurobiSolver(OutputFlag=0))
		@variable(m, probMs[1:N] >= 0)
		@constraint(m, sum(probMs) == 1)  #normalization

		#mean
		@constraint(m, dot(t_grid .+ mode_c, probMs) == 2)
	
		#moment 
		@constraint(m, dot(map(H, t_grid), probMs) == 0)

		@objective(m, Max, p * dot(map(t-> G(p, mode_c, t), t_grid), 
									probMs)
					)

		status = solve(m)
		if status != :Optimal
			println("Status: $status for price $p  mode $mode_c")
		elseif getobjectivevalue(m) > r_mp 
			r_mp = getobjectivevalue(m)
			best_price = p
			dist[:] = getvalue(probMs)
		end
	end #end loop over prices
	return 1/r_mp, best_price, dist
end

#conjecture that maximal CV occurs for uniform
function max_cv_unimodal_guess(mu, S, M, mode)
	#convert to standardized distribution
	c = mu * (1 - M)
	Sc = vopp.comp_Sc(S, M)
	mode_c = (mode - c) / (mu - c)

	if Sc < 2  #Unif[2-Sc, Sc]
		return Sc / sqrt(12)
	else
		return 2 / sqrt(12)
	end
end


