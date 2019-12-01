####
# Some Closed-Form Formulae and Helpers
####
#returns max coeff_dev for dist S, M
comp_Sc(S, M) = (S + M - 1)/M
max_dev(S, M) = M * (S - 1) / (S + M - 1)
function max_cv(S, M)
	Sc = comp_Sc(S, M)
	M * sqrt((Sc - 1))
end

#Returns the minimum B
function min_gm(S, M, mu) 
	Sc = comp_Sc(S, M)
	out = (Sc - 1)/Sc * log(1 - M) + 1/Sc * log(M * (S - 1) + 1)
	mu * exp(out)
end

#helper function for a robust (indefinite) quadratic  
#ax^2 + bx +c >= 0 for all x in [l, u]
#adds everything to model m
function addRQ!(m, a, b, c, l, u)
    y1 = @variable(m)  #theta
    @constraint(m, y1 >= 0)
    y = @variable(m)    #t
    @constraint(m, y >= 0)
    t2 = @variable(m)
    @constraint(m, 4 * y1 * y >= t2^2)
    @constraint(m, t2 == b - (l + u) * (y1 - a))
    @constraint(m, y1 - a >= 0)
    @constraint(m, c + l * u * (y1 - a) - y >= 0)
end


#Based on Theorem 3
function lower_bound_unimodal(D, M)
    if D/M <= 1/3
        return 1 / (1 - D/M)
    end

    return 8D / M / (1 + D/M )^2
end


#conjecture that maximal CV occurs for uniform distribution among unimodals
#Not formally proven
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