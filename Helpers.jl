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
function _addRQ!(m, a, b, c, l, u)
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

#conjecture that maximal CV occurs for uniform distribution among unimodals
#Not formally proven
function max_std_cv_unimodal_guess(mu, S, M, isSym)
    #convert to standardized distribution
    c = mu * (1 - M)
    Sc = vopp.comp_Sc(S, M)

    #adjust for symmetry if needed
    if isSym 
        if Sc > 2 
            Sc = 2
        end
        if Sc < 2 
            lb = 2 - Sc
        else
            lb = 0
        end
    end

    return 2 / sqrt(12) * (Sc - lb) / (Sc + lb)
end


#for the standardized distribution
#document indexes 0 to N
function geom_price_ladder(S, delta)
    N = ceil(Int, 1 + log(S/delta) / log(1 + delta))
    ps = zeros(N + 1)
    ps[0 + 1] = 0.
    ps[1 + 1] = delta
    ps[N + 1] = S
    for i = 2:N-1
        ps[i + 1] = ps[i - 1 + 1] * (1 + delta)
    end
    ps
end


vopp_ub_scale(S, M) = -lambertw(-M / exp(1) / (S + M - 1))

#assumes standardized
function tight_dist_ub_scale(x, S)
    alpha = 1/vopp_ub_scale(S, 1)
    if 0 <= x <= alpha
        return 1.
    elseif x <= S
        return alpha / x
    else
        return 0.
    end
    return -1.
end




