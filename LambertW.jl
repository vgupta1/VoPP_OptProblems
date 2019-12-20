#Lambert-W

## Lambert-W package on github does not work with julia 0.7
## Native implementations
function lambertw_lb(z)
	@assert -exp(-1) <= z <= 0 "Lambertw only defined on [-1/e, 0] : $z"
	x = -exp(1) * z
	-1 - sqrt(2 * log(1/x)) - log(1/x)
end

function lambertw_ub(z)
	@assert -exp(-1) <= z <= 0 "Lambertw only defined on [-1/e, 0] : $z"
	x = -exp(1) * z
	-1 - sqrt(2 * log(1/x)) - 2 * log(1/x) / 3
end

function lambertw(x)
	@assert -exp(-1) <= x <= 0 "Lambertw only defined on [-1/e, 0] : $x"
	#use the Chaterizighou bound to get a guesstimate
	lb = lambertw_lb(x)
	f(z) = z * exp(z) - x
	find_zero(f, (lb, -1.)) 
end
