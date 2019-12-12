####
# Comparing Bounds With Various Shape Constraints
# ARGS[1] i either MAD or CV
# ARGS[2] is outPathStub
####

using DelimitedFiles
include("src.jl")  #load up vopp

#the default values were used to generate data for experiment
function runMADExp(file_out_name;
	S=2, M=1, mu=1, mode=1, 
	dev_grid = range(0., stop=.25, length=100))

	#Create file_out and header
	outPath = "$(file_out_name)__$(S)_$(M)_$(mu)_$(mode)"
	f = open("$(outPath).tab", "w")
	writedlm(f, ["Dev" "Shape" "isLB" "Bound" "Time"])

	for dev in dev_grid
		t = @elapsed bound = vopp.vopp_MAD(S, M, dev)[1]
		writedlm(f, [dev "None" false bound t])

		t = @elapsed bound = vopp.vopp_MAD_unimodal(mu, S, M, dev, mode)
		writedlm(f, [dev "Unimodal" false bound t])

		t = @elapsed bound = vopp.vopp_sym_unimodal_MAD(mu, S, M, dev, mode)
		writedlm(f, [dev "Symmetric" false bound t])

		t = @elapsed bound = vopp.lower_bound_unimodal(dev, M)
		writedlm(f, [dev "UniAnalytic" true bound t])

		t = @elapsed bound = vopp.vopp_lb_unimodal_MAD(mu, S, M, dev, mode, safe_fail=true)[1]
		writedlm(f, [dev "Unimodal" true bound t])

		t = @elapsed bound = vopp.vopp_lb_symmetric_MAD(mu, S, M, dev, safe_fail=true)[1]
		writedlm(f, [dev "Symmetric" true bound t])

		flush(f)
	end
	close(f)
end

if ARGS[1] == "MAD"
	runMADExp(ARGS[2])
elseif ARGS[1]== "CV"
	throw("CV Not yet implemented")
end

