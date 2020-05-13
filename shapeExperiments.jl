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
		#Dumb bound for plotting ease.
		writedlm(f, [dev "None" true 1.0 0.])

		t = @elapsed bound = vopp.vopp_ub_MAD(S, M, dev)[1]
		writedlm(f, [dev "None" false bound t])

		t = @elapsed bound = vopp.vopp_ub_unimodal_MAD(mu, S, M, dev, mode)
		writedlm(f, [dev "Unimodal" false bound t])

		# t = @elapsed bound = vopp.vopp_ub_symmetric_MAD(mu, S, M, dev, mode)
		# writedlm(f, [dev "Symmetric" false bound t])

		t = @elapsed bound = vopp.vopp_lb_unimodal_MAD(mu, S, M, dev, mode, method=:Formula)
		writedlm(f, [dev "UniAnalytic" true bound t])

		t = @elapsed bound = vopp.vopp_lb_unimodal_MAD(mu, S, M, dev, mode, safe_fail=false)[1]
		writedlm(f, [dev "Unimodal" true bound t])

		# t = @elapsed bound = vopp.vopp_lb_symmetric_MAD(mu, S, M, dev, safe_fail=false)[1]
		# writedlm(f, [dev "Symmetric" true bound t])

		flush(f)
	end
	close(f)
end


function runCVExp(file_out_name;
	S=2, M=1, mu=1, mode=1, 
	cv_grid = range(0., stop=.5773502691896258, length=100))

	#Create file_out and header
	outPath = "$(file_out_name)__$(S)_$(M)_$(mu)_$(mode)"
	f = open("$(outPath).tab", "w")
	writedlm(f, ["CV" "Shape" "isLB" "Bound" "Time"])

	for C in cv_grid
		#add a dumb bound for ease in plotting
		writedlm(f, [C "None" true 1.0 0.])

		t = @elapsed bound = vopp.vopp_ub_CV(S, M, mu, C)
		writedlm(f, [C "None" false bound t])

		t = @elapsed bound = vopp.vopp_ub_unimodal_CV(mu, S, M, C, mode)
		writedlm(f, [C "Unimodal" false bound t])

		# t = @elapsed bound = vopp.vopp_ub_symmetric_CV(mu, S, M, C, mode)
		# writedlm(f, [C "Symmetric" false bound t])

		t = @elapsed bound = vopp.vopp_lb_unimodal_CV(mu, S, M, C, mode)[1]
		writedlm(f, [C "Unimodal" true bound t])

		# t = @elapsed bound = vopp.vopp_lb_symmetric_CV(mu, S, M, C, safe_fail=true)[1]
		# writedlm(f, [C "Symmetric" true bound t])

		flush(f)
	end
	close(f)
end


function runGMExp(file_out_name;
	S=2, M=.9, mu=1, mode=1, 
	B_grid = range(.81168, stop=1., length=100))

	#Create file_out and header
	outPath = "$(file_out_name)__$(S)_$(M)_$(mu)_$(mode)"
	f = open("$(outPath).tab", "w")
	writedlm(f, ["GM" "Shape" "isLB" "Bound" "Time"])

	for B in B_grid
		#add a dumb bound for ease in plotting
		writedlm(f, [B "None" true 1.0 0.])

		t = @elapsed bound = vopp.vopp_ub_GM(S, M, mu, B, useSep=true)							
		writedlm(f, [B "None" false bound t])

		t = @elapsed bound = vopp.vopp_ub_unimodal_GM(mu, S, M, B, mode)
		writedlm(f, [B "Unimodal" false bound t])

		t = @elapsed bound = vopp.vopp_lb_unimodal_GM(mu, S, M, B, mode)[1]
		writedlm(f, [B "Unimodal" true bound t])

		flush(f)
	end
	close(f)
end

function runICExp(file_out_name;
	S=2, M=.9, mu=1., mode=1., phat = .8)
	qmin, qmax = vopp.min_max_q_IC_uni(S, M, mu, mode, phat, 500)
	q_grid = range(qmin, stop=qmax, length=100)

	#Create file_out and header
	outPath = "$(file_out_name)__$(S)_$(M)_$(mu)_$(mode)_$(phat)"
	f = open("$(outPath).tab", "w")
	writedlm(f, ["IC" "Shape" "isLB" "Bound" "Time"])

	for q in q_grid
		#add a dumb bound for ease in plotting
		writedlm(f, [q "None" true 1.0 0.])

		t = @elapsed bound = vopp.vopp_ub_IC(S, M, mu, phat, q)
		writedlm(f, [q "None" false bound t])

		t = @elapsed bound = vopp.vopp_ub_unimodal_IC(mu, S, M, phat, q, mode)
		writedlm(f, [q "Unimodal" false bound t])

		t = @elapsed bound = vopp.vopp_lb_unimodal_IC(mu, S, M, phat, q, mode)[1]
		writedlm(f, [q "Unimodal" true bound t])

		flush(f)
	end
	close(f)
end


if ARGS[1] == "MAD"
	runMADExp(ARGS[2])
elseif ARGS[1] == "MAD2"
	runMADExp(ARGS[2], M=.9, dev_grid = range(0., stop=.2368245, length=100))
elseif ARGS[1]== "CV"
	runCVExp(ARGS[2])
elseif ARGS[1] == "CV2"
	runCVExp(ARGS[2], M=.9, cv_grid = range(0., stop=.5477125, length=100))
elseif ARGS[1] == "GM"
	runGMExp(ARGS[2])
elseif ARGS[1] == "IC"
	runICExp(ARGS[2])
end

