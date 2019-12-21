##Simple Plots of MAD UB bound
using DelimitedFiles
include("src.jl")  #load up vopp

function simplePlot(file_out_name; 
					S=4, M=1, N=100)
	#Create file_out and header
	outPath = "$(file_out_name)__$(S)_$(M)_$(N)"
	f = open("$(outPath).tab", "w")
	writedlm(f, ["D" "Bound"])

	for  D in range(1e-10, stop=vopp.delta_h(S) - 1e-10, length=N)
		#add a dumb bound for ease in plotting
		bound = vopp.vopp_ub_MAD(S, 1., D)	
		writedlm(f, [D bound])
	end
	close(f)
end

function contour_plot(file_out_name;
				Smin = 1.1, Smax=4, N=100, TOL = 1e-10)
	#Create file_out and header
	outPath = "$(file_out_name)__$(Smin)_$(Smax)_$(N)"
	f = open("$(outPath).tab", "w")
	writedlm(f, ["S" "D" "Bound"])

	S_grid = range(Smin, stop=Smax, length=N)
	D_grid = range(1e-10, stop=1, length=N)

	for S in S_grid
		for D in D_grid
			if D > vopp.delta_h(S) - TOL
				continue
			end
			bound = vopp.vopp_ub_MAD(S, 1., D)
			writedlm(f, [S D bound])
		end
		flush(f)
	end
	close(f)
end

function scale_bound_dist(file_out_name; S=4, N=100)
	#Create file_out and header
	outPath = "$(file_out_name)__$(S)_$(N)"
	f = open("$(outPath).tab", "w")
	writedlm(f, ["x" "Fbar"])

	x_grid = range(0, stop=S, length=N)
	alpha = 1/vopp_ub_scale(S, 1)
	push!(x_grid, alpha)
	sort!(x_grid)

	for x in x_grid
		writedlm(f, [x vopp.tight_dist_ub_scale(x, S)])
	end
	close(f)
end

if ARGS[1] == "Simple"
	simplePlot(ARGS[2])
elseif ARGS[1] == "Contour"
	contour_plot(ARGS[2])
elseif ARGS[1] == "scale_tight_dist"
	scale_bound_dist(ARGS[2])
end
