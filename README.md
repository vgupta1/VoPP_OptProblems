# VoPP_OptProblems

In the spirit of reproducible research, this repository contains all the code necessary for running the experiments and creating graphs in the paper:  
> "The Vaue of Personalized Pricing" by Adam Elmachtoub, Vishal Gupta and Michael Hamilton.

The full-text of the paper is available at the author's [website](http://faculty.marshall.usc.edu/Vishal-Gupta/research.html).

If you find this code or the paper useful, ***please consider citing it***.

## Overview
All of the source code for computing solutions by various methods can be found in src.jl, written in Julia. This code leverages the JuMP Modeling Language and the Gurobi Solver.   

The files:
 - TestMAD.jl
 - TestCV.jl

call the relevant workhorse functions to generate the principal plots of the paper.  

Finally the folder **plotting** contains functions used to generate plots for the paper, written in R and exported into Tikz.

## Licensing

This code is available under the MIT License.  
Copyright (c) 2019 Vishal Gupta
