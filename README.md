# VoPP_OptProblems

In the spirit of reproducible research, this repository contains all the code necessary for running the experiments and creating graphs in the paper:  
> "The Vaue of Personalized Pricing" by Adam Elmachtoub, Vishal Gupta and Michael Hamilton.

The full-text of the paper is available at the author's [website](http://faculty.marshall.usc.edu/Vishal-Gupta/research.html).

If you find this code or the paper useful, ***please consider citing it***.

## Overview
All of the source code for computing solutions by various methods can be found in src.jl, written in Julia. This code leverages the JuMP Modeling Language and the Gurobi Solver.  The file src.jl includes the following subfiles with various routines
 - closedFormUB.jl  -- computes the closed-form upper bounds and tight distributions given coefficient of deviation 
 - lowerBoundsUnimodal -- mathematical optimization based lower bounds for unimodal distributions
 - UpperBoundsUnimodal -- mathematical optimization based upper bounds for unimodal distributions
 - LambertW.jl  -- A custom implementation of the Lambert W function (-1 branch)
 - helpers.jl -- small helper functions of various sorts.  
 - debugging.jl -- a collection of small scripts that do not yield proper bounds

The files 
 - lowerBoundsSymmetric.jl
 - upperBoundsSymmetric.jl
develop bounds for unimodal, symmetric distributions using mathematical optimization.  These bounds were not ultimately used in the final paper.  


The following are used to generate the data behind the principal plots of the paper 
 - simpleExperiments.jl -- generates the data for the various other plots (lambert W, lots of bounds with varying deviation, etc.)
 - shapeExperiments.jl  -- generates the upper and lower bounds via math programmming for coefficient of deviation and coefficient of variation

These files call workhorse functions from src. jl

Finally the R-Project **VOPP_OptProblems** contains various scripts used to read data from previous experiments and plot appropriately.  These scripts are written in R and exported to Tikz.   

## Licensing

This code is available under the MIT License.  
Copyright (c) 2019 Vishal Gupta
