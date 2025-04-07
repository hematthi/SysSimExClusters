#Pkg.add("ParallelDataTransfer")
using Distributed

addprocs(4, exeflags="--project=$(Base.active_project())") # add number of additional processors

@everywhere include("optimize.jl")

run_ids = 31:50
max_evals = 5000

#@distributed for i=run_ids
#    println("Worker $(myid()) is running i=$i")
#    setup_and_run_optimizer(i)
#end

t_runs = pmap(i->setup_and_run_optimizer(i, max_evals=max_evals), run_ids)
