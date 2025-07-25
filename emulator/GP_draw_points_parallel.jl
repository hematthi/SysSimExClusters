using Distributed

addprocs(4, exeflags="--project=..") #number of additional processors

using SharedArrays
#@everywhere using DistributedArrays
#@everywhere using ParallelDataTransfer

@everywhere include("GP_functions.jl")
include("GP_emulator.jl")





"""
Evaluates 'predict_model_from_uniform_prior_until_accept_n_points_parallel' given a dataset and GP model (i.e. kernel and hyperparameters) to draw 'n_points' with GP mean and std better than 'max_mean' and 'max_std' and save the results as a table in a CSV file.
"""
function draw_points_parallel_with_GP_and_save(n_points::Int64; params_names::Array{String,1}, xdata::Array{Float64,2}, mean_f::Float64, ydata::Vector{Float64}, ydata_err::Vector{Float64}, kernel::Function, hparams::Vector{Float64}, prior_bounds::Union{Array{Tuple{Float64,Float64},1}, Nothing}=nothing, max_mean::Float64=Inf, max_std::Float64=Inf, max_post::Float64=Inf, save_path::String="", run_number::Int64=1)
    @assert size(xdata, 2) == length(params_names)
    @assert size(xdata, 1) == length(ydata) == length(ydata_err)
    @assert n_points > 0
    dims = length(params_names)
    n_train = length(ydata)

    if prior_bounds == nothing
        xtrain_mins, xtrain_maxs = findmin(xdata, dims=1)[1], findmax(xdata, dims=1)[1]
        prior_bounds = [(xtrain_mins[i], xtrain_maxs[i]) for i in 1:dims] # To prevent predicting "outside" of the n-dim box of training points
    end
    vol = prod([x[2]-x[1] for x in prior_bounds])

    println("Drawing $n_points points satisfying GP mean < $max_mean, GP std < $max_std, GP posterior draw < $max_post...")
    @time prior_draws_GP_table = predict_model_from_uniform_prior_until_accept_n_points_parallel(params_names, xdata, ydata, kernel, hparams, ydata_err, n_points; prior_bounds=prior_bounds, max_mean=max_mean, max_std=max_std, max_post=max_post)

    file_name = joinpath(save_path, "GP_train$(n_train)_meanf$(mean_f)_sigmaf$(hparams[1])_lscales$(round(sum(hparams[2:end]), digits=2))_vol$(round(vol, digits=2))_points$(n_points)_mean$(max_mean)_std$(max_std)_post$(max_post)_run$(run_number).csv")
    f = open(file_name, "w")
    println(f, "# n_data = $n_train, mean_f = $mean_f")
    println(f, "# hparams = $hparams")
    println(f, "# prior_box = $prior_bounds, box_vol = $vol")
    println(f, "# n_accept = $n_points, max_mean = $max_mean, max_std = $max_std, max_post = $max_post")
    CSV.write(f, prior_draws_GP_table; append=true, writeheader=true)
    close(f)
end





# To use the optimized GP model to predict at a large number of points drawn from the prior:

Random.seed!()

n_points = 10000
max_mean, max_std, max_post = Inf, Inf, -12.
save_path = data_path
draw_points_parallel_with_GP_and_save(n_points; params_names=GP_model[:params_names], xdata=GP_model[:xtrain], mean_f=GP_model[:mean_f], ydata=GP_model[:ytrain], ydata_err=GP_model[:ytrain_err], kernel=kernel_SE_ndims, hparams=GP_model[:hparams_best], prior_bounds=prior_bounds, max_mean=max_mean, max_std=max_std, max_post=max_post, save_path=save_path) # run_number=parse(Int64, ARGS[1])
