include("clusters.jl")
include("planetary_catalog.jl")
include("optimization.jl")
using Random

sim_param = setup_sim_param_model()





##### To start saving the model iterations in the optimization into a file:

model_name = "Hybrid1"
run_number = ARGS[1] # if want to run on the cluster as part of a job array: ARGS[1]
AD_mod = true
num_targs = 86760
max_evals = 10 #5000
dists_include = ["delta_f", "mult_CRPD_r", "periods_KS", "period_ratios_KS", "durations_KS", "duration_ratios_KS", "depths_KS", "radii_KS", "radius_ratios_KS", "radii_partitioning_KS", "radii_monotonicity_KS", "gap_complexity_KS"]
Pop_per_param = 4

file_name = model_name*"_targs$(num_targs)_evals$(max_evals)_run$(run_number).txt"
f = open(file_name, "w")
println(f, "# All initial parameters:")
write_model_params(f, sim_param)





##### To load the Kepler catalog:

ssk = calc_summary_stats_Kepler(stellar_catalog, planets_cleaned)

##### To load a file with the weights:

# To set the number of targets in the simulated catalog for each subsequent iteration:
add_param_fixed(sim_param,"num_targets_sim_pass_one", num_targs)

# To load and compute the weights, target distance, and target distance std from a precomputed file:
#active_param_true, weights, target_fitness, target_fitness_std = compute_weights_target_fitness_std_from_file("Maximum_AMD_model_split_stars_weights_ADmod_$(AD_mod)_targs86760_evals100_all_pairs.txt", 4950, sim_param; dists_include=dists_include, f=f)
active_param_true, weights, target_fitness, target_fitness_std = compute_weights_target_fitness_std_from_file_split_samples("Maximum_AMD_model_split_stars_weights_ADmod_$(AD_mod)_targs86760_evals100_all_pairs.txt", 4950, sim_param; names_samples=String[], dists_include_samples=[String[]], dists_include_all=dists_include, f=f)





##### To draw the initial values of the active parameters randomly within a search range:

active_param_keys = ["break1_mass", "log_rate_clusters", "log_rate_planets_per_cluster", "mean_ln_mass", "norm_radius", "power_law_P", "power_law_γ0", "power_law_γ1", "power_law_σ0", "power_law_σ1", "sigma_ln_mass", "sigma_logperiod_per_pl_in_cluster", "α_pret"]
active_params_box = [(1., 100.), (log(0.2), log(10.)), (log(0.2), log(10.)), (log(1.), log(100.)), (1., 5.), (-2., 2.), (0., 1.), (0., 1.), (0., 0.5), (0., 0.5), (log(1.), log(100.)), (0., 0.5), (1., 100.)] #search ranges for all of the active parameters

Random.seed!() # to have a random set of initial parameters and optimization run

active_param_start, active_param_transformed_start = draw_random_active_params(active_param_keys, active_params_box, sim_param)



PopSize = length(active_param_true)*Pop_per_param

println("# Active parameters: ", make_vector_of_active_param_keys(sim_param))
println(f, "# Active parameters: ", make_vector_of_active_param_keys(sim_param))
println(f, "# Starting active parameter values: ", active_param_start)
println(f, "# Optimization active parameters search bounds: ", active_params_box)
println(f, "# Method: adaptive_de_rand_1_bin_radiuslimited")
println(f, "# PopulationSize: ", PopSize)
println(f, "# AD_mod: ", AD_mod)
println(f, "# Distances used: ", dists_include)
println(f, "#")
println(f, "# Format: Active_params: [active parameter values]")
println(f, "# Format: Counts: [observed multiplicities][total planets, total planet pairs]")
println(f, "# Format: d_used_keys: [names of distance terms]")
println(f, "# Format: d_used_vals: [distance terms][sum of distance terms]")
println(f, "# Format: d_used_vals_w: [weighted distance terms][sum of weighted distance terms]")
println(f, "#")

target_function(active_param_start, sim_param; ss_fit=ssk, dists_include=dists_include, weights=weights["all"], AD_mod=AD_mod, f=f)




##### To use automated optimization routines to optimize the active model parameters:

# Pkg.add("BlackBoxOptim")       # only need to do these once
# Pkg.checkout("BlackBoxOptim")  # needed to get the lastest version
using BlackBoxOptim              # see https://github.com/robertfeldt/BlackBoxOptim.jl for documentation

t_elapsed = @elapsed begin
    opt_result = bboptimize(active_params -> target_function(active_params, sim_param; ss_fit=ssk, dists_include=dists_include, weights=weights["all"], AD_mod=AD_mod, f=f); SearchRange = active_params_box, NumDimensions = length(active_param_true), Method = :adaptive_de_rand_1_bin_radiuslimited, PopulationSize = PopSize, MaxFuncEvals = max_evals, TargetFitness = target_fitness, FitnessTolerance = target_fitness_std, TraceMode = :verbose)
end

println(f, "# best_candidate: ", best_candidate(opt_result))
println(f, "# best_fitness: ", best_fitness(opt_result))
println(f, "# elapsed time: ", t_elapsed, " seconds")
close(f)

