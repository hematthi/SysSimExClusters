include("clusters.jl")
include("planetary_catalog.jl")
include("optimization.jl")
using Random

# Pkg.add("BlackBoxOptim")       # only need to do these once
# Pkg.checkout("BlackBoxOptim")  # needed to get the lastest version
using BlackBoxOptim              # see https://github.com/robertfeldt/BlackBoxOptim.jl for documentation





function setup_and_run_optimizer(run_number::Int64=parse(Int64,ARGS[1]); max_evals::Int64=10)

    ##### Settings:
    
    model_name = "Hybrid1"
    AD_mod = true
    num_targs = 86760
    dists_include = ["delta_f", "mult_CRPD_r", "periods_KS", "depths_KS", "radii_KS", "radius_ratios_KS", "radii_partitioning_KS", "radii_monotonicity_KS"]
    pop_per_param = 4
    
    active_param_keys = ["log_rate_clusters", "log_rate_planets_per_cluster", "log_α_pret", "mean_ln_mass", "norm_radius", "power_law_P", "power_law_γ0", "power_law_σ0", "sigma_ln_mass"]
    active_params_box = [(log(0.2), log(5.)), (log(0.2), log(5.)), (log(0.01), log(1000.)), (log(1.), log(100.)), (1., 5.), (-2., 2.), (0., 1.), (0., 0.5), (log(1.), log(100.))] # search ranges for all of the active parameters
    
    #####
    
    
    
    ##### Set up the simulation parameters:
    
    sim_param = setup_sim_param_model()
    add_param_fixed(sim_param, "num_targets_sim_pass_one", num_targs)



    ##### Open a file for saving the model iterations in the optimization run:

    file_name = model_name*"_targs$(num_targs)_evals$(max_evals)_run$(run_number).txt"
    f = open(file_name, "w")
    println(f, "# All initial parameters:")
    write_model_params(f, sim_param)



    ##### To load the Kepler catalog:

    ssk = calc_summary_stats_Kepler(stellar_catalog, planets_cleaned)



    ##### Load and compute the weights, target distance, and target distance std from a precomputed file:
    #active_param_true, weights, target_fitness, target_fitness_std = compute_weights_target_fitness_std_from_file("Maximum_AMD_model_split_stars_weights_ADmod_$(AD_mod)_targs86760_evals100_all_pairs.txt", 4950, sim_param; dists_include=dists_include, f=f)
    active_param_true, weights, target_fitness, target_fitness_std = compute_weights_target_fitness_std_from_file_split_samples("Maximum_AMD_model_split_stars_weights_ADmod_$(AD_mod)_targs86760_evals100_all_pairs.txt", 4950, sim_param; names_samples=String[], dists_include_samples=[String[]], dists_include_all=dists_include, f=f)

    pop_size = length(active_param_true)*pop_per_param
    
    
    
    ##### Draw the initial values of the active parameters randomly within the search box:
    
    Random.seed!() # to have a random set of initial parameters and optimization run

    active_param_start, active_param_transformed_start = draw_random_active_params(active_param_keys, active_params_box, sim_param)

    println("# Active parameters: ", make_vector_of_active_param_keys(sim_param))
    println(f, "# Active parameters: ", make_vector_of_active_param_keys(sim_param))
    println(f, "# Starting active parameter values: ", active_param_start)
    println(f, "# Optimization active parameters search bounds: ", active_params_box)
    println(f, "# Method: adaptive_de_rand_1_bin_radiuslimited")
    println(f, "# PopulationSize: ", pop_size)
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



    ##### Run the DE optimizer:

    t_elapsed = @elapsed begin
        opt_result = bboptimize(active_params -> target_function(active_params, sim_param; ss_fit=ssk, dists_include=dists_include, weights=weights["all"], AD_mod=AD_mod, f=f); SearchRange = active_params_box, NumDimensions = length(active_param_true), Method = :adaptive_de_rand_1_bin_radiuslimited, PopulationSize = pop_size, MaxFuncEvals = max_evals, TargetFitness = target_fitness, FitnessTolerance = target_fitness_std, TraceMode = :verbose)
    end

    println(f, "# best_candidate: ", best_candidate(opt_result))
    println(f, "# best_fitness: ", best_fitness(opt_result))
    println(f, "# elapsed time: ", t_elapsed, " seconds")
    close(f)
    
    return t_elapsed
end





##### To run the optimizer:

#t_elapsed = setup_and_run_optimizer(0; max_evals=10)
