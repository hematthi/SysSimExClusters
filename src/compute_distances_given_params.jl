include("clusters.jl")
include("planetary_catalog.jl")
include("optimization.jl")
using Random

sim_param = setup_sim_param_model()





##### To start saving the model iterations in the optimization into a file:

model_name = "Hybrid1"
AD_mod = true
num_targs = 86760
dists_include = ["delta_f", "mult_CRPD_r", "periods_KS", "depths_KS", "radii_KS", "radius_ratios_KS", "radii_partitioning_KS", "radii_monotonicity_KS"]
#dists_include = ["delta_f", "mult_CRPD_r", "periods_KS", "period_ratios_KS", "durations_KS", "duration_ratios_KS", "depths_KS", "radii_KS", "radius_ratios_KS", "radii_partitioning_KS", "radii_monotonicity_KS", "gap_complexity_KS"]

data_table = CSV.read("../emulator/GP_files/Active_params_distances_table_best10000_every1.txt", DataFrame)
params_keys = names(data_table)[1:end-1]
@assert all(make_vector_of_active_param_keys(sim_param) .== String.(params_keys))

params_array = Matrix(data_table[1:end, params_keys])

run_number, runs = 1, 1 #parse(Int64, ARGS[1]), parse(Int64, ARGS[2])
evals = Int(size(params_array,1)/runs)
start, stop = 1+(run_number-1)*evals, run_number*evals

file_name = model_name*"_recompute_optim_best10000_every1_evals$(start)to$(stop)_targs$(num_targs).txt"
f = open(file_name, "w")
println(f, "# All initial parameters:")
write_model_params(f, sim_param)





##### To load the Kepler catalog:

ssk = calc_summary_stats_Kepler(stellar_catalog, planets_cleaned)

##### To load a file with the weights:

# To set the number of targets in the simulated catalog for each subsequent iteration:
add_param_fixed(sim_param,"num_targets_sim_pass_one", num_targs)

# To load and compute the weights, target distance, and target distance std from a precomputed file:
active_param_true, weights, target_fitness, target_fitness_std = compute_weights_target_fitness_std_from_file_split_samples("Maximum_AMD_model_split_stars_weights_ADmod_$(AD_mod)_targs86760_evals100_all_pairs.txt", 4950, sim_param; names_samples=String[], dists_include_samples=[String[]], dists_include_all=dists_include, f=f)





##### To recompute the model with the parameters in the table:

println(f, "# Active parameters: ", String.(params_keys))
println(f, "# AD_mod: ", AD_mod)
println(f, "# Distances used: ", dists_include)
println(f, "#")
println(f, "# Format: Active_params: [active parameter values]")
println(f, "# Format: Counts: [observed multiplicities][total planets, total planet pairs]")
println(f, "# Format: d_used_keys: [names of distance terms]")
println(f, "# Format: d_used_vals: [distance terms][sum of distance terms]")
println(f, "# Format: d_used_vals_w: [weighted distance terms][sum of weighted distance terms]")
println(f, "#")

Random.seed!()

t_elapsed = @elapsed begin
    for i in start:stop #1:size(params_array,1)
        target_function(params_array[i,:], sim_param; ss_fit=ssk, dists_include=dists_include, weights=weights["all"], AD_mod=AD_mod, f=f)
    end
end

println(f, "# elapsed time: ", t_elapsed, " seconds")
close(f)
