dir_path = dirname(@__FILE__)

include(joinpath(dir_path, "../src/clusters.jl"))
include(joinpath(dir_path, "../src/planetary_catalog.jl"))
include(joinpath(dir_path, "../src/optimization.jl"))

##### To load model parameters found using the GP emulator and simulate catalogs if they pass a distance threshold:

# 12 params:
#save_path = "/Users/hematthi/Documents/GradSchool/Research/SysSim/Simulated_catalogs/Hybrid_NR20_AMD_model1/Fit_all_KS/Params12/GP_best_models_100"

#GP_data_path = "/Users/hematthi/Documents/NotreDame_Postdoc/CRC/Files/SysSim/Model_Optimization/Hybrid_NR20_AMD_model1/Fit_all_KS/Params12/GP_files"
#GP_file_name = "GP_train2000_meanf35.0_sigmaf2.7_lscales37.65_vol1425.6_points100000_meanInf_stdInf_post-10.0.csv"

# 8 params:
#save_path = "/Users/hematthi/Documents/GradSchool/Research/SysSim/Simulated_catalogs/Hybrid_NR20_AMD_model1/Fit_all_KS/Params8/GP_best_models_100"

#GP_data_path = "/Users/hematthi/Documents/NPP_ARC_Modernize_Kepler/Personal_research/SysSim/Model_Optimization/Hybrid_NR20_AMD_model1/Fit_all_KS/Params8/GP_files"
#GP_file_name = "GP_train2000_meanf35.0_sigmaf2.7_lscales16.05_vol14.26_points100000_meanInf_stdInf_post-10.0.csv"

# Fit some KS, 8 params:
#save_path = "/Users/hematthi/Documents/GradSchool/Research/SysSim/Simulated_catalogs/Hybrid_NR20_AMD_model1/Fit_some_KS/Params8_fix_highM/GP_best_models_100"

#GP_data_path = "/Users/hematthi/Documents/NPP_ARC_Modernize_Kepler/Personal_research/SysSim/Model_Optimization/Hybrid_NR20_AMD_model1/Fit_some_KS/Params8_fix_highM/GP_files"
#GP_file_name = "GP_train2000_meanf25.0_sigmaf2.7_lscales2.45_vol112.9_points100000_meanInf_stdInf_post-13.0.csv"

# Fit some KS, 9 params:
#save_path = "/Users/hematthi/Documents/GradSchool/Research/SysSim/Simulated_catalogs/Hybrid_NR20_AMD_model1/Fit_some8_KS/Params9_fix_highM/GP_best_models_100"

#GP_data_path = "/Users/hematthi/Documents/NPP_ARC_Modernize_Kepler/Personal_research/SysSim/Model_Optimization/Hybrid_NR20_AMD_model1/Fit_some8_KS/Params9_fix_highM/GP_files"
#GP_file_name = "GP_train2000_meanf30.0_sigmaf2.7_lscales3.03_vol240.84_points10000_meanInf_stdInf_post-13.0.csv"

# Fit some+1 KS, 9 params:
save_path = "/Users/hematthi/Documents/GradSchool/Research/SysSim/Simulated_catalogs/Hybrid_NR20_AMD_model1/Fit_some8p1_KS/Params9_fix_highM/GP_best_models_100"

GP_data_path = "/Users/hematthi/Documents/NPP_ARC_Modernize_Kepler/Personal_research/SysSim/Model_Optimization/Hybrid_NR20_AMD_model1/Fit_some8p1_KS/Params9_fix_highM/GP_files"
GP_file_name = "GP_train2000_meanf30.0_sigmaf2.7_lscales2.22_vol48.52_points10000_meanInf_stdInf_post-12.0.csv"

GP_points = CSV.read(joinpath(GP_data_path, GP_file_name), DataFrame, comment="#")
active_params_names = names(GP_points)[1:end-3]
active_params_best_all = GP_points[!,active_params_names]

# If transformed:
#=
i_tf, j_tf = 2,3 # indices of transformed parameters
active_params_names[i_tf:j_tf] = ["log_rate_clusters", "log_rate_planets_per_cluster"]
active_params_best_all[!,i_tf], active_params_best_all[!,j_tf] = (active_params_best_all[!,i_tf] .- active_params_best_all[!,j_tf])/2., (active_params_best_all[!,i_tf] .+ active_params_best_all[!,j_tf])/2.
rename!(active_params_best_all, Symbol.(active_params_names))
=#

model_name = "Hybrid1"
AD_mod = true
num_targs = 86760
dists_include = ["delta_f", "mult_CRPD_r", "depths_KS", "radii_KS", "radius_ratios_KS", "radii_partitioning_KS", "radii_monotonicity_KS"]
#dists_include = ["delta_f", "mult_CRPD_r", "periods_KS", "period_ratios_KS", "durations_KS", "duration_ratios_KS", "depths_KS", "radii_KS", "radius_ratios_KS", "radii_partitioning_KS", "radii_monotonicity_KS", "gap_complexity_KS"]

d_threshold, mean_f = 18., 30.
n_pass = 100 # number of simulations we want to pass the distance threshold
n_save = 100 # number of simulations we want to pass the distance threshold and also save (choose a small number or else requires a lot of storage space); must not be greater than n_pass!

file_name = model_name*"_pass_GP_meanf$(mean_f)_thres$(d_threshold)_pass$(n_pass)_targs$(num_targs).txt"
f = open(joinpath(save_path, file_name), "w")





##### To load the Kepler catalog:

ssk = calc_summary_stats_Kepler(stellar_catalog, planets_cleaned; sim_param=sim_param)





##### To load a file with the weights:

active_param_true, weights, target_fitness, target_fitness_std = compute_weights_target_fitness_std_from_file_split_samples("../src/Maximum_AMD_model_split_stars_weights_ADmod_$(AD_mod)_targs86760_evals100_all_pairs.txt", 4950, sim_param; names_samples=String[], dists_include_samples=[String[]], dists_include_all=dists_include, f=f)

##### Also include the radii difference from the radius valley:
# NOTE: this is included ad hoc here since the weights file does not have a weight for this distance term yet; use the same weight as the radii distribution for now

push!(dists_include, "radii_delta_valley_KS")
weights["all"]["radii_delta_valley_KS"] = weights["all"]["radii_KS"]
#####





##### To simulate a catalog with each set of params from the GP emulator that passed the threshold, saving all the distances, and saving the catalog if the true distance also passes the threshold:

sim_param = setup_sim_param_model()
add_param_fixed(sim_param,"num_targets_sim_pass_one", num_targs)

println(f, "# Active parameters: ", String.(active_params_names))
println(f, "# AD_mod: ", AD_mod)
println(f, "# Distances used: ", dists_include)
println(f, "#")
println(f, "# Format: Active_params: [active parameter values]")
println(f, "# Format: Counts: [observed multiplicities][total planets, total planet pairs]")
println(f, "# Format: d_used_keys: [names of distance terms]")
println(f, "# Format: d_used_vals: [distance terms][sum of distance terms]")
println(f, "# Format: d_used_vals_w: [weighted distance terms][sum of weighted distance terms]")
println(f, "#")

sim_count = 0
pass_count = 0
save_count = 0
summary_array = Array{Float64,2}(undef, 0, size(GP_points,2)+1)

t_elapsed = @elapsed begin
    while pass_count < n_pass && sim_count < size(active_params_best_all,1)
        global sim_count, pass_count, save_count, summary_array
        sim_count += 1
        println("# Generating simulated catalog ", sim_count)

        # To set up the model parameters:
        for (i,param_name) in enumerate(active_params_names)
            add_param_active(sim_param, string(param_name), active_params_best_all[sim_count, param_name])
        end
        println(f, "Active_params: ", collect(active_params_best_all[sim_count,:])) # to write the params to file

        # Generate a simulated catalog:
        cat_phys = generate_kepler_physical_catalog(sim_param)
        cat_phys_copy = deepcopy(cat_phys) # need to deepcopy to save later, since cat_phys_cut overwrites cat_phys too
        cat_phys_cut = ExoplanetsSysSim.generate_obs_targets(cat_phys,sim_param)
        cat_obs = observe_kepler_targets_single_obs(cat_phys_cut,sim_param)

        # Compute the summary statistics:
        summary_stat = calc_summary_stats_model(cat_obs,sim_param)

        # Compute the individual and total weighted distances:
        dists, counts = calc_all_distances_dict(sim_param, summary_stat, ssk; AD_mod=AD_mod)
        dists_used, dists_used_w = Dict{String,Float64}(), Dict{String,Float64}()
        for (i,key) in enumerate(dists_include)
            dists_used[key] = dists[key]
            dists_used_w[key] = dists[key]*weights["all"][key]
        end
        dists_used_keys = keys(dists_used_w)
        dists_used_vals, dists_used_vals_w = values(dists_used), values(dists_used_w)

        # Write the distances to file:
        println(f, "Counts: ", counts["Nmult1"], [counts["n_pl1"], counts["n_pairs1"]])
        println(f, "d_used_keys: ", dists_used_keys)
        println(f, "d_used_vals: ", dists_used_vals, [sum(dists_used_vals)])
        println(f, "d_used_vals_w: ", dists_used_vals_w, [sum(dists_used_vals_w)])

        d_tot_w = sum([sum(x) for x in dists_used_vals_w])

        # Write the total distances to file:
        println(f, "Total_dist_w: ", [d_tot_w])
        println(f, "#")



        # To save the catalog if the weighted distance passes the distance threshold:
        if d_tot_w <= d_threshold
            pass_count += 1
            if save_count < n_save
                save_count += 1
                println("$sim_count: d_tot_w = $d_tot_w; catalog saved (pass_count = $pass_count, save_count = $save_count)")

                save_physical_catalog_given_cat_phys(cat_phys_copy, sim_param; save_path=save_path, run_number=save_count)
                save_observed_catalog_given_cat_phys_obs(cat_phys, cat_obs, summary_stat, sim_param; save_path=save_path, run_number=save_count)
            else
                println("$sim_count: d_tot_w = $d_tot_w; (pass_count = $pass_count, max save_count reached)")
            end
        else
            println("$sim_count: d_tot_w = $d_tot_w")
        end

        summary_array = vcat(summary_array, reshape([[GP_points[sim_count,j] for j in 1:size(GP_points,2)]; d_tot_w - mean_f], (1,size(GP_points,2)+1)))
    end
end

println(f, "# elapsed time: ", t_elapsed, " seconds")
close(f)

summary_table = DataFrame(summary_array, [names(GP_points); "dist_tot_weighted"])
CSV.write(joinpath(save_path, "Simulate_GP_points_summary.txt"), summary_table)
