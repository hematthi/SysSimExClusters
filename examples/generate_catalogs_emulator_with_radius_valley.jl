dir_path = dirname(@__FILE__)

include(joinpath(dir_path, "../src/clusters.jl"))
include(joinpath(dir_path, "../src/planetary_catalog.jl"))
include(joinpath(dir_path, "../src/optimization.jl"))



##### To define a function for measuring the depth of the radius valley:

using KernelDensity, PyPlot, Statistics

function measure_and_plot_radius_valley_depth_using_kde(
    radii_sample::Vector{Float64};
    radius_valley_bounds::Tuple{Float64,Float64} = (1.8, 2.2),
    x_min::Float64 = 0.5,
    x_max::Float64 = 5.5,
    n_bins::Int = 100,
    bw::Union{String,Float64} = "Scotts",
    bw_scotts_factor::Float64 = 1.0,
    fractional_depth::Bool = true,
    xlabel_text::String = "Planet radius, Rₚ [R_⊕]",
    verbose::Bool = false,
    plot_fig::Bool = false,
    save_name::String = "no_name_fig.pdf",
    save_fig::Bool = false)
    
    # Bandwidth selection:
    if bw == "Scotts"
        bw = length(radii_sample)^(-1/5) * bw_scotts_factor
    end

    # Fit KDE:
    kd = kde(radii_sample, bandwidth=bw)
    x_evals = range(x_min, x_max, length=1001)
    kde_evals = pdf(kd, collect(x_evals))

    # Find indices in the valley region:
    i_valley = findall(x -> (x >= radius_valley_bounds[1] && x <= radius_valley_bounds[2]), x_evals)

    # Find minimum in the valley, and maxima on either side:
    i_min_valley = i_valley[argmin(kde_evals[i_valley])]
    i_max_before_valley = argmax(kde_evals[1:i_valley[1]])
    i_max_after_valley = argmax(kde_evals[i_valley[end]:end]) + i_valley[end] - 1

    min_valley = kde_evals[i_min_valley]
    max_before_valley = kde_evals[i_max_before_valley]
    max_after_valley = kde_evals[i_max_after_valley]

    # Compute the "depth" of the radius valley:
    height = min(max_before_valley, max_after_valley)
    depth = fractional_depth ? (height - min_valley)/height : height - min_valley
    if depth < 0
        println("WARNING: the depth ($depth) is negative! There is no valley in the given radius valley bounds.")
    end

    if verbose
        println("min_valley = $min_valley at radius = $(x_evals[i_min_valley])")
        println("max_before_valley = $max_before_valley at radius = $(x_evals[i_max_before_valley])")
        println("max_after_valley = $max_after_valley at radius = $(x_evals[i_max_after_valley])")
        println("depth = $depth")
    end

    if plot_fig
        fig, ax = subplots()
        counts, bins, patches = ax.hist(radii_sample; bins=n_bins, density=true, range=(x_min, x_max), alpha=0.5, label="Data")
        ax.plot(x_evals, kde_evals, label="KDE fit (bw = $(round(bw, sigdigits=2)))", color="blue")
        ax.axvline(radius_valley_bounds[1], linestyle=":", color="black", lw=1)
        ax.axvline(radius_valley_bounds[2], linestyle=":", color="black", lw=1)
        ax.hlines(height, radius_valley_bounds[1], radius_valley_bounds[2], linestyles="dashed", color="gray", lw=1, label="Height")
        ax.scatter([x_evals[i_min_valley]], [min_valley], color="red", label="Valley min")
        ax.set_xlabel(xlabel_text)
        ax.set_ylabel("Density")
        ax.set_title("Radius Valley Depth using KDE")
        ax.legend(loc="upper right")
        if save_fig
            fig.savefig(save_name)
            close(fig)
        end
    end

    return depth
end





##### To load model parameters found using the GP emulator and simulate catalogs if they pass a distance threshold:

# Clustered initial masses, fit some+1 KS, 10 params:
run_path = "Hybrid_NR20_AMD_model1/clustered_initial_masses/Fit_some8p1_KS/Params10_fix_highM/"

GP_data_path = joinpath("/Users/hematthi/Documents/NPP_ARC_Modernize_Kepler/Personal_research/SysSim/Model_Optimization", run_path, "GP_files/")
GP_file_name = "GP_train2000_meanf30.0_sigmaf2.7_lscales3.85_vol9.72_points100000_meanInf_stdInf_post-15.0.csv"



GP_points = CSV.read(joinpath(GP_data_path, GP_file_name), DataFrame, comment="#")
active_params_names = names(GP_points)[1:end-3]
active_params_best_all = GP_points[!,active_params_names]

model_name = "Hybrid1"
AD_mod = true
num_targs = 86760
dists_include = ["delta_f", "mult_CRPD_r", "periods_KS", "depths_KS", "radii_KS", "radius_ratios_KS", "radii_partitioning_KS", "radii_monotonicity_KS"]
#dists_include = ["delta_f", "mult_CRPD_r", "periods_KS", "period_ratios_KS", "durations_KS", "duration_ratios_KS", "depths_KS", "radii_KS", "radius_ratios_KS", "radii_partitioning_KS", "radii_monotonicity_KS", "gap_complexity_KS"]

d_threshold, mean_f = 15., 30.
radii_delta_depth_threshold = 0.29
bw_factor = 0.25

n_pass = 1000 # number of simulations we want to pass the distance threshold
n_save = n_pass # number of simulations we want to pass the distance threshold and also save (choose a small number or else requires a lot of storage space); must not be greater than n_pass!

save_path = joinpath("/Users/hematthi/Documents/GradSchool/Research/SysSim/Simulated_catalogs", run_path, "GP_dtotmax$(d_threshold)_depthmin$(radii_delta_depth_threshold)_models/")

file_name = model_name*"_pass_GP_meanf$(mean_f)_thres$(d_threshold)_mindepth$(radii_delta_depth_threshold)_pass$(n_pass)_targs$(num_targs).txt"
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
summary_array = Array{Float64,2}(undef, 0, size(GP_points,2)+2)

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
        
        # Compute the measured depth of the radius valley:
        # NOTE: this is of the gap-substracted radii, for the restricted sample
        radii_delta_depth = measure_and_plot_radius_valley_depth_using_kde(summary_stat.stat["radii_delta_valley"], radius_valley_bounds=(-0.2,0.2), x_min=-1.5, x_max=2.5, bw_scotts_factor=bw_factor)

        # Write the distances to file:
        println(f, "Counts: ", counts["Nmult1"], [counts["n_pl1"], counts["n_pairs1"]])
        println(f, "d_used_keys: ", dists_used_keys)
        println(f, "d_used_vals: ", dists_used_vals, [sum(dists_used_vals)])
        println(f, "d_used_vals_w: ", dists_used_vals_w, [sum(dists_used_vals_w)])

        d_tot_w = sum([sum(x) for x in dists_used_vals_w])

        # Write the total distance and measured depth to file:
        println(f, "Total_dist_w: ", [d_tot_w])
        println(f, "radii_delta_depth: ", [radii_delta_depth])
        println(f, "#")



        # To save the catalog if the weighted distance passes the distance threshold:
        if (d_tot_w <= d_threshold) && (radii_delta_depth >= radii_delta_depth_threshold)
            pass_count += 1
            if save_count < n_save
                save_count += 1
                println("$sim_count: d_tot_w = $d_tot_w, depth = $radii_delta_depth; catalog saved (pass_count = $pass_count, save_count = $save_count)")

                save_physical_catalog_given_cat_phys(cat_phys_copy, sim_param; save_path=save_path, run_number=save_count)
                save_observed_catalog_given_cat_phys_obs(cat_phys, cat_obs, summary_stat, sim_param; save_path=save_path, run_number=save_count)
            else
                println("$sim_count: d_tot_w = $d_tot_w; (pass_count = $pass_count, max save_count reached)")
            end
        else
            println("$sim_count: d_tot_w = $d_tot_w, depth = $radii_delta_depth")
        end

        summary_array = vcat(summary_array, reshape([[GP_points[sim_count,j] for j in 1:size(GP_points,2)]; d_tot_w - mean_f; radii_delta_depth], (1,size(GP_points,2)+2)))
    end
end

println(f, "# elapsed time: ", t_elapsed, " seconds")
close(f)

summary_table = DataFrame(summary_array, [names(GP_points); "dist_tot_weighted"; "radii_delta_depth"])
CSV.write(joinpath(save_path, "Simulate_GP_points_summary.txt"), summary_table)
