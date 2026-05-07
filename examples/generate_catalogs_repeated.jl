dir_path = dirname(@__FILE__)

include(joinpath(dir_path, "../src/clusters.jl"))
include(joinpath(dir_path, "../src/planetary_catalog.jl"))





##### To set the active parameters for repeated catalog simulation:

save_path = "/Users/hematthi/Documents/GradSchool/Research/SysSim/Simulated_catalogs/Hybrid_NR20_AMD_model1/Fit_some8_KS/Params9_fix_highM/Radius_valley_model89_repeated_100"

n_catalogs = 100 # number of catalogs to simulate

#####
sim_param = setup_sim_param_model()
add_param_fixed(sim_param,"num_targets_sim_pass_one", 86760)

add_param_active(sim_param, "log_rate_clusters", log(0.32))
add_param_active(sim_param, "log_rate_planets_per_cluster", log(3.11))
add_param_active(sim_param, "power_law_P", 0.94)
add_param_active(sim_param, "mean_ln_mass", 0.26) # ln(Earth masses)
add_param_active(sim_param, "sigma_ln_mass", 1.92) # ln(Earth masses)
add_param_active(sim_param, "norm_radius", 2.28) # Earth radii
add_param_active(sim_param, "power_law_γ0", 0.03)
add_param_active(sim_param, "power_law_σ0", 0.10)
add_param_active(sim_param, "log_α_pret", log(1.31))
#####





##### To simulate a catalog for each set of best active parameters:

for num_cat in 1:n_catalogs
    println("Generating simulated catalog ", num_cat)

    @time cat_phys, cat_phys_cut, cat_obs, summary_stat = generate_and_save_physical_and_observed_catalogs(sim_param; save_path=save_path, run_number=num_cat)
end
