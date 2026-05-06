dir_path = dirname(@__FILE__)

include(joinpath(dir_path, "../src/clusters.jl"))
include(joinpath(dir_path, "../src/planetary_catalog.jl"))

sim_param = setup_sim_param_model()
add_param_fixed(sim_param,"num_targets_sim_pass_one", 86760)





##### To set the save path and the grid of values for the parameter to vary:

save_path = "/Users/hematthi/Documents/GradSchool/Research/SysSim/Simulated_catalogs/Hybrid_NR20_AMD_model1/clustered_initial_masses/Fit_some8p1_KS/Params10_fix_highM/vary_muM"

param_vary = "mean_ln_mass"
param_grid = range(-0.65, 1.35; length=9)





##### To simulate a catalog for each set of best active parameters:

for run_number in 1:length(param_grid)
    println("Generating simulated catalog ", run_number)

    add_param_active(sim_param, param_vary, param_grid[run_number])

    @time cat_phys, cat_phys_cut, cat_obs, summary_stat = generate_and_save_physical_and_observed_catalogs(sim_param; save_path=save_path, run_number=run_number)
end
