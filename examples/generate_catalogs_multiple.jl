dir_path = dirname(@__FILE__)

include(joinpath(dir_path, "../src/clusters.jl"))
include(joinpath(dir_path, "../src/planetary_catalog.jl"))

sim_param = setup_sim_param_model()
add_param_fixed(sim_param,"num_targets_sim_pass_one", 86760)





##### To load a file with all the best active parameters from a set of optimization runs:

save_path = "/Users/hematthi/Documents/GradSchool/Research/SysSim/Simulated_catalogs/Hybrid_NR20_AMD_model1/clustered_initial_masses/Fit_some8p1_KS/Params10_fix_highM/optim_best_100"
file_path = "/Users/hematthi/Documents/NPP_ARC_Modernize_Kepler/Personal_research/SysSim/Model_Optimization/Hybrid_NR20_AMD_model1/clustered_initial_masses/Fit_some8p1_KS/Params10_fix_highM/best_N"

file_name = "Active_params_table_best100.txt"

optim_points = CSV.read(joinpath(file_path, file_name), DataFrame, delim=" ")
active_params_names = names(optim_points)[2:end]
active_params_best_all = optim_points[!,active_params_names]





##### To simulate a catalog for each set of best active parameters:

for run_number in 1:size(active_params_best_all)[1]
    println("Generating simulated catalog ", run_number)

    # To set up the model parameters:
    for (i,param_name) in enumerate(active_params_names)
        add_param_active(sim_param, string(param_name), active_params_best_all[run_number, param_name])
    end

    @time cat_phys, cat_phys_cut, cat_obs, summary_stat = generate_and_save_physical_and_observed_catalogs(sim_param; save_path=save_path, run_number=run_number)
end
