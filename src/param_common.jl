if !@isdefined ExoplanetsSysSim
    using ExoplanetsSysSim
end

import Compat #: UTF8String, ASCIIString

function setup_sim_param_common()
    sim_param = SimParam()

    # For setting up how many targets to generate:
    add_param_fixed(sim_param, "num_targets_sim_pass_one", 86760) # the number of stars in the simulations, not necessarily related to number of Kepler targets
    add_param_fixed(sim_param, "num_kepler_targets", 86760) # the number of Kepler targets for the observational catalog to be compared

    # For generating target star properties:
    add_param_fixed(sim_param, "star_table_setup", ExoplanetsSysSim.StellarTable.setup_star_table)
    add_param_fixed(sim_param, "stellar_catalog", "q1q17_dr25_gaia_berger_fgk_HFR2020b.jld2") #"q1q17_dr25_gaia_fgk_interpolate_ebprp.jld2"
    add_param_fixed(sim_param, "generate_kepler_target", ExoplanetsSysSim.generate_kepler_target_from_table)
    add_param_fixed(sim_param, "window_function", "DR25topwinfuncs.jld2")
    #add_param_fixed(sim_param, "osd_file","dr25fgk_small_osds.jld2")
    add_param_fixed(sim_param, "osd_file","dr25fgk_relaxcut_osds.jld2") # WARNING: need 8gb of memory to read this file

    # For setting range of system sky-inclinations:
    # NOTE: make sure the difference between this and 90 (deg) is not too small compared to the mutual inclinations of the planets in each system
    add_param_fixed(sim_param,"max_incl_sys", 0.0) # degrees; gives system inclinations from "max_incl_sys" (deg) to 90 (deg), so set to 0 for isotropic distribution of system inclinations

    # For calculating observables from physical system properties:
    add_param_fixed(sim_param, "calc_target_obs_single_obs", ExoplanetsSysSim.calc_target_obs_single_obs)
    add_param_fixed(sim_param, "max_tranets_in_sys", 8) # SysSim ignores some planets in any systems with more than this many transiting planets to avoid wasting time on unphysical parameter values
    add_param_fixed(sim_param, "transit_noise_model", ExoplanetsSysSim.transit_noise_model_diagonal)
    add_param_fixed(sim_param, "read_target_obs", ExoplanetsSysSim.simulated_read_kepler_observations) # Read Kepler observations to compare to from disk
    #add_param_fixed(sim_param, "rng_seed",1234) # If you want to be able to reproduce simulations

    return sim_param
end

function add_sim_param_rates_of_planetary_systems_and_clusters_and_planets!(sim_param::SimParam)
    #add_param_active(sim_param, "f_stars_with_planets_attempted", 0.6)
    add_param_fixed(sim_param, "f_stars_with_planets_attempted_color_slope", 0.9)
    add_param_fixed(sim_param, "f_stars_with_planets_attempted_at_med_color", 0.88)
    add_param_fixed(sim_param, "med_color", 0.81)
    add_param_fixed(sim_param, "generate_num_clusters", generate_num_clusters_ZTP)
    add_param_fixed(sim_param, "generate_num_planets_in_cluster", generate_num_planets_in_cluster_ZTP)
    add_param_active(sim_param, "log_rate_clusters", log(0.9))
    add_param_fixed(sim_param, "max_clusters_in_sys", 20)
    add_param_active(sim_param, "log_rate_planets_per_cluster", log(1.65))
    add_param_fixed(sim_param, "max_planets_in_cluster", 20)
end

function add_sim_param_period_distribution!(sim_param::SimParam)
    add_param_fixed(sim_param, "generate_periods", ExoplanetsSysSim.generate_periods_power_law)
    add_param_fixed(sim_param, "min_period", 3.0)
    add_param_fixed(sim_param, "max_period", 300.0)
    #add_param_fixed(sim_param, "break_period", 10.0)
    add_param_active(sim_param, "power_law_P", 0.)
    #add_param_fixed(sim_param, "power_law_P1", 2.0)
    #add_param_fixed(sim_param, "power_law_P2", 0.5)
    add_param_fixed(sim_param, "sigma_logperiod_per_pl_in_cluster", 0.25)
end

function add_sim_param_radius_distribution!(sim_param::SimParam)
    add_param_fixed(sim_param, "generate_sizes", ExoplanetsSysSim.generate_sizes_broken_power_law) # if "generate_sizes_power_law", then takes "power_law_r"; if "generate_sizes_broken_power_law", then takes "power_law_r1", "power_law_r2", and "break_radius"
    add_param_fixed(sim_param, "min_radius", 0.5)
    add_param_fixed(sim_param, "max_radius", 10.0)
    add_param_fixed(sim_param, "break_radius", 3.0)
    #add_param_fixed(sim_param, "power_law_r", -2.5)
    add_param_active(sim_param, "power_law_r1", -1.4)
    add_param_active(sim_param, "power_law_r2", -5.2)
    add_param_active(sim_param, "sigma_log_radius_in_cluster", 0.3)
end

function add_sim_param_mass_radius_distribution!(sim_param::SimParam)
    #add_param_fixed(sim_param, "generate_planet_mass_from_radius", generate_planet_mass_from_radius_Ning2018_table)
    add_param_fixed(sim_param, "generate_planet_mass_from_radius", generate_planet_mass_from_radius_Ning2018_table_above_lognormal_mass_earthlike_rocky_below)
end

function add_sim_param_mass_and_radius_distribution_NR20!(sim_param::SimParam)
    add_param_fixed(sim_param, "mass_radius_model_name", "mass_radius_NR20")
    add_param_fixed(sim_param, "min_mass", 0.1) # Earth masses
    add_param_fixed(sim_param, "max_mass", 1e3) # Earth masses
    add_param_fixed(sim_param, "min_radius", 0.5) # Earth radii
    add_param_fixed(sim_param, "max_radius", 10.0) # Earth radii
    add_param_active(sim_param, "mean_ln_mass", 0.8) # ln(Earth masses)
    add_param_active(sim_param, "sigma_ln_mass", 1.3) # ln(Earth masses)
    add_param_active(sim_param, "sigma_ln_mass_in_cluster", 0.1) # ln(Earth masses)
    add_param_active(sim_param, "norm_radius", 1.96) # Earth radii
    add_param_fixed(sim_param, "break1_mass", 20.) # Earth masses
    #add_param_fixed(sim_param, "break2_mass", 175.7) # Earth masses
    add_param_active(sim_param, "power_law_γ0", 0.15)
    add_param_fixed(sim_param, "power_law_γ1", 0.5)
    #add_param_fixed(sim_param, "power_law_γ2", 0.04)
    add_param_active(sim_param, "power_law_σ0", 0.15)
    add_param_fixed(sim_param, "power_law_σ1", 0.3)
    #add_param_fixed(sim_param, "power_law_σ2", 0.10)
end

function add_sim_param_photoevaporation_NR20!(sim_param::SimParam)
    add_param_fixed(sim_param, "system_age", 5.) # Gyr
    add_param_active(sim_param, "log_α_pret", log(8.)) # fudge factor for the envelope retention probability
end

function add_sim_param_radius_valley_stats!(sim_param::SimParam)
    add_param_fixed(sim_param, "radius_valley_slope", -0.10) # slope in log(R[Earth radii])/log(P[days])
    add_param_fixed(sim_param, "radius_valley_offset", 2.40) # normalization (Earth radii) at P=1 day
    add_param_fixed(sim_param, "radius_valley_min_radius_include", 0.5) # Earth radii
    add_param_fixed(sim_param, "radius_valley_max_radius_include", 4.0) # Earth radii
    add_param_fixed(sim_param, "radius_valley_min_period_include", 3.0) # days
    add_param_fixed(sim_param, "radius_valley_max_period_include", 100.0) # days
end

function add_sim_param_eccentricity_distribution!(sim_param::SimParam)
    add_param_fixed(sim_param, "generate_e_omega", ExoplanetsSysSim.generate_e_omega_rayleigh)
    add_param_fixed(sim_param, "sigma_hk", 0.25)
end

function add_sim_param_inclination_distribution!(sim_param::SimParam)
    add_param_active(sim_param, "f_high_incl", 0.4) # fraction of systems with higher mutual inclinations
    add_param_active(sim_param, "sigma_incl", 50.) # degrees; 0 = coplanar w/ generate_kepler_target_simple; ignored by generate_planetary_system_uncorrelated_incl
    add_param_active(sim_param, "sigma_incl_near_mmr", 1.3) # degrees; for both the low mutual inclination population and planets near MMRs
end

function add_sim_param_stability_criteria_and_amd!(sim_param::SimParam)
    add_param_fixed(sim_param, "num_mutual_hill_radii", 10.0)
    add_param_fixed(sim_param, "f_amd_crit", 1.0) # fraction of critical AMD to distribute (can be greater than 1)
    #add_param_fixed(sim_param, "distribute_amd", distribute_amd_per_mass) # TODO: add param for defining how to distribute the AMD amongst the planets (per unit mass or randomly per planet)
end

function add_sim_param_resonant_chains!(sim_param::SimParam)
    add_param_fixed(sim_param, "resonance_width", 0.05)
    add_param_fixed(sim_param, "period_ratios_mmr", [2.0, 1.5, 4/3, 5/4])
end

function add_sim_param_conditionals!(sim_param::SimParam)
    add_param_fixed(sim_param, "max_attempts_cond", 10000)
    add_param_fixed(sim_param, "cond_period_min", 10.0)
    add_param_fixed(sim_param, "cond_period_max", 20.0)
    add_param_fixed(sim_param, "cond_radius_min", 1.5)
    add_param_fixed(sim_param, "cond_radius_max", 2.0)
    add_param_fixed(sim_param, "cond_mass_min", 0.0)
    add_param_fixed(sim_param, "cond_mass_max", 10.0)
    add_param_fixed(sim_param, "cond_also_transits", true)
end

function setup_sim_param_clustered_amd_model()
    sim_param = setup_sim_param_common()

    # Define the function name for the relevant model:
    add_param_fixed(sim_param, "model_name", "clustered_amd_model")

    # Add the necessary model parameters:
    add_sim_param_rates_of_planetary_systems_and_clusters_and_planets!(sim_param)
    add_sim_param_period_distribution!(sim_param)
    add_sim_param_radius_distribution!(sim_param)
    add_sim_param_mass_radius_distribution!(sim_param)
    add_sim_param_eccentricity_distribution!(sim_param) # for the single-planets only
    add_sim_param_stability_criteria_and_amd!(sim_param)
    add_sim_param_resonant_chains!(sim_param) # still needed for calculating which planets are near MMRs for now
    
    return sim_param
end

function setup_sim_param_clustered_photoevap_amd_model()
    sim_param = setup_sim_param_common()

    # Define the function name for the relevant model:
    add_param_fixed(sim_param, "model_name", "clustered_photoevap_amd_model")

    # Add the necessary model parameters:
    add_sim_param_rates_of_planetary_systems_and_clusters_and_planets!(sim_param)
    add_sim_param_period_distribution!(sim_param)
    add_sim_param_mass_and_radius_distribution_NR20!(sim_param)
    add_sim_param_photoevaporation_NR20!(sim_param)
    add_sim_param_eccentricity_distribution!(sim_param) # for the single-planets only
    add_sim_param_stability_criteria_and_amd!(sim_param)
    add_sim_param_resonant_chains!(sim_param) # still needed for calculating which planets are near MMRs for now
    
    # Add the parameters for computing the difference in planet radii from the radius valley location in period-radius:
    # NOTE: these parameters are only used for computing the summary statistic, NOT in the model/simulations!
    add_sim_param_radius_valley_stats!(sim_param)
    
    return sim_param
end

function setup_sim_param_resonant_chain_amd_model()
    sim_param = setup_sim_param_common()

    # Define the function name for the relevant model:
    add_param_fixed(sim_param, "model_name", "resonant_chain_amd_model")
    
    # Add the necessary model parameters:
    add_sim_param_rates_of_planetary_systems_and_clusters_and_planets!(sim_param)
    add_sim_param_period_distribution!(sim_param)
    add_sim_param_radius_distribution!(sim_param)
    add_sim_param_mass_radius_distribution!(sim_param)
    add_sim_param_eccentricity_distribution!(sim_param) # for the single-planets only
    add_sim_param_stability_criteria_and_amd!(sim_param)
    add_sim_param_resonant_chains!(sim_param)

    return sim_param
end

function setup_sim_param_clustered_and_resonant_chain_amd_mixture_model()
    sim_param = setup_sim_param_common()

    # Define the function name for the relevant model:
    add_param_fixed(sim_param, "model_name", "clustered_and_resonant_chain_amd_mixture_model")
    
    # Add the necessary model parameters:
    add_sim_param_rates_of_planetary_systems_and_clusters_and_planets!(sim_param)
    add_sim_param_period_distribution!(sim_param)
    add_sim_param_radius_distribution!(sim_param)
    add_sim_param_mass_radius_distribution!(sim_param)
    add_sim_param_eccentricity_distribution!(sim_param) # for the single-planets only
    add_sim_param_stability_criteria_and_amd!(sim_param)
    add_sim_param_resonant_chains!(sim_param)
    add_param_fixed(sim_param, "f_resonant_chains", 0.1)

    return sim_param
end

function setup_sim_param_model()
    # Setup the sim params for the desired model:
    #sim_param = setup_sim_param_clustered_amd_model()
    sim_param = setup_sim_param_clustered_photoevap_amd_model()
    #sim_param = setup_sim_param_resonant_chain_amd_model()
    #sim_param = setup_sim_param_clustered_and_resonant_chain_amd_mixture_model()

    add_param_fixed(sim_param, "generate_planetary_system", draw_system_model)

    # To condition the model on a given planet, uncomment the following:
    #add_sim_param_conditionals!(sim_param::SimParam)

    return sim_param
end

function test_setup_sim_param()
    setup_sim_param_model()
end



function write_model_params(f, sim_param::SimParam)
    # Write the model parameters to a file "f" as a header

    # Common simulation parameters:
    println(f, "# num_targets_sim_pass_one: ", get_int(sim_param, "num_targets_sim_pass_one"))
    println(f, "# stellar_catalog: ", get(sim_param, "stellar_catalog", ""))
    println(f, "# osd_file: ", get(sim_param, "osd_file", ""))
    println(f, "# max_incl_sys: ", get_real(sim_param, "max_incl_sys"))
    if haskey(sim_param, "rng_seed")
        println(f, "# rng_seed: ", get_int(sim_param, "rng_seed"))
    end

    # Model specific parameters:
    model_name = get(sim_param, "model_name", "")
    println(f, "# model_name: ", model_name)
    if model_name == "clustered_and_resonant_chain_amd_mixture_model"
        println(f, "# f_resonant_chains: ", get_real(sim_param, "f_resonant_chains"))
    end

    #println(f, "# f_stars_with_planets_attempted: ", get_real(sim_param, "f_stars_with_planets_attempted"))
    println(f, "# f_stars_with_planets_attempted_at_med_color: ", get_real(sim_param, "f_stars_with_planets_attempted_at_med_color"))
    println(f, "# f_stars_with_planets_attempted_color_slope: ", get_real(sim_param, "f_stars_with_planets_attempted_color_slope"))
    println(f, "# med_color: ", get_real(sim_param, "med_color"))
    println(f, "# generate_num_clusters: ", string(get_function(sim_param, "generate_num_clusters")))
    println(f, "# log_rate_clusters: ", get_real(sim_param, "log_rate_clusters"))
    println(f, "# max_clusters_in_sys: ", get_int(sim_param, "max_clusters_in_sys"))
    println(f, "# generate_num_planets_in_cluster: ", string(get_function(sim_param, "generate_num_planets_in_cluster")))
    println(f, "# log_rate_planets_per_cluster: ", get_real(sim_param, "log_rate_planets_per_cluster"))
    println(f, "# max_planets_in_clusters: ", get_int(sim_param, "max_planets_in_cluster"))

    generate_periods_func = string(get_function(sim_param, "generate_periods"))
    println(f, "# generate_periods: ", generate_periods_func)
    println(f, "# min_period: ", get_real(sim_param, "min_period"))
    println(f, "# max_period: ", get_real(sim_param, "max_period"))
    if generate_periods_func == "generate_periods_power_law"
        println(f, "# power_law_P: ", get_real(sim_param, "power_law_P"))
    elseif generate_periods_func == "generate_periods_broken_power_law"
        println(f, "# break_period: ", get_real(sim_param, "break_period"))
        println(f, "# power_law_P1: ", get_real(sim_param, "power_law_P1"))
        println(f, "# power_law_P2: ", get_real(sim_param, "power_law_P2"))
    end
    println(f, "# sigma_logperiod_per_pl_in_cluster: ", get_real(sim_param, "sigma_logperiod_per_pl_in_cluster"))

    ##### TODO: clean this up:
    mass_radius_model_name = get(sim_param, "mass_radius_model_name", "")
    if mass_radius_model_name == "mass_radius_NR20"
        println(f, "# min_mass (M_earth): ", get_real(sim_param, "min_mass"))
        println(f, "# max_mass (M_earth): ", get_real(sim_param, "max_mass"))
        println(f, "# min_radius (R_earth): ", get_real(sim_param, "min_radius"))
        println(f, "# max_radius (R_earth): ", get_real(sim_param, "max_radius"))
        println(f, "# mean_ln_mass (ln M_earth): ", get_real(sim_param, "mean_ln_mass"))
        println(f, "# sigma_ln_mass (ln M_earth): ", get_real(sim_param, "sigma_ln_mass"))
        println(f, "# sigma_ln_mass_in_cluster (ln M_earth): ", get_real(sim_param, "sigma_ln_mass_in_cluster"))
        println(f, "# norm_radius (R_earth): ", get_real(sim_param, "norm_radius"))
        println(f, "# break_mass (M_earth): ", get_real(sim_param, "break1_mass"))
        println(f, "# power_law_γ0: ", get_real(sim_param, "power_law_γ0"))
        println(f, "# power_law_γ1: ", get_real(sim_param, "power_law_γ1"))
        println(f, "# power_law_σ0: ", get_real(sim_param, "power_law_σ0"))
        println(f, "# power_law_σ1: ", get_real(sim_param, "power_law_σ1"))
        
        println(f, "# system_age (Gyr): ", get_real(sim_param, "system_age"))
        println(f, "# log_α_pret: ", get_real(sim_param, "log_α_pret"))
    else
        #generate_sizes_func = string(get_function(sim_param, "generate_sizes"))
        generate_sizes_func = "" # avoid writing M-R params for now
        println(f, "# generate_sizes: ", generate_sizes_func)
        println(f, "# min_radius (R_earth): ", get_real(sim_param, "min_radius"))
        println(f, "# max_radius (R_earth): ", get_real(sim_param, "max_radius"))
        if generate_sizes_func == "generate_sizes_power_law"
            println(f, "# power_law_r: ", get_real(sim_param, "power_law_r"))
        elseif generate_sizes_func == "generate_sizes_broken_power_law"
            println(f, "# break_radius (R_earth): ", get_real(sim_param, "break_radius"))
            println(f, "# power_law_r1: ", get_real(sim_param, "power_law_r1"))
            println(f, "# power_law_r2: ", get_real(sim_param, "power_law_r2"))
        end
        #println(f, "# sigma_log_radius_in_cluster: ", get_real(sim_param, "sigma_log_radius_in_cluster"))

        #generate_masses_func = string(get_function(sim_param, "generate_planet_mass_from_radius"))
        generate_masses_func = "" # avoid writing M-R params for now
        if generate_masses_func == "generate_planet_mass_from_radius_powerlaw"
            println(f, "# mr_model: ", generate_masses_func)
            println(f, "# mr_power_index: ", get_real(sim_param, "mr_power_index"))
            println(f, "# mr_max_mass (M_earth): ", get_real(sim_param, "mr_max_mass")/ExoplanetsSysSim.earth_mass)
        elseif generate_masses_func == "generate_planet_mass_from_radius_Ning2018"
            println(f, "# mr_model: Ning2018")
        elseif generate_masses_func == "generate_planet_mass_from_radius_Ning2018_table"
            println(f, "# mr_model: Ning2018_table")
        elseif generate_masses_func == "generate_planet_mass_from_radius_Ning2018_table_above_normal_density_earthlike_rocky_below"
            println(f, "# radius_switch (R_earth): ", radius_switch)
            println(f, "# mr_model: Ning2018_table (above radius_switch)")
            println(f, "# mr_model: Normal density around Earthlike rocky (below radius_switch)")
        elseif generate_masses_func == "generate_planet_mass_from_radius_Ning2018_table_above_lognormal_mass_earthlike_rocky_below"
            println(f, "# radius_switch (R_earth): ", radius_switch)
            println(f, "# mr_model: Ning2018_table (above radius_switch)")
            println(f, "# mr_model: Lognormal mass around Earthlike rocky (below radius_switch)")
        end
    end
    ##### end TODO

    println(f, "# generate_e_omega: ", string(get_function(sim_param, "generate_e_omega")))
    println(f, "# sigma_hk: ", get_real(sim_param, "sigma_hk"))
    #println(f, "# sigma_hk_at_med_color: ", get_real(sim_param, "sigma_hk_at_med_color"))
    #println(f, "# sigma_hk_color_slope: ", get_real(sim_param, "sigma_hk_color_slope"))
    #println(f, "# f_high_incl: ", get_real(sim_param, "f_high_incl"))
    #println(f, "# sigma_incl (deg): ", get_real(sim_param, "sigma_incl"))
    #println(f, "# sigma_incl_near_mmr (deg): ", get_real(sim_param, "sigma_incl_near_mmr"))
    println(f, "# num_mutual_hill_radii: ", get_real(sim_param, "num_mutual_hill_radii"))
    println(f, "# f_amd_crit: ", get_real(sim_param, "f_amd_crit"))
    #println(f, "# distribute_amd: ", string(get_function(sim_param, "distribute_amd"))) # TODO: add param

    println(f, "# resonance_width: ", get_real(sim_param, "resonance_width"))
    println(f, "# period_ratios_mmr: ", get(sim_param, "period_ratios_mmr", Float64[]))

    # Conditionals, if conditioning on a given planet:
    if haskey(sim_param, "max_attempts_cond")
        println(f, "# max_attempts_cond: ", get_int(sim_param, "max_attempts_cond"))
        println(f, "# cond_period_min: ", get_real(sim_param, "cond_period_min"))
        println(f, "# cond_period_max: ", get_real(sim_param, "cond_period_max"))
        println(f, "# cond_radius_min (R_earth): ", get_real(sim_param, "cond_radius_min"))
        println(f, "# cond_radius_max (R_earth): ", get_real(sim_param, "cond_radius_max"))
        println(f, "# cond_mass_min (M_earth): ", get_real(sim_param, "cond_mass_min"))
        println(f, "# cond_mass_max (M_earth): ", get_real(sim_param, "cond_mass_max"))
        println(f, "# cond_also_transits: ", get_bool(sim_param, "cond_also_transits"))
    end

    println(f, "#")
end
