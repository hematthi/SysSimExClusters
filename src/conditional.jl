function has_conditional_planet(ps::PlanetarySystem, sim_param::SimParam)

    # Check if system has a planet in the conditional period and radius range:
    cond_period_min = get_real(sim_param, "cond_period_min")
    cond_period_max = get_real(sim_param, "cond_period_max")
    cond_radius_min = get_real(sim_param, "cond_radius_min")
    cond_radius_max = get_real(sim_param, "cond_radius_max")
    cond_mass_min = get_real(sim_param, "cond_mass_min")
    cond_mass_max = get_real(sim_param, "cond_mass_max")
    cond_also_transits = get_bool(sim_param, "cond_also_transits")

    found_cond_planet = false
    for pl in 1:length(ps.planet)
        period = ps.orbit[pl].P
        radius = ps.planet[pl].radius
        mass = ps.planet[pl].mass
        in_period_range = cond_period_min <= period <= cond_period_max
        in_radius_range = cond_radius_min <= radius <= cond_radius_max
        in_mass_range = cond_mass_min <= mass <= cond_mass_max
        if in_period_range && in_radius_range && in_mass_range
            found_cond_planet = cond_also_transits ? ExoplanetsSysSim.does_planet_transit(ps, pl) : true
        end
    end
    return found_cond_planet
end

function draw_system_model_conditional(star::StarAbstract, sim_param::SimParam; verbose::Bool=false)
    local ps

    model_name = get(sim_param, "model_name", "")
    model_func = models[model_name]
    max_attempts_cond = get_int(sim_param, "max_attempts_cond")

    attempts_cond = 0
    found_cond_planet = false
    while !found_cond_planet && attempts_cond < max_attempts_cond
        attempts_cond += 1
        ps = model_func(star, sim_param; verbose=verbose)
        found_cond_planet = has_conditional_planet(ps, sim_param)
    end
    if found_cond_planet
        @info("Generated system with a conditional planet after $attempts_cond attempts.")
    else
        ps = PlanetarySystem(star)
        @info("Failed to generate a system with a conditional planet after $max_attempts_cond attempts; returning just the star.")
    end

    return ps
end
