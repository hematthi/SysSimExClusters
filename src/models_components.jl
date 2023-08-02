if !@isdefined ExoplanetsSysSim
    using ExoplanetsSysSim
end

function get_fraction_of_stars_with_planets(star::StarAbstract, sim_param::SimParam)
    if haskey(sim_param, "f_stars_with_planets_attempted")
        f_swpa = get_real(sim_param, "f_stars_with_planets_attempted")
        @assert(0 <= f_swpa <= 1)
    else
        f_swpa = 1.
    end

    return f_swpa
end

function calc_fraction_of_stars_with_planets_linear_with_color(star::StarAbstract, sim_param::SimParam)
    # Calculate the fraction of stars with planets as a linear function of the stellar color:
    global stellar_catalog
    star_color = stellar_catalog[star.id, :bp_rp] - stellar_catalog[star.id, :e_bp_min_rp_interp]
    f_swpa_color_slope = get_real(sim_param, "f_stars_with_planets_attempted_color_slope")
    f_swpa_at_med_color = get_real(sim_param, "f_stars_with_planets_attempted_at_med_color")
    med_color = get_real(sim_param, "med_color")
    @assert(0 <= f_swpa_at_med_color <= 1)

    f_swpa = f_swpa_color_slope*(star_color - med_color) + f_swpa_at_med_color
    f_swpa = min(f_swpa, 1.)
    f_swpa = max(f_swpa, 0.)
    @assert(0 <= f_swpa <= 1)

    return f_swpa
end

function generate_num_clusters_poisson(s::Star, sim_param::SimParam)
    lambda::Float64 = exp(get_real(sim_param, "log_rate_clusters"))
    max_clusters_in_sys::Int64 = get_int(sim_param, "max_clusters_in_sys")
    return draw_truncated_poisson(lambda, min=0, max=max_clusters_in_sys, n=1)[1]
    #return ExoplanetsSysSim.generate_num_planets_poisson(lambda, max_clusters_in_sys) ##### Use this if setting max_clusters_in_sys > 20
end

function generate_num_clusters_ZTP(s::Star, sim_param::SimParam)
    lambda::Float64 = exp(get_real(sim_param, "log_rate_clusters"))
    max_clusters_in_sys::Int64 = get_int(sim_param, "max_clusters_in_sys")
    return draw_truncated_poisson(lambda, min=1, max=max_clusters_in_sys, n=1)[1]
end

function generate_num_planets_in_cluster_poisson(s::Star, sim_param::SimParam)
    lambda::Float64 = exp(get_real(sim_param, "log_rate_planets_per_cluster"))
    max_planets_in_cluster::Int64 = get_int(sim_param, "max_planets_in_cluster")
    return draw_truncated_poisson(lambda, min=0, max=max_planets_in_cluster, n=1)[1]
end

function generate_num_planets_in_cluster_ZTP(s::Star, sim_param::SimParam)
    lambda::Float64 = exp(get_real(sim_param, "log_rate_planets_per_cluster"))
    max_planets_in_cluster::Int64 = get_int(sim_param, "max_planets_in_cluster")
    return draw_truncated_poisson(lambda, min=1, max=max_planets_in_cluster, n=1)[1]
end

function generate_num_clusters_and_planets(s::Star, sim_param::SimParam)
    generate_num_clusters = get_function(sim_param, "generate_num_clusters")
    generate_num_planets_in_cluster = get_function(sim_param, "generate_num_planets_in_cluster")

    num_clusters = generate_num_clusters(s, sim_param)::Int64
    num_pl_in_cluster = map(x -> generate_num_planets_in_cluster(s, sim_param)::Int64, 1:num_clusters)
    return (num_clusters, num_pl_in_cluster)
end

function generate_stable_cluster(star::StarT, sim_param::SimParam; n::Int64=1) where {StarT<:StarAbstract}
    @assert(n >= 1)

    # Load functions and model parameters:
    generate_sizes = get_function(sim_param, "generate_sizes")
    min_radius::Float64 = get_real(sim_param, "min_radius")
    max_radius::Float64 = get_real(sim_param, "max_radius")
    sigma_log_radius_in_cluster = get_real(sim_param, "sigma_log_radius_in_cluster")
    generate_planet_mass_from_radius = get_function(sim_param, "generate_planet_mass_from_radius")

    min_period::Float64 = get_real(sim_param, "min_period")
    max_period::Float64 = get_real(sim_param, "max_period")
    max_period_ratio = max_period/min_period
    sigma_logperiod_per_pl_in_cluster = get_real(sim_param, "sigma_logperiod_per_pl_in_cluster")
    
    # Parameters for defining planets near MMRs:
    # TODO: replace with a definition using the zeta statistic
    w_mmrs = get_real(sim_param, "resonance_width")
    mmrs = get_any(sim_param, "period_ratios_mmr", Vector{Float64})

    # If single planet in cluster:
    if n==1
        R = generate_sizes(star, sim_param)
        mass = [generate_planet_mass_from_radius(R[1], sim_param)]
        P = [1.0] # generate_periods_power_law(star,sim_param)
        return (P, R, mass)
    end

    # If reach here, then at least 2 planets in cluster

    # Draw radii and masses:
    mean_R = generate_sizes(star,sim_param)[1]
    Rdist = Truncated(LogNormal(log(mean_R),sigma_log_radius_in_cluster), min_radius, max_radius) # clustered planet sizes
    R = rand(Rdist, n)
    mass = map(r -> generate_planet_mass_from_radius(r, sim_param), R)
    #println("# Rp = ", R)
    #println("# mass = ", mass)

    # Draw unscaled periods, checking for mutual-Hill stability (assuming circular orbits) of the entire cluster as a precondition:
    log_mean_P = 0.0
    Pdist = Truncated(LogNormal(log_mean_P,sigma_logperiod_per_pl_in_cluster*n), 1/sqrt(max_period_ratio), sqrt(max_period_ratio)) # truncated unscaled period distribution to ensure that the cluster can fit in the period range [min_period, max_period] after scaling by a period scale
    local P

    # Draw unscaled periods first, checking for mutual Hill separation stability assuming circular and coplanar orbits

    #= Allowed regions sampling:
    P = zeros(n)
    for i in 1:n # Draw periods one at a time
        if any(isnan.(P))
            P[i:end] .= NaN
            #println("Cannot fit any more planets in cluster.")
            break
        end
        P[i] = draw_period_lognormal_allowed_regions_mutualHill(P[1:i-1], mass[1:i-1], mass[i], star.mass, sim_param; μ=log_mean_P, σ=n*sigma_logperiod_per_pl_in_cluster, x_min=1/sqrt(max_period_ratio), x_max=sqrt(max_period_ratio))
    end
    found_good_periods = all(isnan.(P)) ? false : true
    @assert(test_stability(P, mass, star.mass, sim_param)) # should always be true if our unscaled period draws are correct
    =#

    # Rejection sampling:
    found_good_periods = false # will be true if entire cluster is likely to be stable assuming circular and coplanar orbits (given sizes/masses and periods)
    found_no_near_below_MMRs = true # will be true if there are no period ratios just within any MMRs in the cluster
    max_attempts = 1000
    attempts_periods = 0
    while (!found_good_periods || !found_no_near_below_MMRs) && attempts_periods < max_attempts
        attempts_periods += 1
        
        P = rand(Pdist, n)
        Pratios = P[2:end]./P[1:end-1]
        
        # Check the stability of the cluster:
        if test_stability(P, mass, star.mass, sim_param)
            # If pass mutual Hill criteria, also check circular MMR overlap criteria:
            a = semimajor_axis.(P, mass .+star.mass)
            μ = mass ./star.mass
            if test_stability_circular_MMR_overlap(μ, a)
                found_good_periods = true
            else
                @info("Found set of periods passing mutual Hill criteria but not circular MMR overlap criteria.")
            end
        end
        
        # Check if there are any planets just within MMRs:
        # TODO: replace with a definition using the zeta statistic
        #=
        has_near_below_MMRs = false
        for pr_mmr in mmrs
            has_near_below_MMRs = has_near_below_MMRs || any(pr_mmr*(1-w_mmrs) .< Pratios .< pr_mmr)
        end
        found_no_near_below_MMRs = !has_near_below_MMRs
        =#
    end # while trying to draw periods
    
    if attempts_periods == max_attempts
        P[1:end] .= NaN
        @info("Max attempts ($attempts_periods) reached. Returning an empty cluster.")
    end
    #

    return (P, R, mass) # NOTE: can also return earlier if only one planet in cluster; also, the planets are NOT sorted at this point
end

function move_period_ratios_from_narrow_to_wide_of_mmrs_inside_out!(Plist::Vector{Float64}, sim_param::SimParam; truncate_max_period::Bool=true)
    # NOTE: 'Plist' must be sorted!

    # Load functions and model parameters:
    # TODO: replace with a definition using the zeta statistic
    # TODO: add a parameter for the probability/fraction of such period ratios to move
    max_period::Float64 = get_real(sim_param, "max_period")
    w_mmrs = get_real(sim_param, "resonance_width")
    mmrs = get_any(sim_param, "period_ratios_mmr", Vector{Float64})
    @assert(all(mmrs .> 1))
    
    Pratios = Plist[2:end]./Plist[1:end-1]
    for (i,pr) in enumerate(Pratios)
        for pr_mmr in mmrs
            if pr_mmr*(1-w_mmrs) < pr < pr_mmr
                Plist[i+1:end] *= pr_mmr/pr # shift all the periods beyond the i^th planet to longer periods such that the period ratio of the i^th pair is set to 'pr_mmr'
            end
        end
    end
    
    if truncate_max_period
        Plist[Plist .>= max_period] .= NaN
    end
end

function generate_resonant_chain(star::StarT, sim_param::SimParam; n::Int64=2) where {StarT<:StarAbstract}
    @assert(n >= 2)

    # Load functions and model parameters:
    generate_sizes = get_function(sim_param, "generate_sizes")
    min_radius::Float64 = get_real(sim_param, "min_radius")
    max_radius::Float64 = get_real(sim_param, "max_radius")
    sigma_log_radius_in_cluster = get_real(sim_param, "sigma_log_radius_in_cluster")
    generate_planet_mass_from_radius = get_function(sim_param, "generate_planet_mass_from_radius")

    power_law_P = get_real(sim_param, "power_law_P")
    min_period::Float64 = get_real(sim_param, "min_period")
    max_period::Float64 = get_real(sim_param, "max_period")
    max_period_ratio = max_period/min_period
    mmrs = get_any(sim_param, "period_ratios_mmr", Vector{Float64})
    @assert(all(mmrs .> 1))

    # TODO: consider checking our stability criteria and putting this into a loop to resample until they are met
    # Draw radii and masses:
    mean_R = generate_sizes(star,sim_param)[1]
    Rdist = Truncated(LogNormal(log(mean_R),sigma_log_radius_in_cluster), min_radius, max_radius) # clustered planet sizes
    R = rand(Rdist, n)
    mass = map(r -> generate_planet_mass_from_radius(r, sim_param), R)

    # Draw the period ratios to set the (unscaled) periods:
    P_chain = [1.] # initialize the period list with the first planet, set to 1
    for i in 1:n-1
        period_ratio = rand(mmrs) # randomly select an MMR from the list of MMRs
        append!(P_chain, P_chain[i]*period_ratio)
    end

    # Discard any planets that would make the chain too long to fit in our period range:
    # Alternative: can also consider drawing the period of the first planet and then discarding any planets that don't fit within the maximum period
    keep::Array{Bool,1} = P_chain .< max_period_ratio
    R = R[keep]
    mass = mass[keep]
    P_chain = P_chain[keep]
    if length(P_chain) < n
        @info("Truncating resonant chain to fit in our period range ($n ==> $(length(P_chain)) planets).")
    end

    # Draw a period scale such that the chain fits in our period range:
    min_period_scale::Float64 = min_period/minimum(P_chain)
    max_period_scale::Float64 = max_period/maximum(P_chain)
    period_scale = draw_power_law(power_law_P, min_period_scale, max_period_scale, 1)
    P_chain .*= period_scale

    return (P_chain, R, mass) # the orbital periods are sorted, by design
end

function draw_planetary_system_orbits_and_sky_orientation_by_distributing_amd(Plist::Vector{Float64}, masslist::Vector{Float64}, star::Star, sim_param::SimParam; verbose::Bool=false)

    num_pl = length(Plist)
    @assert(num_pl == length(masslist))

    # Compute the critical AMD for the system and distribute it between the planets (to draw eccentricities and inclinations):
    f_amd_crit = get_real(sim_param, "f_amd_crit")
    f_amd_crit = f_amd_crit + (1 - f_amd_crit)*rand() # uniform between f_amd_crit and 1
    @assert(0 <= f_amd_crit <= 1)

    ecclist::Array{Float64,1} = Array{Float64}(undef, num_pl)
    ωlist::Array{Float64,1} = Array{Float64}(undef, num_pl)
    Ωlist::Array{Float64,1} = Array{Float64}(undef, num_pl)
    meananomlist::Array{Float64,1} = Array{Float64}(undef, num_pl)
    inclmutlist::Array{Float64,1} = Array{Float64}(undef, num_pl)
    AMDlist::Array{Float64,1} = Array{Float64}(undef, num_pl)

    μlist = masslist ./ star.mass # mass ratios
    alist = map(P -> semimajor_axis(P, star.mass), Plist)
    if num_pl == 1
        generate_e_omega = get_function(sim_param, "generate_e_omega")
        sigma_ecc::Float64 = get_real(sim_param, "sigma_hk")
        ecc::Float64, ω::Float64 = generate_e_omega(sigma_ecc)
        AMDlist, ecclist, ωlist, inclmutlist = μlist .*sqrt.(alist) .*(1 - sqrt(1 - ecc^2)), [ecc], [ω], [0.]
    else
        AMDlist, ecclist, ωlist, inclmutlist = draw_ecc_incl_system_critical_AMD(μlist, alist; f_amd_crit=f_amd_crit, check_stability=false)
    end

    Ωlist, meananomlist = 2π .* rand(num_pl), 2π .* rand(num_pl) # relative to the reference plane

    # Assign a reference (invariant) plane for the system:
    # NOTE: aligning sky plane to x-y plane (z-axis is unit normal)
    vec_z = [0.,0.,1.]

    vec_sys = draw_random_normal_vector() # unit normal of reference plane
    incl_sys = calc_angle_between_vectors(vec_z, vec_sys) # inclination of reference plane (rad) relative to sky plane
    Ω_sys = calc_Ω_in_sky_plane(vec_sys) # argument of ascending node of reference plane (rad) relative to sky plane

    # Calculate the sky orientations of each orbit:
    inclskylist::Array{Float64,1} = Array{Float64}(undef, num_pl)
    Ωskylist::Array{Float64,1} = Array{Float64}(undef, num_pl)

    inclskylist, Ωskylist = calc_sky_incl_Ω_orbits_given_system_vector(inclmutlist, Ωlist, vec_sys)

    # At this point, the "mutual inclinations" (inclmutlist) and ascending nodes (Ωlist) are relative to the reference plane, but that is NOT the same as the system invariable plane
    # To calculate true mutual inclinations relative to the system invariable plane:
    if num_pl > 1
        vec_orb_list = [calc_orbit_vector_given_system_vector(inclmutlist[i], Ωlist[i], vec_sys) for i in 1:num_pl] # unit normals of each planet's orbital plane

        blist = alist .* sqrt.((1 .- ecclist).*(1 .+ ecclist)) # semi-minor axis of each planet's orbit
        Llist = masslist .* blist .* sqrt.(ExoplanetsSysSim.G_mass_sun_in_mks * star.mass ./ alist) # angular momentum (magnitude) of each planet's orbit, as calculated from the Vis-viva equation
        Lvec_sys = sum(Llist .* vec_orb_list) # angular momentum vector of the system
        vec_invariable = Lvec_sys ./ norm(Lvec_sys) # unit normal to system invariable plane

        inclinvariablelist = [calc_angle_between_vectors(vec_invariable, vec_orb) for vec_orb in vec_orb_list] # mutual inclinations relative to system invariable plane
        Ωinvariablelist = [calc_Ω_in_plane(vec_orb, vec_invariable) for vec_orb in vec_orb_list] # ascending nodes relative to system invariable plane
    else
        Lvec_sys = vec_sys
        vec_invariable = vec_sys
        inclinvariablelist = inclmutlist
        Ωinvariablelist = Ωlist
    end

    incl_invariable = calc_angle_between_vectors(vec_z, vec_invariable) # inclination of system invariable plane (rad) relative to sky plane
    Ω_invariable = calc_Ω_in_sky_plane(vec_invariable) # argument of ascending node of system invariable plane (rad) relative to sky plane

    #AMD_crit = critical_AMD_system(μlist, alist)
    #AMD_tot_assigned = sum(AMDlist)
    #AMD_tot = total_AMD_system(μlist, alist, ecclist, inclinvariablelist)[1]

    #= For debugging:
    println("P (d): ", Plist)
    println("R (R⊕): ", Rlist ./ ExoplanetsSysSim.earth_radius)
    println("M (M⊕): ", masslist ./ ExoplanetsSysSim.earth_mass)
    println("e: ", ecclist)
    println("ω (deg): ", ωlist .* 180/π)
    println("i_m (deg): ", inclmutlist .* 180/π)
    println("i (deg): ", inclskylist .* 180/π)
    println("AMD: ", AMDlist)
    =#

    return (ecclist, inclskylist, ωlist, Ωskylist, meananomlist, incl_invariable, Ω_invariable)
end
