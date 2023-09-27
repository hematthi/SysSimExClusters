##### Temporary file for coding up some functions for computing models from Neil & Rogers (2020)
##### TODO: move to other files/organize these functions
##### TODO: write better docstrings
##### TODO: turn some parameters into sim params

using Distributions
using ExoplanetsSysSim





##### Functions for computing the components of the NR20 model:

function mean_radius_and_scatter_given_mass_neil_rogers2020(M::Real; C::Real, M_break1::Real, M_break2::Real, γ0::Real, γ1::Real, γ2::Real, σ0::Real, σ1::Real, σ2::Real)

    # M: planet mass (Earth masses)
    # C: normalization (Earth radii)
    # M_break1: break point 1 in planet mass (between low and intermediate masses)
    # M_break2: break point 2 in planet mass (between intermediate and high masses)
    # γ0: power-law slope for low mass range
    # γ1: power-law slope for intermediate mass range
    # γ2: power-law slope for high mass range
    # σ0: scatter in radius for low mass range
    # σ1: scatter in radius for intermediate mass range
    # σ2: scatter in radius for high mass range
    
    # Power-laws in each mass regime (low/intermediate/high):
    μ0 = C*M^γ0
    μ1 = C*M_break1^(γ0-γ1)*M^γ1
    μ2 = C*M_break1^(γ0-γ1)*M_break2^(γ1-γ2)*M^γ2
    
    # Smooth between the power-laws using a logistic function:
    S1 = 1/(1 + exp(-5*(log(M)-log(M_break1))))
    S2 = 1/(1 + exp(-5*(log(M)-log(M_break2))))
    
    μ = (1-S1)*μ0 + S1*(1-S2)*μ1 + S1*S2*μ2
    σ = (1-S1)*σ0 + S1*(1-S2)*σ1 + S1*S2*σ2
    
    return (μ, σ)
end

function mean_radius_and_scatter_given_mass_neil_rogers2020_one_break(M::Real; C::Real, M_break::Real, γ0::Real, γ1::Real, σ0::Real, σ1::Real)

    # M: planet mass (Earth masses)
    # C: normalization (Earth radii)
    # M_break: break point in planet mass
    # γ0: power-law slope for masses below break
    # γ1: power-law slope for masses above break
    # σ0: scatter in radius for masses below break
    # σ1: scatter in radius for masses above break
    
    # Power-laws in each mass regime:
    μ0 = C*M^γ0
    μ1 = C*M_break^(γ0-γ1)*M^γ1
    
    # Smooth between the power-laws using a logistic function:
    S = 1/(1 + exp(-5*(log(M)-log(M_break))))
    
    μ = (1-S)*μ0 + S*μ1
    σ = (1-S)*σ0 + S*σ1
    
    return (μ, σ)
end


function draw_radius_given_mass_neil_rogers2020(μ::Real, σ::Real)
    Rdist = Normal(μ, σ*μ) # normal distribution with a fractional scatter
    return rand(Rdist)
end

function draw_radius_given_mass_neil_rogers2020(M::Real; C::Real, M_break1::Real, M_break2::Real, γ0::Real, γ1::Real, γ2::Real, σ0::Real, σ1::Real, σ2::Real)
    # Wrapper function for the above
    
    μ, σ = mean_radius_and_scatter_given_mass_neil_rogers2020(M; C=C, M_break1=M_break1, M_break2=M_break2, γ0=γ0, γ1=γ1, γ2=γ2, σ0=σ0, σ1=σ1, σ2=σ2)
    return draw_radius_given_mass_neil_rogers2020(μ, σ)
end

function draw_radius_given_mass_neil_rogers2020_one_break(M::Real; C::Real, M_break::Real, γ0::Real, γ1::Real, σ0::Real, σ1::Real)
    
    μ, σ = mean_radius_and_scatter_given_mass_neil_rogers2020_one_break(M; C=C, M_break=M_break, γ0=γ0, γ1=γ1, σ0=σ0, σ1=σ1)
    return draw_radius_given_mass_neil_rogers2020(μ, σ)
end


function radius_given_mass_analytical_relation_seager2007(M::Real; R1::Real, M1::Real, k1::Real, k2::Real, k3::Real)
    R = R1*10^(k1 + (1/3)*log10(M/M1) - k2*(M/M1)^k3)
    return R
end

radius_given_mass_pure_iron_fit_seager2007(M::Real) = radius_given_mass_analytical_relation_seager2007(M; R1=2.52, M1=5.80, k1=-0.20949, k2=0.0804, k3=0.394)

radius_given_mass_pure_silicate_fit_seager2007(M::Real) = radius_given_mass_analytical_relation_seager2007(M; R1=3.90, M1=10.55, k1=-0.209594, k2=0.0799, k3=0.413)


σ_cgs = 5.670374419e−5 # Stefan-Boltzmann constant (erg cm^-2 s^-1 K^-4)

function luminosity_star(R_star::Real, T_eff::Real)

    # R_star: stellar radius (Solar radii)
    # T_eff: stellar effective temperature (K)
    
    R_star_in_cm = R_star*ExoplanetsSysSim.sun_radius_in_m_IAU2015 * 100.
    L = 4π*R_star_in_cm^2. * σ_cgs*T_eff^4. # luminosity (erg s^-1)
    return L
end

function bolometric_flux_at_planet_period(P::Real; Mstar::Real=1., L::Real=luminosity_star(1., 5780.))

    # P: planet orbital period (days)
    # Mstar=1.: stellar mass (Solar masses)
    # L=luminosity_star(1., 5780.): stellar luminosity (erg s^-1)
    
    sma = ExoplanetsSysSim.semimajor_axis(P, Mstar) * ExoplanetsSysSim.AU_in_m_IAU2012 * 100. # cm
    F = L/(4π*sma^2) # erg s^-1 cm^-2
    return F
end

F_oplus = bolometric_flux_at_planet_period(365.25, Mstar=1., L=luminosity_star(1., 5780.)) # bolometric flux on Earth at present (1.369 erg s^-1 cm^-2)

function mass_loss_timescale_lopez2012(M_env::Real, R_prim::Real, F_p::Real, F_oplus::Real=F_oplus, F_XUV_E100::Real=504., ϵ::Real=0.1)

    # M_env: planet envelope mass (Earth masses)
    # R_prim: planet primordial radius (before any envelope stripping; Earth radii)
    # F_p: incident bolometric flux on the planet (erg s^-1 cm^-2)
    # F_oplus=F_oplus: incident bolometric flux on the Earth at present (erg s^-1 cm^-2)
    # F_XUV_E100=504.: XUV flux at Earth at 100 Myr (erg s^-1 cm^-2)
    # ϵ=0.1: mass-loss efficiency
    
    G_in_cgs = ExoplanetsSysSim.G_in_mks_IAU2015 * (100^3)/1000 # cm^3 g^-1 s^-2
    M_env_in_g = M_env * ExoplanetsSysSim.earth_mass * ExoplanetsSysSim.sun_mass_in_kg_IAU2010 * 1000. # g
    R_prim_in_cm = R_prim * ExoplanetsSysSim.earth_radius_eq_in_m_IAU2015 * 100. # cm
    t_loss_in_s = ((G_in_cgs*M_env_in_g^2)/(π*ϵ*R_prim_in_cm^3*F_XUV_E100)) * (F_oplus/F_p) # s
    t_loss = t_loss_in_s / (ExoplanetsSysSim.sec_in_day * ExoplanetsSysSim.day_in_year * 1e9) # Gyr
    return t_loss
end

function prob_retain_envelope_neil_rogers2020(τ::Real, t_loss::Real, α::Real)

    # τ: age of the star (Gyr)
    # t_loss: mass-loss time scale (Gyr) given by `mass_loss_timescale_lopez2012`
    # α: prefactor/fudge factor to scale

    p_ret = min(α*t_loss/τ, 1)
    return p_ret
end

function α_prefactor_scaled_to_nominal_values(ϵ::Real, F_XUV_E100::Real, τ::Real)
    # NOTE: while Neil & Rogers provide an equation (what is computed by this function) showing how this "α" prefactor is scaled to nominal values of the parameters for which these values are assumed, they treat "α" as a single free parameter and find α ~ 8 for all of their models
    
    # ϵ: mass-loss efficiency
    # F_XUV_E100: XUV flux at Earth at 100 Myr (erg s^-1 cm^-2) # typo in Eq. 15 of Neil & Rogers (2020)? Equation does not include cm^-2 but the text does
    # τ: age of the star (Gyr)
    
    α = (0.1/ϵ)*(504. /F_XUV_E100)*(5. /τ)
end


function envelope_mass_smoothed_low_high_neil_rogers2020(M::Real, M_transition::Real=20.)
    
    # M: planet mass (Earth masses)
    # M_transition=20: transition planet mass based on where the scaling from Thorngren et al. (2016) is valid (set to 20 Earth masses)
    
    # The envelope mass scalings given the total mass
    M_env_low = 0.1*M
    M_env_high = M - sqrt(M)
    
    # Smooth between the envelope mass scalings using a logistic function:
    # NOTE: this is interpreted based on the text from Neil & Rogers (2020); they do not explicitly state the equations used for this part of the model
    S = 1/(1 + exp(-5*(log(M)-log(M_transition))))
    
    M_env = (1-S)*M_env_low + S*M_env_high # Earth masses
    return M_env
end





##### Functions for simulating a multiplanet model with the NR20 components:

##### NOTE: will need to generate stable clusters and draw the period scales (using initial masses) BEFORE implementing the photoevaporation component (since photoevaporation depends on the period/flux)

function generate_stable_cluster_radius_from_mass(star::StarT, sim_param::SimParam; n::Int64=1) where {StarT<:StarAbstract}
    @assert(n >= 1)

    # Load functions and model parameters:
    min_period = get_real(sim_param, "min_period")
    max_period = get_real(sim_param, "max_period")
    max_period_ratio = max_period/min_period
    sigma_logperiod_per_pl_in_cluster = get_real(sim_param, "sigma_logperiod_per_pl_in_cluster")
    
    # Note: these params are in Earth units
    min_mass = get_real(sim_param, "min_mass")
    max_mass = get_real(sim_param, "max_mass")
    min_radius = get_real(sim_param, "min_radius")
    max_radius = get_real(sim_param, "max_radius")
    μ_mass = get_real(sim_param, "mean_ln_mass")
    σ_mass = get_real(sim_param, "sigma_ln_mass")
    C = get_real(sim_param, "norm_radius")
    M_break1 = get_real(sim_param, "break1_mass")
    #M_break2 = get_real(sim_param, "break2_mass")
    γ0 = get_real(sim_param, "power_law_γ0")
    γ1 = get_real(sim_param, "power_law_γ1")
    #γ2 = get_real(sim_param, "power_law_γ2")
    σ0 = get_real(sim_param, "power_law_σ0")
    σ1 = get_real(sim_param, "power_law_σ1")
    #σ2 = get_real(sim_param, "power_law_σ2")

    # Draw initial masses and radii:
    M_init_dist = Truncated(LogNormal(μ_mass, σ_mass), min_mass, max_mass)
    
    #M_init = rand(M_init_dist, n)
    #R_init = map(M -> draw_radius_given_mass_neil_rogers2020(M; C=C, M_break1=M_break1, M_break2=M_break2, γ0=γ0, γ1=γ1, γ2=γ2, σ0=σ0, σ1=σ1, σ2=σ2), M_init)
    
    M_init::Array{Float64,1} = Array{Float64}(undef, n)
    R_init::Array{Float64,1} = Array{Float64}(undef, n)
    for i in 1:n # rejection sampling so initial radii are also between the bounds
        radius_in_bounds = false
        while !radius_in_bounds
            M_init[i] = rand(M_init_dist)
            R_init[i] = draw_radius_given_mass_neil_rogers2020_one_break(M_init[i]; C=C, M_break=M_break1, γ0=γ0, γ1=γ1, σ0=σ0, σ1=σ1)
            radius_in_bounds = true ? (min_radius < R_init[i] < max_radius) : false
        end
    end
    
    # Convert masses and radii to Solar units at this point:
    M_init *= ExoplanetsSysSim.earth_mass
    R_init *= ExoplanetsSysSim.earth_radius
    
    #println("# Rp_init = ", R_init)
    #println("# Mp_init = ", M_init)
    
    # If single planet in cluster:
    if n==1
        P = [1.0] # unscaled period
        return (P, R_init, M_init)
    end
    
    # If reach here, then at least 2 planets in cluster
    # Draw unscaled periods, checking for mutual-Hill stability (assuming circular orbits) of the entire cluster as a precondition:
    log_mean_P = 0.0
    Pdist = Truncated(LogNormal(log_mean_P,sigma_logperiod_per_pl_in_cluster*n), 1/sqrt(max_period_ratio), sqrt(max_period_ratio)) # truncated unscaled period distribution to ensure that the cluster can fit in the period range [min_period, max_period] after scaling by a period scale
    local P

    # Rejection sampling:
    found_good_periods = false # will be true if entire cluster is likely to be stable assuming circular and coplanar orbits (given sizes/masses and periods)
    max_attempts = 1000
    attempts_periods = 0
    while !found_good_periods && attempts_periods < max_attempts
        attempts_periods += 1
        
        P = rand(Pdist, n)
        Pratios = P[2:end]./P[1:end-1]
        
        # Check the stability of the cluster:
        if test_stability(P, M_init, star.mass, sim_param)
            # If pass mutual Hill criteria, also check circular MMR overlap criteria:
            a = semimajor_axis.(P, M_init .+ star.mass)
            μ = M_init ./star.mass
            if test_stability_circular_MMR_overlap(μ, a)
                found_good_periods = true
            else
                @info("Found set of periods passing mutual Hill criteria but not circular MMR overlap criteria.")
            end
        end
    end # while trying to draw periods
    
    if attempts_periods == max_attempts
        P[1:end] .= NaN
        #@info("Max attempts ($attempts_periods) reached. Returning an empty cluster.")
    end
    #

    return (P, R_init, M_init) # NOTE: can also return earlier if only one planet in cluster; also, the planets are NOT sorted at this point
end


function generate_planetary_system_clustered_periods_and_sizes_photoevap_distribute_amd(star::StarAbstract, sim_param::SimParam; verbose::Bool=false)

    # First, generate number of clusters (to attempt) and planets (to attempt) in each cluster:
    num_clusters::Int64, num_pl_in_cluster::Array{Int64,1} = generate_num_clusters_and_planets(star, sim_param)
    num_pl = sum(num_pl_in_cluster)
    num_pl_in_cluster_true = zeros(Int64, num_clusters) # true numbers of planets per cluster, after subtracting the number of NaNs
    #println("num_clusters: ", num_clusters, " ; num_pl_in_clusters", num_pl_in_cluster)

    # Return early if zero planets are drawn (NOTE: should never happen if draw clusters and planets per cluster using ZTPs):
    if num_pl==0
        return PlanetarySystem(star)
    end

    # Generate a set of periods, planet radii, and planet masses:
    power_law_P = get_real(sim_param, "power_law_P")
    min_period = get_real(sim_param, "min_period")
    max_period = get_real(sim_param, "max_period")

    clusteridlist::Array{Int64,1} = Array{Int64}(undef, num_pl)
    Plist::Array{Float64,1} = Array{Float64}(undef, num_pl)
    R_init_list::Array{Float64,1} = Array{Float64}(undef, num_pl)
    M_init_list::Array{Float64,1} = Array{Float64}(undef, num_pl)

    @assert(num_pl_in_cluster[1] >= 1)
    pl_start = 1
    pl_stop = 0
    for c in 1:num_clusters
        n = num_pl_in_cluster[c]
        pl_stop += n

        # Draw a stable cluster (with unscaled periods):
        Plist_tmp::Array{Float64,1}, R_init_list_tmp::Array{Float64,1}, M_init_list_tmp::Array{Float64,1} = generate_stable_cluster_radius_from_mass(star, sim_param, n=num_pl_in_cluster[c])

        clusteridlist[pl_start:pl_stop] = ones(Int64, num_pl_in_cluster[c])*c
        R_init_list[pl_start:pl_stop], M_init_list[pl_start:pl_stop] = R_init_list_tmp, M_init_list_tmp

        #= Allowed regions sampling:
        idx = .!isnan.(Plist[1:pl_stop-n])
        idy = .!isnan.(Plist_tmp)
        if any(idy)
            min_period_scale::Float64 = min_period/minimum(Plist_tmp[idy])
            max_period_scale::Float64 = max_period/maximum(Plist_tmp[idy])

            if min_period_scale < max_period_scale
                period_scale = draw_periodscale_power_law_allowed_regions_mutualHill(num_pl_in_cluster_true[1:c-1], Plist[1:pl_stop-n][idx], M_init_list[1:pl_stop-n][idx], Plist_tmp[idy], M_init_list_tmp[idy], star.mass, sim_param; x0=min_period_scale, x1=max_period_scale, α=power_law_P)
            else # cluster cannot fit at all
                period_scale = NaN
            end
        else # void cluster; all NaNs
            period_scale = NaN
        end
        Plist[pl_start:pl_stop] = Plist_tmp .* period_scale
        if isnan(period_scale)
            Plist[pl_stop:end] .= NaN
            #println("Cannot fit cluster into system; returning clusters that did fit.")
            break
        end
        @assert(test_stability(view(Plist,1:pl_stop), view(M_init_list,1:pl_stop), star.mass, sim_param)) # should always be true if our period scale draws are correct
        =#

        # Rejection sampling:
        valid_cluster = !any(isnan.(Plist_tmp)) # if the cluster has any nans, the whole cluster is discarded
        valid_period_scale = false
        max_attempts_period_scale = 100
        attempts_period_scale = 0
        min_period_scale::Float64 = min_period/minimum(Plist_tmp)
        max_period_scale::Float64 = max_period/maximum(Plist_tmp)
        while !valid_period_scale && attempts_period_scale<max_attempts_period_scale && valid_cluster
            attempts_period_scale += 1

            period_scale::Array{Float64,1} = draw_power_law(power_law_P, min_period_scale, max_period_scale, 1)

            Plist[pl_start:pl_stop] = Plist_tmp .* period_scale

            if test_stability(view(Plist,1:pl_stop), view(M_init_list,1:pl_stop), star.mass, sim_param)
                valid_period_scale = true
                #= If pass mutual Hill criteria, also check circular MMR overlap criteria: ##### WARNING: currently bugged (need to ignore NaNs)
                alist = semimajor_axis.(view(Plist,1:pl_stop), view(M_init_list,1:pl_stop) .+star.mass)
                μlist = view(M_init_list,1:pl_stop) ./star.mass
                if test_stability_circular_MMR_overlap(μlist, alist)
                    valid_period_scale = true
                else
                    @info("Found period scale passing mutual Hill criteria but not circular MMR overlap criteria.")
                end
                =#
            end
        end  # while !valid_period_scale...

        #if attempts_period_scale > 1
            #println("attempts_period_scale: ", attempts_period_scale)
        #end

        if !valid_period_scale
            Plist[pl_start:pl_stop] .= NaN
        end
        #

        num_pl_in_cluster_true[c] = sum(.!isnan.(Plist[pl_start:pl_stop]))
        pl_start += n
    end # for c in 1:num_clusters

    # Discard failed planets and sort the remaining planets by period:

    isnanPlist::Array{Bool,1} = isnan.(Plist::Array{Float64,1})
    if any(isnanPlist) # if any loop failed to generate valid planets, it should set a NaN in the period list
        keep::Array{Bool,1} = .!(isnanPlist) # currently, keeping clusters that could be fit, rather than throwing out entire systems and starting from scratch
        num_pl = sum(keep)

        clusteridlist = clusteridlist[keep]
        Plist = Plist[keep]
        R_init_list = R_init_list[keep]
        M_init_list = M_init_list[keep]
    end

    # Return the system (star) early if failed to draw any planets:
    if num_pl==0
        return PlanetarySystem(star)
    end

    idx = sortperm(Plist)
    clusteridlist = clusteridlist[idx]
    Plist = Plist[idx]
    R_init_list = R_init_list[idx]
    M_init_list = M_init_list[idx]
    
    # Implement envelope mass loss via photoevaporation:
    
    t_age = get_real(sim_param, "system_age") # age of system (Gyr)
    α = get_real(sim_param, "α_pret") # fudge factor for envelope retention probability/mass-loss timescale
    
    # NOTE/TODO: the functions for implementing photoevaporation use Earth units but the masses and radii are in Solar units at this point; figure out a way to do this more consistently and conveniently
    
    # Convert from Solar units to Earth units for this part:
    M_init_list /= ExoplanetsSysSim.earth_mass
    R_init_list /= ExoplanetsSysSim.earth_radius
    
    M_env_list::Array{Float64,1} = Array{Float64}(undef, num_pl)
    F_p_list::Array{Float64,1} = Array{Float64}(undef, num_pl)
    t_loss_list::Array{Float64,1} = Array{Float64}(undef, num_pl)
    p_ret_list::Array{Float64,1} = Array{Float64}(undef, num_pl)
    M_final_list::Array{Float64,1} = deepcopy(M_init_list)
    R_final_list::Array{Float64,1} = deepcopy(R_init_list)
    
    M_env_list = map(M -> envelope_mass_smoothed_low_high_neil_rogers2020(M), M_init_list) # initial envelope masses (Earth masses)
    Lstar = luminosity_star(star.radius, stellar_catalog[star.id,"teff"]) # WARNING: "star" struct does not have an effective temperature? Either need to get from stellar table or find another way to compute L (e.g. empirically with M^4)
    F_p_list = map(P -> bolometric_flux_at_planet_period(P, Mstar=star.mass, L=Lstar), Plist) # erg s^-1 cm^-2
    
    t_loss_list = map(i -> mass_loss_timescale_lopez2012(M_env_list[i], R_init_list[i], F_p_list[i]), 1:num_pl) # mass-loss timescales (Gyrs)
    p_ret_list = map(t -> prob_retain_envelope_neil_rogers2020(t_age, t, α), t_loss_list) # probabilities of retaining envelope
    bools_ret_list = rand(num_pl) .< p_ret_list # 1 = retain envelope, 0 = lose envelope
    
    M_final_list[.!bools_ret_list] = M_final_list[.!bools_ret_list] .- M_env_list[.!bools_ret_list] # final planet masses (Earth masses)
    R_final_list[.!bools_ret_list] = map(M -> draw_radius_given_mass_neil_rogers2020(radius_given_mass_pure_silicate_fit_seager2007(M), 0.05), M_final_list[.!bools_ret_list]) # final planet radii (Earth radii)
    
    # Convert back from Earth units to Solar units:
    M_init_list *= ExoplanetsSysSim.earth_mass
    R_init_list *= ExoplanetsSysSim.earth_radius
    M_env_list *= ExoplanetsSysSim.earth_mass
    M_final_list *= ExoplanetsSysSim.earth_mass
    R_final_list *= ExoplanetsSysSim.earth_radius
    
    
    # Draw eccentricities and inclinations and orient the system and the orbits of the planets:
    # NOTE: choosing to compute these using the initial planet masses
    ecclist, inclskylist, ωlist, Ωskylist, meananomlist, incl_invariable, Ω_invariable = draw_planetary_system_orbits_and_sky_orientation_by_distributing_amd(Plist, M_init_list, star, sim_param; verbose=verbose)

    # This final mutual-Hill stability test should be unnecessary if the period scales were drawn properly
    # NOTE: NOT including eccentricities in final mutual-Hill stability test since they were drawn by distributing AMD after the periods were set
    @assert(test_stability(Plist, M_init_list, star.mass, sim_param))

    # Packing and returning the system in structs:
    pl = Array{Planet}(undef, num_pl)
    orbit = Array{Orbit}(undef, num_pl)
    for i in 1:num_pl
        pl[i] = Planet(R_final_list[i], M_final_list[i], clusteridlist[i])
        orbit[i] = Orbit(Plist[i], ecclist[i], inclskylist[i], ωlist[i], Ωskylist[i], meananomlist[i])
    end
    sys_ref_plane = SystemPlane(incl_invariable, Ω_invariable)

    return PlanetarySystem(star, pl, orbit, sys_ref_plane)
end
