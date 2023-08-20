include("models_components.jl")
include("models_test.jl")
include("conditional.jl")

function generate_planetary_system_clustered_periods_and_sizes_distribute_amd(star::StarAbstract, sim_param::SimParam; verbose::Bool=false)

    #local num_pl, clusteridlist, Plist, Rlist, masslist, ecclist, ωlist, Ωlist, meananomlist, inclmutlist, inclskylist, Ωskylist

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
    Rlist::Array{Float64,1} = Array{Float64}(undef, num_pl)
    masslist::Array{Float64,1} = Array{Float64}(undef, num_pl)

    @assert(num_pl_in_cluster[1] >= 1)
    pl_start = 1
    pl_stop = 0
    for c in 1:num_clusters
        n = num_pl_in_cluster[c]
        pl_stop += n

        # Draw a stable cluster (with unscaled periods):
        Plist_tmp::Array{Float64,1}, Rlist_tmp::Array{Float64,1}, masslist_tmp::Array{Float64,1} = generate_stable_cluster(star, sim_param, n=num_pl_in_cluster[c])

        clusteridlist[pl_start:pl_stop] = ones(Int64, num_pl_in_cluster[c])*c
        Rlist[pl_start:pl_stop], masslist[pl_start:pl_stop] = Rlist_tmp, masslist_tmp

        #= Allowed regions sampling:
        idx = .!isnan.(Plist[1:pl_stop-n])
        idy = .!isnan.(Plist_tmp)
        if any(idy)
            min_period_scale::Float64 = min_period/minimum(Plist_tmp[idy])
            max_period_scale::Float64 = max_period/maximum(Plist_tmp[idy])

            if min_period_scale < max_period_scale
                period_scale = draw_periodscale_power_law_allowed_regions_mutualHill(num_pl_in_cluster_true[1:c-1], Plist[1:pl_stop-n][idx], masslist[1:pl_stop-n][idx], Plist_tmp[idy], masslist_tmp[idy], star.mass, sim_param; x0=min_period_scale, x1=max_period_scale, α=power_law_P)
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
        @assert(test_stability(view(Plist,1:pl_stop), view(masslist,1:pl_stop), star.mass, sim_param)) # should always be true if our period scale draws are correct
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

            if test_stability(view(Plist,1:pl_stop), view(masslist,1:pl_stop), star.mass, sim_param)
                valid_period_scale = true
                #= If pass mutual Hill criteria, also check circular MMR overlap criteria: ##### WARNING: currently bugged (need to ignore NaNs)
                alist = semimajor_axis.(view(Plist,1:pl_stop), view(masslist,1:pl_stop) .+star.mass)
                μlist = view(masslist,1:pl_stop) ./star.mass
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
        Rlist = Rlist[keep]
        masslist = masslist[keep]
    end

    # Return the system (star) early if failed to draw any planets:
    if num_pl==0
        return PlanetarySystem(star)
    end

    idx = sortperm(Plist)
    clusteridlist = clusteridlist[idx]
    Plist = Plist[idx]
    Rlist = Rlist[idx]
    masslist = masslist[idx]
    
    ##### TESTING NEW MODEL:
    #= Move period ratios just narrow of MMRs to just wide (or at) MMRs:
    
    move_period_ratios_from_narrow_to_wide_of_mmrs_inside_out!(Plist, sim_param)
    
    keep = .!isnan.(Plist)
    num_pl = sum(keep)
    clusteridlist = clusteridlist[keep]
    Plist = Plist[keep]
    Rlist = Rlist[keep]
    masslist = masslist[keep]
    =#
    #####

    # Draw eccentricities and inclinations and orient the system and the orbits of the planets:
    ecclist, inclskylist, ωlist, Ωskylist, meananomlist, incl_invariable, Ω_invariable = draw_planetary_system_orbits_and_sky_orientation_by_distributing_amd(Plist, masslist, star, sim_param; verbose=verbose)

    # This final mutual-Hill stability test should be unnecessary if the period scales were drawn properly
    # NOTE: NOT including eccentricities in final mutual-Hill stability test since they were drawn by distributing AMD after the periods were set
    @assert(test_stability(Plist, masslist, star.mass, sim_param))

    # Packing and returning the system in structs:
    pl = Array{Planet}(undef, num_pl)
    orbit = Array{Orbit}(undef, num_pl)
    for i in 1:num_pl
        pl[i] = Planet(Rlist[i], masslist[i], clusteridlist[i])
        orbit[i] = Orbit(Plist[i], ecclist[i], inclskylist[i], ωlist[i], Ωskylist[i], meananomlist[i])
    end
    sys_ref_plane = SystemPlane(incl_invariable, Ω_invariable)

    return PlanetarySystem(star, pl, orbit, sys_ref_plane)
end

function generate_planetary_system_resonant_chain_clustered_sizes_distribute_amd(star::StarAbstract, sim_param::SimParam; verbose::Bool=false)

    # First, draw the number of planets in the resonant chain:
    # TODO: add to sim_param? Hardcoded for now
    res_chain_lengths = [2, 3, 4, 5, 6]
    @assert(minimum(res_chain_lengths) >= 2)

    num_pl = rand(res_chain_lengths)

    # Generate a set of periods, planet radii, and planet masses:
    Plist, Rlist, masslist = generate_resonant_chain(star, sim_param; n=num_pl)

    # Draw eccentricities and inclinations and orient the system and the orbits of the planets:
    ecclist, inclskylist, ωlist, Ωskylist, meananomlist, incl_invariable, Ω_invariable = draw_planetary_system_orbits_and_sky_orientation_by_distributing_amd(Plist, masslist, star, sim_param; verbose=verbose)

    # NOTE: currently, the resonant chains are not enforced to pass our mutual-Hill stability criterion
    # NOTE: NOT including eccentricities in final mutual-Hill stability test since they were drawn by distributing AMD after the periods were set
    #####@assert(test_stability(Plist, masslist, star.mass, sim_param))

    # Packing and returning the system in structs:
    pl = Array{Planet}(undef, num_pl)
    orbit = Array{Orbit}(undef, num_pl)
    for i in 1:num_pl
        pl[i] = Planet(Rlist[i], masslist[i], 0)
        orbit[i] = Orbit(Plist[i], ecclist[i], inclskylist[i], ωlist[i], Ωskylist[i], meananomlist[i])
    end
    sys_ref_plane = SystemPlane(incl_invariable, Ω_invariable)

    return PlanetarySystem(star, pl, orbit, sys_ref_plane)
end

function draw_system_clustered_amd_model(star::StarAbstract, sim_param::SimParam; verbose::Bool=false)

    # Decide whether to assign a planetary system to the star at all:
    f_swpa = calc_fraction_of_stars_with_planets_linear_with_color(star, sim_param)
    if rand() > f_swpa
        #println("Star not assigned a planetary system.")
        return PlanetarySystem(star)
    end

    # If get to this point, will try to draw a planetary system:
    ps = generate_planetary_system_clustered_periods_and_sizes_distribute_amd(star, sim_param; verbose=verbose)
    return ps
end

function draw_system_clustered_photoevap_amd_model(star::StarAbstract, sim_param::SimParam; verbose::Bool=false)

    # Decide whether to assign a planetary system to the star at all:
    f_swpa = calc_fraction_of_stars_with_planets_linear_with_color(star, sim_param)
    if rand() > f_swpa
        #println("Star not assigned a planetary system.")
        return PlanetarySystem(star)
    end

    # If get to this point, will try to draw a planetary system:
    ps = generate_planetary_system_clustered_periods_and_sizes_photoevap_distribute_amd(star, sim_param; verbose=verbose)
    return ps
end

function draw_system_resonant_chain_amd_model(star::StarAbstract, sim_param::SimParam; verbose::Bool=false)

    # Decide whether to assign a planetary system to the star at all:
    f_swpa = calc_fraction_of_stars_with_planets_linear_with_color(star, sim_param)
    if rand() > f_swpa
        #println("Star not assigned a planetary system.")
        return PlanetarySystem(star)
    end

    # If get to this point, will try to draw a planetary system:
    ps = generate_planetary_system_resonant_chain_clustered_sizes_distribute_amd(star, sim_param; verbose=verbose)
    return ps
end

function draw_system_clustered_and_resonant_chain_amd_mixture_model(star::StarAbstract, sim_param::SimParam; verbose::Bool=false)

    # Decide whether to assign a planetary system to the star at all:
    f_swpa = calc_fraction_of_stars_with_planets_linear_with_color(star, sim_param)
    if rand() > f_swpa
        #println("Star not assigned a planetary system.")
        return PlanetarySystem(star)
    end

    # If get to this point, will try to draw a planetary system:
    f_reschain = get_real(sim_param, "f_resonant_chains")
    @assert(0 <= f_reschain <= 1)
    if rand() > f_reschain
        # Draw from the clustered periods and sizes population:
        ps = generate_planetary_system_clustered_periods_and_sizes_distribute_amd(star, sim_param; verbose=verbose)
    else
        # Draw from the resonant chain population:
        ps = generate_planetary_system_resonant_chain_clustered_sizes_distribute_amd(star, sim_param; verbose=verbose)
    end
    return ps
end

# Dictionary of model names and their functions for drawing from them:
models = Dict{String,Function}()
models["clustered_amd_model"] = draw_system_clustered_amd_model
models["clustered_photoevap_amd_model"] = draw_system_clustered_photoevap_amd_model
models["resonant_chain_amd_model"] = draw_system_resonant_chain_amd_model
models["clustered_and_resonant_chain_amd_mixture_model"] = draw_system_clustered_and_resonant_chain_amd_mixture_model

function draw_system_model(star::StarAbstract, sim_param::SimParam; verbose::Bool=false)
    # Wrapper function that reads the model name from sim_param and draws from the relevant model
    if haskey(sim_param, "max_attempts_cond") # indicates that we want to condition on a given planet
        return draw_system_model_conditional(star, sim_param; verbose=verbose)
    end

    model_name = get(sim_param, "model_name", "")
    return models[model_name](star, sim_param; verbose=verbose)
end
