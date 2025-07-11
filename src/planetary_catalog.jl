using JLD2
using FileIO
using Dates

sim_param = setup_sim_param_model()

include("misc_functions.jl")
include("complexity_stats.jl")





##### To load the DR25 KOI catalog, update the planet radii using Gaia stellar radii, apply the necessary cuts:

planet_catalog = CSV.read(joinpath(dirname(pathof(ExoplanetsSysSim)), "../data/q1_q17_dr25_koi.csv"), DataFrame, header=157)
#planet_catalog = load(joinpath(dirname(pathof(ExoplanetsSysSim)), "../data/q1_q17_dr25_koi.jld2"), "koi")

stellar_catalog = ExoplanetsSysSim.StellarTable.setup_star_table(sim_param)



"""
    keep_planet_candidates_given_sim_param(planet_catalog; sim_param, stellar_catalog, recompute_radii=true)

Reduce a planet catalog to a cleaned sample of planets ("confirmed" and "candidate") around the stars in the stellar catalog, with periods and radii in the range specified in `sim_param`. Optionally recompute the planet radii using the revised stellar radii (assuming `stellar_catalog` has revised stellar radii from Gaia DR2).

# Arguments:
- `planet_catalog::DataFrame`: KOI table.
- `sim_param::SimParam`: a SimParam object containing various simulation parameters.
- `stellar_catalog::DataFrame`: table of stars.
- `recompute_radii::Bool=true`: whether to recompute the planet radii by using the transit depths in `planet_catalog` and the stellar radii in `stellar_catalog` (default true). If true, also replaces the stellar radii and masses in `planet_catalog` with the values in `stellar_catalog`.

# Returns:
- `planets_keep::DataFrame`: table of remaining planets.
"""
function keep_planet_candidates_given_sim_param(planet_catalog::DataFrame; sim_param::SimParam, stellar_catalog::DataFrame, recompute_radii::Bool=true)

    planets_keep = planet_catalog[(planet_catalog[!,:koi_disposition] .== "CONFIRMED") .| (planet_catalog[!,:koi_disposition] .== "CANDIDATE"), :] # table containing only the confirmed and candidate objects
    println("Candidate and confirmed planets: ", size(planets_keep, 1))

    in_stellar_catalog = [] # will be filled with booleans indicating whether each koi in 'planets_keep' is found in the 'stellar_catalog' or not
    for i in 1:length(planets_keep[!,:kepid])
        if any(x->x==planets_keep[i,:kepid], stellar_catalog[!,:kepid])
            push!(in_stellar_catalog, true)
        else
            push!(in_stellar_catalog, false)
        end
    end
    in_stellar_catalog_indices = findall(in_stellar_catalog)

    planets_keep = planets_keep[in_stellar_catalog_indices, :]
    println("After removing planets not around stars in stellar catalog: ", size(planets_keep, 1))

    planets_keep = planets_keep[(ismissing.(planets_keep[!,:koi_duration]) .== false) .& (planets_keep[!,:koi_duration] .>= 0), :]
    planets_keep = planets_keep[(ismissing.(planets_keep[!,:koi_depth]) .== false) .& (planets_keep[!,:koi_depth] .> 0), :]
    println("After removing planets with missing or negative transit durations or depths: ", size(planets_keep, 1))

    if recompute_radii
        for (i,kepid) in enumerate(planets_keep[!,:kepid]) # to replace the stellar and planetary radii in 'planets_keep' with the more reliable values as derived from the stellar properties in 'stellar_catalog'
            stellar_radii_new = stellar_catalog[stellar_catalog[!,:kepid] .== kepid, :radius][1]
            stellar_mass = stellar_catalog[stellar_catalog[!,:kepid] .== kepid, :mass][1]
            planets_keep[i, :koi_srad] = stellar_radii_new
            planets_keep[i, :koi_smass] = stellar_mass # to save stellar masses because for some reason, currently the stellar mass for planets after the first planet in multi-planet systems are set to nan
            planets_keep[i, :koi_prad] = (1 ./ExoplanetsSysSim.earth_radius)*stellar_radii_new*sqrt(planets_keep[i, :koi_depth]/(1e6))
        end
    end

    # To make additional cuts in periods and radii to be comparable to our simulated sample:
    planets_keep = planets_keep[(planets_keep[!,:koi_period] .> get_real(sim_param,"min_period")) .& (planets_keep[!,:koi_period] .< get_real(sim_param,"max_period")), :]
    planets_keep = planets_keep[(planets_keep[!,:koi_prad] .> get_real(sim_param,"min_radius")) .& (planets_keep[!,:koi_prad] .< get_real(sim_param,"max_radius")) .& (.~ismissing.(planets_keep[!,:koi_prad])), :]
    planets_keep[!,:koi_period] = collect(skipmissing(planets_keep[!,:koi_period]))
    planets_keep[!,:koi_prad] = collect(skipmissing(planets_keep[!,:koi_prad]))
    println("After applying our period and radius cuts (final count): ", size(planets_keep, 1))

    return planets_keep
end

@time planets_cleaned = keep_planet_candidates_given_sim_param(planet_catalog; sim_param=sim_param, stellar_catalog=stellar_catalog, recompute_radii=true)

# If we want to write the cleaned planetary catalog and the stellar catalog to a csv file, keeping only the columns we need:
#CSV.write("q1_q17_dr25_gaia_berger_fgk_HFR2020b_koi_cleaned.csv", planets_cleaned[[:kepid, :kepoi_name, :koi_disposition, :koi_pdisposition, :koi_score, :koi_period, :koi_duration, :koi_depth, :koi_prad, :koi_steff, :koi_slogg, :koi_srad, :koi_smass]])
#CSV.write("q1_q17_dr25_gaia_berger_fgk_HFR2020b_cleaned.csv", stellar_catalog[[:kepid, :mass, :radius, :teff, :bp_rp, :lum_val, :e_bp_min_rp_interp, :e_bp_min_rp_val, :rrmscdpp04p5]])





##### To compute arrays of the observables (multiplicities, periods, period ratios, transit durations, transit depths, period-normalized transit duration ratios (xi), and transit depth ratios) from the remaining sample of planets:

"""
    calc_summary_stats_Kepler(stellar_catalog, planets_cleaned)

Compute the summary statistics of a Kepler planet catalog and compile them into dictionary.

# Arguments:
- `stellar_catalog::DataFrame`: table of target stars.
- `planets_cleaned::DataFrame`: table of planet candidates around the stars in `stellar_catalog`.
- `sim_param::SimParam`: a SimParam object containing various simulation parameters.

# Returns:
A `CatalogSummaryStatistics` object, where the `stat` field is a dictionary containing the summary statistics.
"""
function calc_summary_stats_Kepler(stellar_catalog::DataFrame, planets_cleaned::DataFrame; sim_param::SimParam)
    KOI_systems = [x[1:6] for x in planets_cleaned[!,:kepoi_name]]
    checked_bools = zeros(size(planets_cleaned,1)) # 0's denote KOI that were not checked yet; 1's denote already checked KOI

    M_obs = Int64[] # list for the planet multiplicities of the systems
    period_ratios = Float64[] # list for the period ratios
    duration_ratios = Float64[] # list for the period-normalized transit duration ratios
    duration_ratios_nonmmr = Float64[] # list for the period-normalized transit duration ratios not near any resonances
    duration_ratios_mmr = Float64[] # list for the period-normalized transit duration ratios near a resonance
    depth_ratios = Float64[] # list for the transit depth ratios
    periods = collect(skipmissing(planets_cleaned[!,:koi_period])) # list of the periods (days)
    durations = collect(skipmissing(planets_cleaned[!,:koi_duration]./24)) # list of the transit durations (days)
    durations_norm_circ = Float64[] # list for the transit durations normalized by the circular, central durations
    durations_norm_circ_singles = Float64[]
    durations_norm_circ_multis = Float64[]
    depths = collect(skipmissing(planets_cleaned[!,:koi_depth]./(1e6))) # list of the transit depths (fraction)
    radii = collect(skipmissing(planets_cleaned[!,:koi_prad])) # list for the planet radii (Earth radii)
    radii_delta_from_valley_list = Float64[] # list for the planet radii differences from the location of the radius valley (Earth radii)
    
    # Parameters for computing the difference in planet radii from the location of the radius valley:
    radius_valley_slope = get_real(sim_param, "radius_valley_slope")
    radius_valley_offset = get_real(sim_param, "radius_valley_offset")
    radius_valley_min_radius_include = get_real(sim_param, "radius_valley_min_radius_include")
    radius_valley_max_radius_include = get_real(sim_param, "radius_valley_max_radius_include")
    radius_valley_min_period_include = get_real(sim_param, "radius_valley_min_period_include")
    radius_valley_max_period_include = get_real(sim_param, "radius_valley_max_period_include")

    depths_above = Float64[] # list for the transit depths of planets above the photoevaporation boundary in Carrera et al 2018
    depths_below = Float64[] # list for the transit depths of planets below the boundary
    depth_ratios_above = Float64[] # list for the transit depth ratios, both above the boundary
    depth_ratios_below = Float64[] # list for the transit depth ratios, both below the boundary
    depth_ratios_across = Float64[] # list tfor the transit depth ratios, across the boundary

    # System level metrics taken from or inspired by Gilbert & Fabrycky (2020):
    radii_partitioning = Float64[]
    radii_monotonicity = Float64[]
    gap_complexity = Float64[]

    for i in 1:length(KOI_systems)
        if checked_bools[i] == 0 # if the KOI has not been checked (included while looking at another planet in the same system)
            system_i = (1:length(KOI_systems))[KOI_systems .== KOI_systems[i]]
            checked_bools[system_i] .= 1

            star_radius = planets_cleaned[!,:koi_srad][system_i][1]
            star_mass = planets_cleaned[!,:koi_smass][system_i][1]

            # To get the periods and transit durations in this system:
            system_P = planets_cleaned[!,:koi_period][system_i] # periods of all the planets in this system
            system_dur = planets_cleaned[!,:koi_duration][system_i] ./ 24 # transit durations (in days) of all the planets in this system
            system_dep = planets_cleaned[!,:koi_depth][system_i] # transit depths (in ppm) of all the planets in this system
            system_radii = planets_cleaned[!,:koi_prad][system_i] # radii of all the planets in this system

            system_sort_i = sortperm(system_P) # indices that would sort the periods of the planets in this system
            system_P = system_P[system_sort_i] # periods of all the planets in this system, sorted
            system_dur = system_dur[system_sort_i] # transit durations of all the planets in this system, sorted by period
            system_dep = system_dep[system_sort_i]/(1e6) # transit depths of all the planets in this system, sorted by period
            system_radii = system_radii[system_sort_i] # radii of all the planets in this system, sorted by period

            system_dur_norm_circ = system_dur ./ map(P -> calc_transit_duration_central_circ_obs(P; Mstar=star_mass, Rstar=star_radius), system_P)
            append!(durations_norm_circ, system_dur_norm_circ)
            if length(system_dur_norm_circ) == 1
                append!(durations_norm_circ_singles, system_dur_norm_circ)
            elseif length(system_dur_norm_circ) > 1
                append!(durations_norm_circ_multis, system_dur_norm_circ)
            end

            # To count the total number of planets in this system:
            push!(M_obs, length(system_P))

            # To compute the period ratios, period-normalized transit duration ratios, and transit depth ratios in this system:
            system_Pratio = system_P[2:end] ./ system_P[1:end-1] # period ratios of all the adjacent planet pairs in this system
            system_depratio = system_dep[2:end] ./ system_dep[1:end-1] # transit depth ratios of all the adjacent planet pairs in this system
            system_xi = (system_dur[1:end-1] ./ system_dur[2:end]) .* (system_Pratio .^(1//3)) # period-normalized transit duration ratios of all the adjacent planet pairs in this system

            append!(period_ratios, system_Pratio)
            append!(depth_ratios, system_depratio)
            append!(duration_ratios, system_xi)
            for (j,period_ratio) in enumerate(system_Pratio)
                if is_period_ratio_near_resonance(period_ratio, sim_param)
                    append!(duration_ratios_mmr, system_xi[j])
                else
                    append!(duration_ratios_nonmmr, system_xi[j])
                end
            end
            
            # To compute the radii delta to the radius valley, if the planet is within bounds:
            for j in 1:length(system_radii)
                period = system_P[j]
                radius = system_radii[j]
                if (radius_valley_min_radius_include <= radius <= radius_valley_max_radius_include) && (radius_valley_min_period_include <= period <= radius_valley_max_period_include)
                    radius_delta_valley = radius_delta_from_period_radius_gap(radius, period; m=radius_valley_slope, Rgap0=radius_valley_offset)
                    append!(radii_delta_from_valley_list, radius_delta_valley)
                end
            end

            # To separate the planets in the system as above and below the boundary:
            system_above_bools = [photoevap_boundary_Carrera2018(system_radii[x], system_P[x]) for x in 1:length(system_P)]
            #if length(system_above_bools) > 1 println(system_above_bools) end

            # To record the transit depths of the planets above and below the boundary:
            for (j,dep) in enumerate(system_dep)
                if system_above_bools[j] == 1
                    append!(depths_above, dep)
                elseif system_above_bools[j] == 0
                    append!(depths_below, dep)
                end
            end

            # To record the transit depth ratios of the planets above, below, and across the boundary:
            for (j,depratio) in enumerate(system_depratio)
                if system_above_bools[j] + system_above_bools[j+1] == 2 # both planets are above the boundary
                    append!(depth_ratios_above, depratio)
                elseif system_above_bools[j] + system_above_bools[j+1] == 1 # one planet is above, the other is below the boundary
                    append!(depth_ratios_across, depratio)
                elseif system_above_bools[j] + system_above_bools[j+1] == 0 # both planets are below the boundary
                    append!(depth_ratios_below, depratio)
                end
            end

            # To compute system level metrics taken from or inspired by Gilbert & Fabrycky (2020):
            if length(system_radii) >= 2
                append!(radii_partitioning, partitioning(system_radii))
                append!(radii_monotonicity, monotonicity_GF2020(system_radii))
            end
            if length(system_P) >= 3
                append!(gap_complexity, gap_complexity_GF2020(system_P))
            end
        end
    end

    Nmult = [sum(M_obs .== k) for k in 1:maximum(M_obs)]

    # To create a dictionary for the Kepler data:
    stat = Dict{String,Any}()
    stat["num_targets"] = size(stellar_catalog,1)
    stat["num_tranets"] = size(planets_cleaned,1)
    stat["num_n-tranet_systems"] = Nmult
    stat["periods"] = periods
    stat["period_ratios"] = period_ratios
    stat["durations"] = durations
    stat["durations_norm_circ"] = durations_norm_circ
    stat["durations_norm_circ_singles"] = durations_norm_circ_singles
    stat["durations_norm_circ_multis"] = durations_norm_circ_multis
    stat["duration_ratios"] = duration_ratios
    stat["duration_ratios_nonmmr"] = duration_ratios_nonmmr
    stat["duration_ratios_mmr"] = duration_ratios_mmr
    stat["depths"] = depths
    stat["depths_above"] = depths_above
    stat["depths_below"] = depths_below
    stat["radii"] = radii
    stat["radii_delta_valley"] = radii_delta_from_valley_list
    stat["radius_ratios"] = sqrt.(depth_ratios)
    stat["radius_ratios_above"] = sqrt.(depth_ratios_above)
    stat["radius_ratios_below"] = sqrt.(depth_ratios_below)
    stat["radius_ratios_across"] = sqrt.(depth_ratios_across)
    stat["radii_partitioning"] = radii_partitioning
    stat["radii_monotonicity"] = radii_monotonicity
    stat["gap_complexity"] = gap_complexity

    return CatalogSummaryStatistics(stat, Dict{String,Any}())
end

ssk = calc_summary_stats_Kepler(stellar_catalog, planets_cleaned; sim_param=sim_param)



"""
    calc_summary_stats_collection_Kepler(stellar_catalog, planet_catalog, names_samples, star_id_samples, sim_param)

Compute the summary statistics of a Kepler planet catalog given a number of samples (divisions in stellar catalogs) and compile them into dictionaries in a `CatalogSummaryStatisticsCollection` object.

# Arguments:
- `stellar_catalog::DataFrame`: table of stars.
- `planet_catalog::DataFrame`: table of planets.
- `names_samples::Vector{String}`: names of the samples.
- `star_id_samples::Vector{Vector{Int64}}`: list of lists of target stars in each sample.
- `sim_param::SimParam`: a SimParam object containing various simulation parameters.

# Returns:
- `cssc::CatalogSummaryStatisticsCollection`: object containing a CatalogSummaryStatistics object for each stellar sample.
"""
function calc_summary_stats_collection_Kepler(stellar_catalog::DataFrame, planet_catalog::DataFrame, names_samples::Vector{String}, star_id_samples::Vector{Vector{Int64}}, sim_param::SimParam)
    @assert length(names_samples) == length(star_id_samples)

    @time planets_cleaned = keep_planet_candidates_given_sim_param(planet_catalog; sim_param=sim_param, stellar_catalog=stellar_catalog, recompute_radii=true)

    cssc = CatalogSummaryStatisticsCollection()
    cssc.star_id_samples["all"] = collect(1:size(stellar_catalog,1))
    cssc.css_samples["all"] = calc_summary_stats_Kepler(stellar_catalog, planets_cleaned; sim_param=sim_param)

    for (i,name) in enumerate(names_samples)
        stellar_catalog_sample = stellar_catalog[star_id_samples[i],:]
        @time planets_cleaned_sample = keep_planet_candidates_given_sim_param(planet_catalog; sim_param=sim_param, stellar_catalog=stellar_catalog_sample)
        cssc.star_id_samples[name] = star_id_samples[i]
        cssc.css_samples[name] = calc_summary_stats_Kepler(stellar_catalog_sample, planets_cleaned_sample; sim_param=sim_param)
    end
    return cssc
end
