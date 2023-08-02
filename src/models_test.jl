##### Temporary file for coding up some functions for computing models from Neil & Rogers (2020)
##### TODO: move to other files/organize these functions
##### TODO: write better docstrings
##### TODO: turn some parameters into sim params

#using Distributions
using ExoplanetsSysSim



function mean_radius_and_scatter_given_mass_neil_rogers2020(M::Real; M_break1::Real, M_break2::Real, γ0::Real, γ1::Real, γ2::Real, σ0::Real, σ1::Real, σ2::Real)

    # M: planet mass (Earth masses)
    # M_break1: break point 1 in planet mass (between low and intermediate masses)
    # M_break2: break point 2 in planet mass (between intermediate and high masses)
    # γ0: power-law slope for low mass range
    # γ1: power-law slope for intermediate mass range
    # γ2: power-law slope for high mass range
    # σ0: scatter in radius for low mass range
    # σ1: scatter in radius for intermediate mass range
    # σ2: scatter in radius for high mass range
    
    # Power-laws in each mass regime (low/intermediate/high):
    C = 1. # not sure what this normalization should be yet
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


function draw_radius_given_mass_neil_rogers2020(M::Real, μ::Real, σ::Real)
    Rdist = Normal(μ, σ*μ) # normal distribution with a fractional scatter
    return rand(Rdist)
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

function bolometric_flux_at_planet_period(P::Real, Mstar::Real, L::Real)

    # P: planet orbital period (days)
    # Mstar: stellar mass (Solar masses)
    # L: stellar luminosity (erg s^-1)
    
    sma = ExoplanetsSysSim.semimajor_axis(P, Mstar) * ExoplanetsSysSim.AU_in_m_IAU2012 * 100. # cm
    F = L/(4π*sma^2) # erg s^-1 cm^-2
    return F
end

F_oplus = bolometric_flux_at_planet_period(365.25, 1., luminosity_star(1., 5780.)) # bolometric flux on Earth at present (erg s^-1 cm^-2)

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
