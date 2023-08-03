using PyPlot

include("../src/models_test.jl")



##### To make some plots testing some new functions for implementing Neil & Rogers (2020) model ("NR20"):

logrange(x1,x2,n) = (10^y for y in range(log10(x1), log10(x2), length=n)) # useful convenience function

M_array = logrange(0.1, 1e4, 1000) # array of planet masses (Earth masses) given the bounds in NR20



### Testing the functions for computing the mean and scatter in planet radius given planet mass:

# Best-fit parameters (medians of posteriors) for Model 2 from Table 1 of NR20:
μ_M, σ_M = 1.00, 1.65 # ln(Earth masses)
C = 2.37 # Earth radii
M_break1 = 17.4 # Earth masses
M_break2 = 175.7 # Earth masses
γ0, γ1, γ2 = 0.00, 0.74, 0.04
σ0, σ1, σ2 = 0.18, 0.34, 0.10

μσ_R_array = [mean_radius_and_scatter_given_mass_neil_rogers2020(M; C=C, M_break1=M_break1, M_break2=M_break2, γ0=γ0, γ1=γ1, γ2=γ2, σ0=σ0, σ1=σ1, σ2=σ2) for M in M_array]

μ_R_array = [μσ_R[1] for μσ_R in μσ_R_array]
σ_R_array = [μσ_R[2] for μσ_R in μσ_R_array]

# Plot planet radius (mean and scatter) vs. planet mass:
plot = figure(figsize=(8,8))
subplot(1,1,1)
fill_between([log10(M) for M in M_array], log10.(μ_R_array .+ μ_R_array .* σ_R_array), log10.(μ_R_array .- μ_R_array .* σ_R_array), label="NR20, Model 2, scatter")
scatter([log10(M) for M in M_array], log10.(μ_R_array), s=10, label="NR20, Model 2, mean")
scatter([log10(M) for M in M_array], [log10(radius_given_mass_pure_iron_fit_seager2007(M)) for M in M_array], s=10, label="S07, pure-iron")
scatter([log10(M) for M in M_array], [log10(radius_given_mass_pure_silicate_fit_seager2007(M)) for M in M_array], s=10, label="S07, pure-silicate")
xlabel("log10(M [M_earth])", fontsize=20)
ylabel("log10(R [R_earth])", fontsize=20)
legend(fontsize=16)
tight_layout()



### Testing the functions for computing envelope mass given planet mass:

M_env_array = [envelope_mass_smoothed_low_high_neil_rogers2020(M) for M in M_array]

# Plot planet envelope mass vs. planet mass:
plot = figure(figsize=(8,8))
subplot(1,1,1)
scatter([log10(M) for M in M_array], log10.(M_env_array), s=10, label="NR20")
xlabel("log10(M [M_earth])", fontsize=20)
ylabel("log10(M_env [M_earth])", fontsize=20)
legend(fontsize=16)
tight_layout()

# Plot fraction of planet mass in envelope vs. planet mass:
plot = figure(figsize=(8,8))
subplot(1,1,1)
scatter([log10(M) for M in M_array], M_env_array./M_array, s=10, label="NR20")
xlabel("log10(M [M_earth])", fontsize=20)
ylabel("M_env / M", fontsize=20)
legend(fontsize=16)
tight_layout()





### To draw a simple population from the distributions to test the existing models:

R_min = 0.4 # minimum planet radii to truncate at, in Earth radii

Mdist = LogNormal(μ_M, σ_M) # the distribution of initial planet masses

# Draw a simple population:
N_pl = 10000 # number of planets to draw

M_init_all = rand(Mdist, N_pl) # initial planet masses (Earth masses)
R_init_all = [draw_radius_given_mass_neil_rogers2020(M; C=C, M_break1=M_break1, M_break2=M_break2, γ0=γ0, γ1=γ1, γ2=γ2, σ0=σ0, σ1=σ1, σ2=σ2) for M in M_init_all] # initial planet radii (Earth radii)
R_init_all[R_init_all .< R_min] .= R_min
M_env_all = map(M -> envelope_mass_smoothed_low_high_neil_rogers2020(M), M_init_all) # initial envelope masses (Earth masses)
