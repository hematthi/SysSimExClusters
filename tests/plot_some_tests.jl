using PyPlot

include("../src/models_test.jl")



##### To make some plots testing some new functions for implementing Neil & Rogers (2020) model ("NR20"):

logrange(x1,x2,n) = (10^y for y in range(log10(x1), log10(x2), length=n)) # useful convenience function

M_min, M_max = 0.1, 1e4
M_array = logrange(M_min, M_max, 1000) # array of planet masses (Earth masses) given the bounds in NR20



### Best-fit parameters (medians of posteriors) for Model 2 from Table 1 of NR20:
μ_M, σ_M = 1.00, 1.65 # ln(Earth masses)
C = 2.37 # Earth radii
M_break1 = 17.4 # Earth masses
M_break2 = 175.7 # Earth masses
γ0, γ1, γ2 = 0.00, 0.74, 0.04
σ0, σ1, σ2 = 0.18, 0.34, 0.10
R_min, R_max = 0.4, 20. # minimum and maximum planet radii to truncate at (Earth radii)
P_min, P_break, P_max = 0.3, 7.16, 100. # days
β1, β2 = 0.88, -0.76
α = 7.98 # fudge factor for mass-loss timescale



### Testing the functions for computing the mean and scatter in planet radius given planet mass:

μσ_R_array = [mean_radius_and_scatter_given_mass_neil_rogers2020(M; C=C, M_break1=M_break1, M_break2=M_break2, γ0=γ0, γ1=γ1, γ2=γ2, σ0=σ0, σ1=σ1, σ2=σ2) for M in M_array]

μ_R_array = [μσ_R[1] for μσ_R in μσ_R_array]
σ_R_array = [μσ_R[2] for μσ_R in μσ_R_array]

# Plot planet radius (mean and scatter) vs. planet mass:
fig = figure(figsize=(8,8))
subplot(1,1,1)
fill_between([log10(M) for M in M_array], log10.(μ_R_array .+ μ_R_array .* σ_R_array), log10.(μ_R_array .- μ_R_array .* σ_R_array), alpha=0.2, label="NR20, Model 2, scatter")
plot([log10(M) for M in M_array], log10.(μ_R_array), label="NR20, Model 2, mean")
plot([log10(M) for M in M_array], [log10(radius_given_mass_pure_iron_fit_seager2007(M)) for M in M_array], label="S07, pure-iron")
plot([log10(M) for M in M_array], [log10(radius_given_mass_pure_silicate_fit_seager2007(M)) for M in M_array], label="S07, pure-silicate")
xlabel("log10(M [M_earth])", fontsize=20)
ylabel("log10(R [R_earth])", fontsize=20)
legend(fontsize=16)
tight_layout()



### Testing the functions for computing envelope mass given planet mass:

M_env_array = [envelope_mass_smoothed_low_high_neil_rogers2020(M) for M in M_array]

# Plot planet envelope mass vs. planet mass:
fig = figure(figsize=(8,8))
subplot(1,1,1)
plot([log10(M) for M in M_array], log10.(M_env_array), label="NR20")
xlabel("log10(M [M_earth])", fontsize=20)
ylabel("log10(M_env [M_earth])", fontsize=20)
legend(fontsize=16)
tight_layout()

# Plot fraction of planet mass in envelope vs. planet mass:
fig = figure(figsize=(8,8))
subplot(1,1,1)
plot([log10(M) for M in M_array], M_env_array./M_array, label="NR20")
xlabel("log10(M [M_earth])", fontsize=20)
ylabel("M_env / M", fontsize=20)
legend(fontsize=16)
tight_layout()





### To draw a simple population from the distributions to test the existing models:

M_init_dist = LogNormal(μ_M, σ_M) # the distribution of initial planet masses

# Draw a simple population:
N_pl = 10000 # number of planets to draw

M_init_all = rand(M_init_dist, N_pl) # initial planet masses (Earth masses)
R_init_all = map(M -> draw_radius_given_mass_neil_rogers2020(M; C=C, M_break1=M_break1, M_break2=M_break2, γ0=γ0, γ1=γ1, γ2=γ2, σ0=σ0, σ1=σ1, σ2=σ2), M_init_all) # initial planet radii (Earth radii)
R_init_all[R_init_all .< R_min] .= R_min
M_env_all = map(M -> envelope_mass_smoothed_low_high_neil_rogers2020(M), M_init_all) # initial envelope masses (Earth masses)

P_all = ExoplanetsSysSim.draw_broken_power_law(β1, β2, P_min, P_max, P_break, N_pl) # orbital periods (days)
F_p_all = map(P -> bolometric_flux_at_planet_period(P), P_all) # assuming all planets are around Solar mass and luminosity stars (erg s^-1 cm^-2)

t_loss_all = map(i -> mass_loss_timescale_lopez2012(M_env_all[i], R_init_all[i], F_p_all[i]), 1:N_pl) # mass-loss timescales (Gyrs)
p_ret_all = map(t -> prob_retain_envelope_neil_rogers2020(5., t, α), t_loss_all) # probabilities of retaining envelope
bools_ret_all = rand(N_pl) .< p_ret_all # 1 = retain envelope, 0 = lose envelope

M_final_all = deepcopy(M_init_all)
M_final_all[.!bools_ret_all] = M_final_all[.!bools_ret_all] .- M_env_all[.!bools_ret_all] # final planet masses (Earth masses)

R_final_all = deepcopy(R_init_all)
R_final_all[.!bools_ret_all] = map(M -> draw_radius_given_mass_neil_rogers2020(M, radius_given_mass_pure_silicate_fit_seager2007(M), 0.05), M_final_all[.!bools_ret_all]) # final planet radii (Earth radii)



### Plot the population:

# Plot initial planet radius vs. mass vs. envelope mass fraction:
fig = figure(figsize=(10,8))
subplot(1,1,1)
fill_between(log10.(M_array), log10.(μ_R_array .+ μ_R_array .* σ_R_array), log10.(μ_R_array .- μ_R_array .* σ_R_array), alpha=0.2, label="NR20, Model 2, scatter")
plot(log10.(M_array), log10.(μ_R_array), label="NR20, Model 2, mean")
plot(log10.(M_array), [log10(radius_given_mass_pure_iron_fit_seager2007(M)) for M in M_array], label="S07, pure-iron")
plot(log10.(M_array), [log10(radius_given_mass_pure_silicate_fit_seager2007(M)) for M in M_array], label="S07, pure-silicate")
sc = scatter(log10.(M_init_all), log10.(R_init_all), s=1, c=M_env_all./M_init_all, label="Sample population") #c=log10.(M_env_all)
cbar = colorbar(sc)
cbar.set_label("M_env_init / M_init")
xlim([log10(M_min), log10(M_max)])
xlabel("log10(M_init [M_earth])", fontsize=20)
ylabel("log10(R_init [R_earth])", fontsize=20)
legend(fontsize=16)
tight_layout()

# Plot distributions of initial planet radius, mass, and envelope mass fraction:
fig = figure(figsize=(16,5))
subplot(1,3,1)
hist(log10.(R_init_all), bins=range(log10(R_min), log10(R_max), 51))
xlabel("log10(R_init [R_earth])", fontsize=20)
ylabel("Count", fontsize=20)
subplot(1,3,2)
hist(log10.(M_init_all), bins=range(log10(M_min), log10(M_max), 51))
xlabel("log10(M_init [M_earth])", fontsize=20)
subplot(1,3,3)
hist(log10.(M_env_all./M_init_all), bins=range(-1., 0., 51))
xlabel("log10(M_env_init/M_init [M_earth])", fontsize=20)
tight_layout()

# Plot distributions of orbital period, mass-loss timescale, and probability of envelope retention:
fig = figure(figsize=(16,5))
subplot(1,3,1)
hist(log10.(P_all), bins=range(log10(P_min), log10(P_max), 51))
xlabel("log10(P [days])", fontsize=20)
ylabel("Count", fontsize=20)
subplot(1,3,2)
hist(log10.(t_loss_all), bins=range(-6, 6, 51))
xlabel("log10(t_loss [Gyr])", fontsize=20)
subplot(1,3,3)
hist(log10.(p_ret_all), bins=range(-6, 0, 51))
xlabel("log10(Prob_ret)", fontsize=20)
tight_layout()

# Plot period vs. initial radius vs. mass-loss timescale:
fig = figure(figsize=(10,8))
subplot(1,1,1)
sc = scatter(log10.(P_all), log10.(R_init_all), s=1, c=log10.(t_loss_all))
cbar = colorbar(sc)
cbar.set_label("log10(t_loss [Gyr])")
xlim([log10(P_min), log10(P_max)])
xlabel("log10(P [days])", fontsize=20)
ylabel("log10(R_init [R_earth])", fontsize=20)
tight_layout()

# Plot period vs. initial envelope mass vs. mass-loss timescale:
fig = figure(figsize=(10,8))
subplot(1,1,1)
sc = scatter(log10.(P_all), log10.(M_env_all), s=1, c=log10.(t_loss_all))
cbar = colorbar(sc)
cbar.set_label("log10(t_loss [Gyr])")
xlim([log10(P_min), log10(P_max)])
xlabel("log10(P [days])", fontsize=20)
ylabel("log10(M_env_init [M_earth])", fontsize=20)
tight_layout()

#

# Plot distributions of final planet radius, mass, and envelope mass fraction (also overplotted with initial distributions for comparison):
fig = figure(figsize=(16,5))
subplot(1,3,1)
hist(log10.(R_init_all), bins=range(log10(R_min), log10(R_max), 51), alpha=0.2, label="Initial")
hist(log10.(R_final_all), bins=range(log10(R_min), log10(R_max), 51), alpha=0.2, label="Final")
xlabel("log10(R_final [R_earth])", fontsize=20)
ylabel("Count", fontsize=20)
legend(fontsize=16)
subplot(1,3,2)
hist(log10.(M_init_all), bins=range(log10(M_min), log10(M_max), 51), alpha=0.2, label="Initial")
hist(log10.(M_final_all), bins=range(log10(M_min), log10(M_max), 51), alpha=0.2, label="Final")
xlabel("log10(M_final [M_earth])", fontsize=20)
subplot(1,3,3)
hist(log10.(M_env_all./M_init_all), bins=range(-1., 0., 51), alpha=0.2, label="Initial")
hist(log10.(M_env_all./M_final_all)[bools_ret_all], bins=range(-1., 0., 51), alpha=0.2, label="Final")
xlabel("log10(M_env_init/M_final [M_earth])", fontsize=20)
tight_layout()

# Plot period vs. final radius (split into evaporated/retained envelopes):
fig = figure(figsize=(10,8))
subplot(1,1,1)
scatter(log10.(P_all[.!bools_ret_all]), log10.(R_final_all[.!bools_ret_all]), s=1, color="r", label="Lost envelope")
scatter(log10.(P_all[bools_ret_all]), log10.(R_final_all[bools_ret_all]), s=1, color="b", label="Retained envelope")
xlim([log10(P_min), log10(P_max)])
xlabel("log10(P [days])", fontsize=20)
ylabel("log10(R_final [R_earth])", fontsize=20)
legend(fontsize=16)
tight_layout()

# Plot final planet radius vs. mass (split into evaporated/retained envelopes):
fig = figure(figsize=(10,8))
subplot(1,1,1)
fill_between(log10.(M_array), log10.(μ_R_array .+ μ_R_array .* σ_R_array), log10.(μ_R_array .- μ_R_array .* σ_R_array), alpha=0.2, label="NR20, Model 2, scatter")
plot(log10.(M_array), log10.(μ_R_array), label="NR20, Model 2, mean")
plot(log10.(M_array), [log10(radius_given_mass_pure_iron_fit_seager2007(M)) for M in M_array], label="S07, pure-iron")
plot(log10.(M_array), [log10(radius_given_mass_pure_silicate_fit_seager2007(M)) for M in M_array], label="S07, pure-silicate")
scatter(log10.(M_final_all[.!bools_ret_all]), log10.(R_final_all[.!bools_ret_all]), s=1, color="r", label="Lost envelope")
scatter(log10.(M_final_all[bools_ret_all]), log10.(R_final_all[bools_ret_all]), s=1, color="b", label="Retained envelope")
xlim([log10(M_min), log10(M_max)])
xlabel("log10(M_final [M_earth])", fontsize=20)
ylabel("log10(R_final [R_earth])", fontsize=20)
legend(fontsize=16)
tight_layout()

# Plot initial vs. final planet mass vs. initial envelope mass:
fig = figure(figsize=(10,8))
subplot(1,1,1)
sc = scatter(log10.(M_init_all), log10.(M_final_all), s=1, c=log10.(M_env_all))
cbar = colorbar(sc)
cbar.set_label("log10(M_env_init [M_earth])")
xlim([log10(M_min), log10(M_max)])
ylim([log10(M_min), log10(M_max)])
xlabel("log10(M_init [M_earth])", fontsize=20)
ylabel("log10(M_final [M_earth])", fontsize=20)
tight_layout()

# Plot initial vs. final planet radius vs. initial envelope mass:
fig = figure(figsize=(10,8))
subplot(1,1,1)
sc = scatter(log10.(R_init_all), log10.(R_final_all), s=1, c=log10.(M_env_all))
cbar = colorbar(sc)
cbar.set_label("log10(M_env_init [M_earth])")
xlim([log10(R_min), log10(R_max)])
ylim([log10(R_min), log10(R_max)])
xlabel("log10(R_init [R_earth])", fontsize=20)
ylabel("log10(R_final [R_earth])", fontsize=20)
tight_layout()
