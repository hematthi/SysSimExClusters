using DataFrames, JLD2, ExoplanetsSysSim, PyPlot

data_path = joinpath(dirname(pathof(ExoplanetsSysSim)), "../data")

fname_orig = "q1q17_dr25_gaia_berger_fgk_HFR2020b.jld2"
fname_regen = "q1q17_dr25_gaia_berger_fgk_HFR2020b_regen.jld2"

data_orig = load(joinpath(data_path, fname_orig))
data_regen = load(joinpath(data_path, fname_regen))

sc_orig = data_orig["stellar_catalog"]
sc_regen = data_regen["stellar_catalog"]



##### To compare histograms and other plots of the stellar parameters:

# Basic parameters:
params_to_plot = ["mass", "radius", "dens", "logg", "teff", "feh", "bp_rp", "lum_val"]

for param in params_to_plot
    plot = figure(figsize=(8,5))
    subplot(1,1,1)
    hist([sc_orig[!,param], sc_regen[!,param]], bins=100, histtype="step", label=["Original", "Regenerated"])
    xlabel(param, fontsize=20)
    ylabel("Counts", fontsize=20)
    legend(fontsize=16)
    tight_layout()
end

# Extinction-corrected colors:
plot = figure(figsize=(8,5))
subplot(1,1,1)
hist([sc_orig[!,"bp_rp"] .- sc_orig[!,"e_bp_min_rp_interp"], sc_regen[!,"bp_rp"] .- sc_regen[!,"e_bp_min_rp_interp"]], bins=100, histtype="step", label=["Original", "Regenerated"])
#xlim([0.5, 1.7])
xlabel("bp-rp - E*(bp-rp)", fontsize=20)
ylabel("Counts", fontsize=20)
legend(fontsize=16)
tight_layout()

# Color-magnitude (HR) diagram:
plot = figure(figsize=(8,8))
subplot(1,1,1)
scatter(sc_orig[!,"bp_rp"] .- sc_orig[!,"e_bp_min_rp_interp"], log10.(sc_orig[!,"lum_val"]), s=1, zorder=2, label="Original")
scatter(sc_regen[!,"bp_rp"] .- sc_regen[!,"e_bp_min_rp_interp"], log10.(sc_regen[!,"lum_val"]), s=1, label="Regenerated")
xlabel("Bp-Rp", fontsize=20)
ylabel("log(L)", fontsize=20)
legend(fontsize=16)
tight_layout()

# Color-temperature diagram:
plot = figure(figsize=(8,8))
subplot(1,1,1)
scatter(sc_orig[!,"bp_rp"] .- sc_orig[!,"e_bp_min_rp_interp"], sc_orig[!,"teff"], s=1, zorder=2, label="Original")
scatter(sc_regen[!,"bp_rp"] .- sc_regen[!,"e_bp_min_rp_interp"], sc_regen[!,"teff"], s=1, label="Regenerated")
xlabel("Bp-Rp", fontsize=20); ylabel("T_eff (K)", fontsize=20)
legend(fontsize=16)
tight_layout()
