
using BNP_WMReg_Joint

using BSON
import Clustering.hclust, Clustering.cutree
import Distances.pairwise, Distances.Euclidean
import BayesInference.embed, BayesInference.etr, BayesInference.logsumexp
using PDMats
using StatsBase
using Statistics
using Random
using Dates

using LinearAlgebra
BLAS.set_num_threads(1)

const whichsim = ARGS[1]
const whichfit = ARGS[2]

BSON.@load "data/$(whichsim).bson" y X y_valid X_valid Y_valid ###### set directory correctly
BSON.@load "postsimB/postsim_$(whichsim)_$(whichfit).bson"  model sims mesg ###### set directory correctly

outfile = "KL/KL_BNP_WMAR_$(mesg).txt"

nprime = length(y_valid)
all( y[1:model.n] .== model.y ) || throw("Data mismatch.")

nrep = size(Y_valid, 1)
nsim = length(sims)

L = deepcopy(model.K)
H = deepcopy(model.H)

### Calculate llik ordinate for each Y_valid

using Distributions

function β_y_modify_γ(β_y::Array{T, 2}, γ::BitArray{1}) where T <: Real
    modify_indx = findall(.!(γ))
    βout = deepcopy(β_y)

    for k in modify_indx # this works even if no modifications are necessary
        βout[:,k] *= 0
    end
    return βout
end
function β_y_modify_γ(β_y::Array{T, 2}, γ::BitArray{2}) where T <: Real
    size(β_y) == size(γ) || throw("β and γ size mismatch.")    
    return β_y .* γ
end

function ldnorm(y::Float64, mu::Float64, sig2::Float64)
    -0.5 * log(2.0π * sig2) - (0.5 * (y - mu)^2 / sig2)
end

ldw_mat = ldensweight_mat(X_valid[1:nprime, 1:L], sims[1:nsim], model.state.γδc)

ldens00 = Array{Float64,2}(undef, nrep, H)
ldens0 = Array{Float64,1}(undef, nprime)
fill!(ldens0, 0.0)

for j = 1:nprime

    for ii = 1:nsim
        βγ_y = β_y_modify_γ(sims[ii][:β_y], sims[ii][:γ])

        for h = 1:H
            Ey_h = sims[ii][:μ_y][h] - sum( βγ_y[h,:] .* (X_valid[j,1:L] - sims[ii][:μ_x][h,:]) )
            for jj = 1:nrep
                # ldens00[jj,h] = ldw_mat[ii, j, h] + logpdf(Normal(Ey_h, sqrt(sims[ii][:δ_y][h])), Y_valid[jj,j])
                ldens00[jj,h] = ldw_mat[ii, j, h] + ldnorm(Y_valid[jj,j], Ey_h, sims[ii][:δ_y][h])
            end
        end

        ldens0[j] += mean( logsumexp(ldens00, 2) ) # mean over replicates jj

    end
    if j % 50 == 0
        report_file = open(outfile, "a+")
        write(report_file, "Calculating density ordinates. Obs $(j) of $(nprime) at $(Dates.now()).\n")
        close(report_file)
    end
    # println("j = ", j, " of ", nprime)
end

B_vec = ldens0 ./ float(nsim) # posterior mean

### Compute and aggregate KL
using Distributions

if whichsim == "Ricker_single"
    function ldens_true(y::Float64, x::Vector{Float64}, theta::Float64, sig::Float64)
            μ = x[2] * exp( theta - x[2] )
            return logpdf(Normal(μ, sig), y)
    end

elseif whichsim == "Ricker_single_lognormal"
    function ldens_true(y::Float64, x::Vector{Float64}, theta::Float64, sig::Float64)
            μ = log(x[2]) + theta - x[2]
            return logpdf(LogNormal(μ, sig), y)
    end

elseif whichsim == "Ricker_twolag_lognormal"
    function ldens_true(y::Float64, x::Vector{Float64}, theta::Float64, sig::Float64)
            μ = log(x[2]) + theta - x[2]
            σ = sig * x[1]
            return logpdf(LogNormal(μ, σ), y)
    end

end

ld_true = [ ldens_true(Y_valid[jj,j], X_valid[j,:], 2.6, 0.09) for jj = 1:nrep, j = 1:nprime ]

A_vec = mean(ld_true, dims=1)[1,:]
# plot(A_vec)

# using Plotly
# trace_true = scatter(Dict("x" => 1:nprime,
#     "y" => A_vec, "name" => "data",
#     "mode" => "markers", "type" => "scatter", "marker" => Dict("color" => "blue")))
# trace_mod = scatter(Dict("x" => 1:nprime,
#     "y" => B_vec, "name" => "validation",
#     "mode" => "markers", "type" => "scatter", "marker" => Dict("color" => "red")))
# plot([trace_true, trace_mod])

## sanity check
# trace = scatter(Dict("x" => X_valid[1:nprime,2],
#     "y" => mean(Y_valid[:,1:nprime], dims=1)[1,:], "name" => "data",
#     "mode" => "markers", "type" => "scatter", "marker" => Dict("color" => "blue")))
# trace = scatter(Dict("x" => X_valid[1:nprime,2],
#     "y" => Y_valid[1,1:nprime], "name" => "data",
#     "mode" => "markers", "type" => "scatter", "marker" => Dict("color" => "blue")))
# plot([trace])


# plot(A_vec - B_vec)

KL_valid = A_vec .- B_vec

KL_out = mean(KL_valid)

report_file = open(outfile, "a+")

### Break out over certain regions by indexing KL_valid
if whichsim == "Ricker_single_lognormal"

    x2_cutoff = 2.0

    low_x2_bool = X_valid[1:nprime, 2] .<= x2_cutoff
    low_x2fit_bool = X[1:model.n, 2] .<= x2_cutoff

    n_low2 = sum(low_x2_bool)
    n_high2 = sum(.!low_x2_bool)

    n_low2_fit = sum(low_x2fit_bool)
    n_high2_fit = sum(.!low_x2fit_bool)

    write(report_file, "\n\nX2 below $(x2_cutoff), n_fit: $(n_low2_fit), n_valid: $(n_low2), mean KL: $(mean(KL_valid[low_x2_bool]))\n")
    write(report_file, "\nX2 above $(x2_cutoff), n_fit: $(n_high2_fit), n_valid: $(n_high2), mean KL: $(mean(KL_valid[.!low_x2_bool]))\n")

elseif whichsim == "Ricker_twolag_lognormal"

    x1_cutoff = 2.0
    x2_cutoff = 3.0

    low_x1_bool = X_valid[1:nprime, 1] .<= x1_cutoff
    low_x2_bool = X_valid[1:nprime, 2] .<= x2_cutoff

    low_x1fit_bool = X[1:model.n, 1] .<= x1_cutoff
    low_x2fit_bool = X[1:model.n, 2] .<= x2_cutoff


    low1_low2_bool = [ low_x1_bool[j] && low_x2_bool[j] for j = 1:nprime ] 
    low1_high2_bool = [ low_x1_bool[j] && !low_x2_bool[j] for j = 1:nprime ] 
    high1_low2_bool = [ !low_x1_bool[j] && low_x2_bool[j] for j = 1:nprime ] 
    high1_high2_bool = [ !low_x1_bool[j] && !low_x2_bool[j] for j = 1:nprime ] 

    low1_low2_fit_bool = [ low_x1fit_bool[j] && low_x2fit_bool[j] for j = 1:model.n ] 
    low1_high2_fit_bool = [ low_x1fit_bool[j] && !low_x2fit_bool[j] for j = 1:model.n ] 
    high1_low2_fit_bool = [ !low_x1fit_bool[j] && low_x2fit_bool[j] for j = 1:model.n ] 
    high1_high2_fit_bool = [ !low_x1fit_bool[j] && !low_x2fit_bool[j] for j = 1:model.n ] 


    n_low1 = sum(low_x1_bool)
    n_high1 = sum(.!low_x1_bool)

    n_low2 = sum(low_x2_bool)
    n_high2 = sum(.!low_x2_bool)

    n_low1_low2 = sum( low1_low2_bool )
    n_low1_high2 = sum( low1_high2_bool )
    n_high1_low2 = sum( high1_low2_bool )
    n_high1_high2 = sum( high1_high2_bool )


    n_low1_fit = sum(low_x1fit_bool)
    n_high1_fit = sum(.!low_x1fit_bool)

    n_low2_fit = sum(low_x2fit_bool)
    n_high2_fit = sum(.!low_x2fit_bool)

    n_low1_low2_fit = sum( low1_low2_fit_bool )
    n_low1_high2_fit = sum( low1_high2_fit_bool )
    n_high1_low2_fit = sum( high1_low2_fit_bool )
    n_high1_high2_fit = sum( high1_high2_fit_bool )


    write(report_file, "\n\nX1 below $(x1_cutoff), n_fit: $(n_low1_fit), n_valid: $(n_low1), mean KL: $(mean(KL_valid[low_x1_bool]))\n")
    write(report_file, "\nX1 above $(x1_cutoff), n_fit: $(n_high1_fit), n_valid: $(n_high1), mean KL: $(mean(KL_valid[.!low_x1_bool]))\n")

    write(report_file, "\n\nX2 below $(x2_cutoff), n_fit: $(n_low2_fit), n_valid: $(n_low2), mean KL: $(mean(KL_valid[low_x2_bool]))\n")
    write(report_file, "\nX2 above $(x2_cutoff), n_fit: $(n_high2_fit), n_valid: $(n_high2), mean KL: $(mean(KL_valid[.!low_x2_bool]))\n")

    write(report_file, "\n\nX1 below $(x1_cutoff) and X2 below $(x2_cutoff), n_fit: $(n_low1_low2_fit), n_valid: $(n_low1_low2), mean KL: $(mean(KL_valid[low1_low2_bool]))\n")
    write(report_file, "\nX1 below $(x1_cutoff) and X2 above $(x2_cutoff), n_fit: $(n_low1_high2_fit), n_valid: $(n_low1_high2), mean KL: $(mean(KL_valid[low1_high2_bool]))\n")
    write(report_file, "\nX1 above $(x1_cutoff) and X2 below $(x2_cutoff), n_fit: $(n_high1_low2_fit), n_valid: $(n_high1_low2), mean KL: $(mean(KL_valid[high1_low2_bool]))\n")
    write(report_file, "\nX1 above $(x1_cutoff) and X2 above $(x2_cutoff), n_fit: $(n_high1_high2_fit), n_valid: $(n_high1_high2), mean KL: $(mean(KL_valid[high1_high2_bool]))\n")
    
end


write(report_file, "\n\nFinal KL: $(KL_out)\n")
close(report_file)
