using BNP_WMReg_Joint

using BSON
import Clustering.hclust, Clustering.cutree
import Distances.pairwise, Distances.Euclidean
import BayesInference.embed, BayesInference.etr
using PDMats
using StatsBase
using Statistics
using Random
using Dates

using LinearAlgebra
BLAS.set_num_threads(1)

## Settings from shell script
const srnd = parse(Int64, ARGS[1])
const n = parse(Int64, ARGS[2])
const K = parse(Int64, ARGS[3])
const H = parse(Int64, ARGS[4])
const Σx_type = Symbol(ARGS[5])
const γ_type = Symbol(ARGS[6])
const γ_init = ARGS[7] # "all" or "none"
const SNR = parse(Float64, ARGS[8])
const whichdat = ARGS[9]

## MCMC Settings (used in article)
# n_burn = Int(300e3)
# n_keep = 5000
# n_adapt = 500
# n_thin = 100
# centerscale = false

## Alternate settings (for shorter runs)
n_burn = Int(200e3)
n_keep = 1000
n_adapt = 500
n_thin = 50
centerscale = false

const mesg = "$(whichdat)_n$(n)_K$(K)_SigX$(Σx_type)_VS$(γ_type)_gaminit$(γ_init)_snr$(Int(floor(SNR)))_$(srnd)"
Random.seed!(srnd)

dat = BSON.load("data/$(whichdat).bson")
y = deepcopy(dat[:y])
X = deepcopy(dat[:X])

y_fit = deepcopy(y[1:n])
X_fit = deepcopy(X[1:n, 1:K])

if centerscale
    X_fit = (X_fit .- mean(y_fit)) ./ std(y_fit) # time series
    y_fit = (y_fit .- mean(y_fit)) ./ std(y_fit)
end

## Customize the prior
range_y = maximum(y_fit) - minimum(y_fit)
center_y = mean(y_fit)
range_X = fill(deepcopy(range_y), K) # assuming time series application
center_X = fill(deepcopy(center_y), K)

prior = Prior_BNP_WMReg_Joint(K, H, range_y = range_y, center_y = center_y,
    range_X = range_X, center_X = center_X, snr=SNR, Σx_type=Σx_type, γ_type=γ_type)

if γ_type == :local
    prior.π_ξ = (0.1 .+ (0.5 .^ collect(1:K) .* 0.4/0.5 ))
end

init = init_state_BNP_WMReg_Joint(n, K, H, prior, random=1, Σx_type=Σx_type, γ_type=γ_type) # random: 0 for none, 1 for lower level params, 2 for all

if γ_init == "all"
    ## default is for all to initialize on
    all(init.γ) || throw("Reset default lag initialization to all on.")
elseif γ_init == "none"
    all(init.γ) || throw("Reset default lag initialization to all on.")
    init.γ = .!init.γ
else
    throw("Invalid initialization for lag selectors.")
end

if γ_type == :global
    init.π_γ = (0.1 .+ (0.5 .^ collect(1:K) .* 0.4/0.5 ))
end

## Clustering algorithm to initialize allocation
if γ_init == "all"
    all(init.γ) || throw("Reset default lag initialization to all on.")

    clust = hclust(pairwise(Euclidean(), transpose(hcat(y_fit, X_fit)), dims=2), linkage=:ward)
    init.S = cutree(clust, k=H)
    init.n_occup = length(unique(init.S))
    println(counts(init.S, 1:H))

elseif γ_init == "none"
    all(.!init.γ) || throw("Reset default lag initialization to all on.")

    clust = hclust(abs.(broadcast(-, y_fit, permutedims(y_fit))), linkage=:ward)
    init.S = cutree(clust, k=H)
    init.n_occup = length(unique(init.S))
    println(counts(init.S, 1:H))

else
    throw("Invalid initialization for lag selectors.")
end


## log file
logfilename = "postsim_progress/out_prog_$(mesg).txt"

## initialize model
model = Model_BNP_WMReg_Joint(y_fit, X_fit, H, prior, init, Σx_type=Σx_type, γ_type=γ_type)

monitor = Monitor_BNP_WMReg_Joint(true, # ηlω
                              false, # γ
                              false, # S
                              false) # G0

updatevars = Updatevars_BNP_WMReg_Joint(true,  # η
                                    true,  # lω
                                    false, # γ
                                    true,  # α
                                    true,  # S
                                    false)  # G0

@time sims, accptr = mcmc!(model, 1000,
    updatevars, monitor, logfilename, 1, 100)
println(model.state.iter)
println(accptr)

adapt!(model, n_iter_collectSS=n_adapt, n_iter_scale=n_adapt,
    adapt_thin=10, updatevars=updatevars,
    report_filename=logfilename, maxtries=10,
    accptr_bnds=[0.02, 0.15], tune_only=true)

if model.γ_type != :fixed
    updatevars.γ = true
end
monitor.γ = true # always on
updatevars.G0 = true

adapt!(model, n_iter_collectSS=n_adapt, n_iter_scale=n_adapt,
    adapt_thin=10, updatevars=updatevars,
    report_filename=logfilename, maxtries=10,
    accptr_bnds=[0.02, 0.15], tune_only=true)

## time to explore posterior
sims, accptr = mcmc!(model, Int(floor(n_burn/10)),
    updatevars, monitor, logfilename, 2, 30000)

report_file = open(logfilename, "a+")
if model.γ_type == :global
    write(logfilename, "\n\nCurrent variable selection rates: $([ mean( [ 1*sims[i][:γ][k] for i=1:length(sims) ] ) for k = 1:model.K ])\n\n")
    println([ mean( [ 1*sims[i][:γ][k] for i=1:length(sims) ] ) for k = 1:model.K ])
elseif model.γ_type == :local
    γ_occupied = counts(model.state.S, 1:model.H)'model.state.γ ./ float(model.n)
    write(logfilename, "Current variable selection of occupied clusters: $(γ_occupied)\n\n")
    println( mean( [ sims[i][:Scounts]'sims[i][:γ] ./ float(model.n) for i=1:length(sims) ] ) )
end
write(report_file, "Pre-burn-in complete. Now adapting with variable selection on.\n\n")
close(report_file)

for ii = 1:2
    ## adapt with variable selection on
    adapt!(model, n_iter_collectSS=n_adapt, n_iter_scale=n_adapt,
        adapt_thin=10, updatevars=updatevars,
        report_filename=logfilename, maxtries=10,
        accptr_bnds=[0.02, 0.20], tune_only=true)

    ## time to explore posterior
    global sims, accptr = mcmc!(model, Int(floor(n_burn/10)),
        updatevars, monitor, logfilename, 1, 10000)

    report_file = open(logfilename, "a+")
    if model.γ_type == :global
        write(logfilename, "\n\nCurrent variable selection rates: $([ mean( [ 1*sims[i][:γ][k] for i=1:length(sims) ] ) for k = 1:model.K ])\n\n")
        println([ mean( [ 1*sims[i][:γ][k] for i=1:length(sims) ] ) for k = 1:model.K ])
    elseif model.γ_type == :local
        γ_occupied = counts(model.state.S, 1:model.H)'model.state.γ ./ float(model.n)
        write(logfilename, "Current variable selection of occupied clusters: $(γ_occupied)\n\n")
        println( mean( [ sims[i][:Scounts]'sims[i][:γ] ./ float(model.n) for i=1:length(sims) ] ) )
    end
    write(report_file, "Adaptation $(ii) of 3 complete.\n\n")
    close(report_file)
end

## burn-in
timestartburn = Dates.now()
sims, accptr = mcmc!(model, Int(floor(n_burn/10)),
    updatevars, monitor, logfilename, 10, 50000)
etr(timestartburn, n_keep, n_thin, logfilename)

report_file = open(logfilename, "a+")
if model.γ_type == :global
    write(logfilename, "\n\nCurrent variable selection rates: $([ mean( [ 1*sims[i][:γ][k] for i=1:length(sims) ] ) for k = 1:model.K ])\n\n")
    println([ mean( [ 1*sims[i][:γ][k] for i=1:length(sims) ] ) for k = 1:model.K ])
elseif model.γ_type == :local
    γ_occupied = counts(model.state.S, 1:model.H)'model.state.γ ./ float(model.n)
    write(logfilename, "Current variable selection of occupied clusters: $(γ_occupied)\n\n")
    println( mean( [ sims[i][:Scounts]'sims[i][:γ] ./ float(model.n) for i=1:length(sims) ] ) )
end
write(report_file, "Burn-in complete. Beginning inference samples.\n\n")
close(report_file)

## inference runs
@time global sims, accptr = mcmc!(model, n_keep,
    updatevars, monitor, logfilename, n_thin, 50000)

## save samples
using BSON

if haskey(dat, :truth)

    bson("postsim/postsim_$(mesg).bson",
        model=deepcopy(model), sims=deepcopy(sims),
        init=deepcopy(init), mesg=deepcopy(mesg),
        truth=deepcopy(dat[:truth]) )

else

    bson("postsim/postsim_$(mesg).bson",
        model=deepcopy(model), sims=deepcopy(sims),
        init=deepcopy(init), mesg=deepcopy(mesg) )

end
