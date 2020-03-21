
fileid = parse(Int64, ARGS[1])
simregexp = Regex("postsim.*_.*$(fileid).*.bson")


using BNP_WMReg_Joint
using RCall

using BSON
using BSON: @load

using Statistics
using PDMats
using StatsBase
using LinearAlgebra


files = filter(x -> occursin(simregexp, x), readdir("./postsim"))
nfiles = length(files)

i = 1
for i in 1:nfiles
    file_now = files[i]
    println("starting file: $(file_now)")

    BSON.@load "./postsim/$(file_now)" sims mesg model init
    nsim = length(sims)

    ### Send to R
    R"rm(list=ls())"
    iter = model.state.iter
    accptr =  model.state.accpt ./ model.state.iter
    K = model.K
    H = model.H
    n = model.n
    priorinfo = "$(model.prior)"
    initinfo = "$(init)"

    Sigx_type = "$(model.Σx_type)"
    gam_type = "$(model.γ_type)"

    sims_llik = [ sims[i][:llik] for i=1:nsim ]
    sims_noccup = [ sims[i][:n_occup] for i=1:nsim ]
    sims_Scounts = permutedims(hcat([ sims[i][:Scounts] for i=1:length(sims) ]...))

    if model.γ_type in (:global, :fixed)

        sims_gam = permutedims(hcat([ sims[i][:γ] for i=1:length(sims) ]...))
        sims_gam_occ = "NA"
        sims_gam_w = "NA"
        sims_gam_glob = "NA"

        if haskey(sims[1], :lfc_on)
            sims_fc_on = permutedims(hcat([ exp.(sims[i][:lfc_on]) for i=1:length(sims) ]...))
            sims_pigam = deepcopy(model.state.π_γ)
        else
            sims_fc_on = "NA"
            sims_pigam = deepcopy(model.state.π_γ)
        end

        b = 0.05
        sims_gam_alt = [ ( sims[ii][:γ][k] && abs(sims[ii][:β_y][h,k]) > b )  for ii = 1:length(sims), h = 1:model.H, k = 1:model.K ]
        sims_gam_alt_occup = [ sims[ii][:Scounts]'sims_gam_alt[ii,:,k] / model.n  for ii = 1:length(sims), k = 1:model.K ]

    elseif model.γ_type == :local
        
        sims_gam = "NA"
        sims_gam_occ = permutedims(hcat([ vec(sims[i][:γ_occupied]) for i=1:length(sims) ]...))
        sims_gam_w = [ sum( exp.(sims[i][:lω][sims[i][:γ][:,k]]) ) for i=1:length(sims), k = 1:model.K ]
        sims_gam_glob = [ any( sims[i][:γ][ 1:maximum( findall( sims[i][:Scounts] .> 0 ) ) , k] ) for i=1:length(sims), k = 1:model.K ]

        sims_fc_on = "NA"
        sims_pigam = permutedims(hcat([ sims[i][:π_γ] for i=1:length(sims) ]...))

        b = 0.05
        sims_gam_alt = [ ( sims[ii][:γ][h,k] && abs(sims[ii][:β_y][h,k]) > b )  for ii = 1:length(sims), h = 1:model.H, k = 1:model.K ]
        sims_gam_alt_occup = [ sims[ii][:Scounts]'sims_gam_alt[ii,:,k] / model.n  for ii = 1:length(sims), k = 1:model.K ]

    end
    
    sims_w = permutedims(hcat([ exp.(sims[i][:lω]) for i=1:length(sims) ]...))
    sims_alpha = [ sims[i][:α] for i=1:nsim ]

    sims_mu = permutedims(hcat([ sims[i][:μ_y] for i=1:length(sims) ]...))
    sims_beta = permutedims(cat([ sims[i][:β_y] for i=1:length(sims) ]..., dims=3), (3,2,1))
    sims_sig2 = permutedims(hcat([ sims[i][:δ_y] for i=1:length(sims) ]...))

    if model.γ_type in (:global, :fixed)
        sims_int = [ sims[i][:μ_y][h] + sum( sims[i][:γ] .* sims[i][:β_y][h,:] .* sims[i][:μ_x][h,:]) for i=1:length(sims), h=1:H ]
    elseif model.γ_type == :local
        sims_int = [ sims[i][:μ_y][h] + sum( sims[i][:γ][h,:] .* sims[i][:β_y][h,:] .* sims[i][:μ_x][h,:]) for i=1:length(sims), h=1:H ]
    end


    @rput mesg nsim iter accptr K H n priorinfo initinfo Sigx_type gam_type sims_llik sims_noccup sims_Scounts sims_gam sims_gam_occ sims_gam_w sims_gam_glob sims_fc_on sims_pigam sims_gam_alt_occup sims_w sims_alpha sims_mu sims_beta sims_sig2 sims_int

    R"ls()"
    R"save.image(file=paste0('postsim/mcmc_', $(mesg), '.rda'))"

    cp("./postsim/$(file_now)", "./postsimB/$(file_now)")
    rm("./postsim/$(file_now)")

    println("completed file number: $(i) of $(nfiles) \n\n")
end
