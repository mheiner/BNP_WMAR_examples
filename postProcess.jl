using BNP_WMReg_Joint

using Statistics
using Plotly

using PDMats
using StatsBase

using BSON
using BSON: @load

println(pwd())

## load
using LinearAlgebra # along with all others above

whichdat = "pinkSalmon_first30_log_n25_K5_SigXdiag_VSlocal_gaminitall_snr5_1121191"

whichdat = "Ricker_single_..."
whichdat = "Ricker_single_lognormal_..."
whichdat = "Ricker_twolag_lognormal_..."

BSON.@load "postsimB/postsim_$(whichdat).bson" sims model init mesg




### posterior checks
println(model.state.llik)

println(model.state.iter)
println(model.state.accpt / float(model.state.iter))
println(counts(model.state.S, 1:model.H))

println(exp.(model.state.lω))
println(model.state.α)

iter_use = 1:length(sims) # all simulations
iter_use = 1:1000 # first thousand
iter_use = range(1, stop=length(sims), step=5) # by 5
length(collect(iter_use))

plot( [ sims[i][:llik] for i in iter_use ] )

plot( [ sims[i][:n_occup] for i in iter_use ] )

println( round.( mean( hcat([ sims[i][:Scounts] for i in iter_use ]...), dims=2 ), digits=1 ) )
println( round.( mean( hcat([ sims[i][:Scounts] .> 0 for i in iter_use ]...), dims=2 ), digits=1 ) )
println( sum( mean( hcat([ sims[i][:Scounts] .> 0 for i in iter_use ]...), dims=2 ) .> 0.0 ) )

h = 1
plot( [ sims[i][:Scounts][h] for i in iter_use ] )



### lag selection

## global
k = 1
plot( [ 1*sims[i][:γ][k] for i in iter_use ] )
println( mean( [ 1*sims[i][:γ][k] for i in iter_use, k = 1:model.K ], dims=1 ) )


## local
γoccup = permutedims( hcat( [ vec(sims[ii][:γ_occupied]) for ii in iter_use ]... ) )
pm_γoccup = mean( γoccup, dims=1 )[1,:]
println(pm_γoccup)

b = 0.05
γ_alt = [ ( sims[iter_use[ii]][:γ][h,k] && abs(sims[iter_use[ii]][:β_y][h,k]) > b )  for ii = 1:length(iter_use), h = 1:model.H, k = 1:model.K ]
γ_alt_occup = [ sims[iter_use[ii]][:Scounts]'γ_alt[ii,:,k] / model.n  for ii = 1:length(iter_use), k = 1:model.K ]
pm_γ_alt_occup = mean( γ_alt_occup, dims=1 )[1,:]
println(pm_γ_alt_occup)

sims_gam_glob = [ any( sims[i][:γ][ 1:maximum( findall( sims[i][:Scounts] .> 0 ) ) , k] ) for i in iter_use, k = 1:model.K ]
pm_gam_glob = mean(sims_gam_glob, dims=1)[1,:]

k = 1
plot( sims_gam_glob[:,k] )

mean( [ sims[i][:π_γ] for i in iter_use ] )

k = 3
plot( [ sims[i][:π_γ][k] for i in iter_use ] )

mean( [ 1.0*sims[i][:ξ] for i in iter_use ] )

k = 2
plot( [ sims[i][:ξ][k] for i in iter_use ] )

sims_gam_w = [ sum( exp.(sims[i][:lω][sims[i][:γ][:,k]]) ) for i=1:length(sims), k = 1:model.K ]

k = 3
plot( sims_gam_w[iter_use, k] )
plot( γoccup[iter_use, k] )


h = 1
k = 2
println( round.( mean( hcat([ sims[i][:Scounts] for i in iter_use ]...), dims=2 ), digits=1 ) )

plot( [ 1*sims[i][:γ][h,k] for i in iter_use ] )
for h = 1:model.H
    println(mean( hcat([ 1*sims[i][:γ][h,:] for i in iter_use ]...), dims=2 ))
end

maxScount_indx = [ findmax(sims[i][:Scounts])[2] for i in 1:length(sims) ]
plot( [ 1*sims[i][:γ][maxScount_indx[i],k] for i in iter_use ] )





### component kernel parameters

model.state.δ_x
h = 1
k = 1
plot( [sims[i][:δ_x][h,k] for i in iter_use] )
Statistics.median([sims[i][:δ_x][h,k] for i in iter_use])


model.state.μ_x
h = 2
k = 1
plot([sims[i][:μ_x][h,k] for i in iter_use])
mean([sims[i][:μ_x][h,k] for i in iter_use])
std(sims.μ_x, dims=1)[1,:,:]


model.state.β_x
h = 1
kk = 2
k = 1
plot( [sims[i][:β_x][kk][h,k] for i in iter_use] )


model.state.μ_y
h = 2
plot([sims[i][:μ_y][h] for i in iter_use])


model.state.β_y
h = 1
k = 1
plot( [sims[i][:β_y][h,k] for i in iter_use] )
Statistics.median([sims[i][:β_y][h,k] for i in iter_use])
Statistics.mean([sims[i][:β_y][h,k] for i in iter_use])
println(Statistics.std([sims[i][:β_y][h,k] for i in iter_use]))

Statistics.mean([sims[i][:β_y][h,k] for i in iter_use])
println(Statistics.std([sims[i][:β_y][h,k] for i in iter_use]))

## prior covariance for beta_y
inv(model.state.Λ0star_ηy.mat) * model.state.δ_y[h]
model.prior.s0_δy_s0 * model.prior.Λ0star_ηy_S0.mat

model.state.δ_y
h = 2
plot( [sims[i][:δ_y][h] for i in iter_use] )
median([sims[i][:δ_y][h] for i in iter_use])


# intercepts
h = 3

# local lag selection
plot( [sims[i][:μ_y][h] + sum( sims[i][:γ][h,:] .* sims[i][:β_y][h,:] .* sims[i][:μ_x][h,:]) for i in iter_use ] )

# global lag selection
plot( [sims[i][:μ_y][h] + sum( sims[i][:γ] .* sims[i][:β_y][h,:] .* sims[i][:μ_x][h,:]) for i in iter_use ] )
Statistics.mean([sims[i][:μ_y][h] + sum( sims[i][:γ] .* sims[i][:β_y][h,:] .* sims[i][:μ_x][h,:]) for i in iter_use ])
println(Statistics.std([sims[i][:μ_y][h] + sum( sims[i][:γ] .* sims[i][:β_y][h,:] .* sims[i][:μ_x][h,:]) for i in iter_use ]))



h = 4
plot(([sims[i][:lω][h] for i in iter_use]))
plot(exp.([sims[i][:lω][h] for i in iter_use]))

vec(mean( exp.( permutedims(hcat([sims[i][:lω] for i in iter_use]...)) ) , dims=1 ))
vec(std( exp.( permutedims(hcat([sims[i][:lω] for i in iter_use]...)) ) , dims=1 ))
counts(model.state.S, 1:model.H)

plot(exp.([sims[i][:lω][h] for i in iter_use]))
plot(exp.([sims[i][:lω][model.H] for i in iter_use]))


plot([sims[i][:α] for i in iter_use])








### Functionals

include("Rfuns.jl")

npred = 25

minx = minimum(model.X)
maxx = maximum(model.X)
ranx = maxx - minx

X_pred_vals = collect(range(minx - 0.01*ranx, length=npred,
                        stop=(maxx + 0.01*ranx)))

## γ_use is a boolean vector specifying which lags we want to examine
γ_use = deepcopy(model.state.γ) # global only
γ_use = [false, true, false, false, false]
γ_use = pm_γoccup .> 0.11 # local
γ_use = pm_gam_glob .> 0.4 # local 
println(1*γ_use)

iter_use = sort(findall([ sims[ii][:γ] == γ_use for ii in 1:length(sims)])) # global only
iter_use = range(1, stop=length(sims), step=10)
iter_use = 1:(length(sims))

length(iter_use)

γindx_use = findall(γ_use)
Kuse = sum(γ_use)
X_pred00 = rcopy(xpndgrid([X_pred_vals for k = 1:Kuse]...))

X_pred0 = hcat([ fill(mean(model.X[:,k]), size(X_pred00,1)) for k = 1:model.K ]...)
X_pred0 = hcat([ rand(size(X_pred00,1)) .* (maximum(model.X[:,k]) - minimum(model.X[:,k])) .+ minimum(model.X[:,k])  for k = 1:model.K ]...) # draw uniformly over the range

X_pred0[:,findall(γ_use)] = deepcopy(X_pred00)
X_pred = deepcopy(X_pred0)

# X_pred = deepcopy(model.X) # for residual analysis

## use posterior samples to recreate mixture weights
ldw = ldensweight_mat(X_pred, sims[iter_use], model.state.γδc)
dw = exp.(ldw)
sum(dw[1,1,:])

i = 1
ℓ = 2
h = 2
plot([scatter(x=X_pred[:,ℓ], y=dw[i,:,h], mode="scatter")])

## examine weight function
mean_dw = mean(dw, dims=1)[1,:,:]

println( round.( mean( hcat([ sims[i][:Scounts] for i in iter_use ]...), dims=2 ), digits=1 ) )

h = 2
plot( [ sims[i][:Scounts][h] for i in iter_use ] )

mean_dw_mat = Matrix(reshape(mean_dw[:,h], fill(npred,Kuse)...)) # weight "surface"
trace1 = surface(Dict(
  :x => X_pred_vals,
  :y => X_pred_vals,
  :z => mean_dw_mat,
  # :z => q025_Q_mat_plot,
  :colorscale => "Viridis",
  :opacity => 0.7,
  :showscale => false,
  :type => "surface"
))
data = [trace1]
plot(data)

ℓ = 1
h = 3
plot([scatter(x=X_pred[:,ℓ], y=mean_dw[:,h], mode="scatter")]) # 1D

plot([scatter(x=X_pred[:,ℓ], y=mean_dw[1,:,h], mode="scatter")])
plot(scatter(x=X[:,1], y=y, mode="markers"))


## use weight function to create transition mean function (surface)

Ey = getEy(X_pred, dw, sims[iter_use])

# ii = 500
# plot([scatter(x=X_pred[:,1], y=Ey[ii,:], mode="scatter")])

mean_Ey = mean(Ey, dims=1)[1,:]

q025_Ey = [quantile(Ey[:,j], 0.025) for j = 1:size(Ey,2)]
q50_Ey = [quantile(Ey[:,j], 0.5) for j = 1:size(Ey,2)]
q975_Ey = [quantile(Ey[:,j], 0.975) for j = 1:size(Ey,2)]

resid = model.y - mean_Ey # ONLY use this if X_pred is the X matrix from the model
plotR(mean_Ey, resid, ylab="residual", xlab="posterior mean predicted value")


## transition quantile function (surface)

quant = 0.2
Q = getQuant(quant, X_pred, dw, sims[iter_use])

mean_Q = mean(Q, dims=1)[1,:]

q025_Q = [quantile(Q[:,j], 0.025) for j = 1:size(Q,2)]
q975_Q = [quantile(Q[:,j], 0.975) for j = 1:size(Q,2)]


### plot transition surfaces

### two dimensions

mean_Ey_mat = Matrix(reshape(mean_Ey, fill(npred,Kuse)...))
maximum(mean_Ey_mat)
minimum(mean_Ey_mat)
mean_Ey_mat_plot = [ abs(mean_Ey_mat[i,j]) < maximum(abs.([minx, maxx])) + 0.6*ranx ? mean_Ey_mat[i,j] : NaN for i = 1:npred, j=1:npred ]
mean_Ey_mat_plot = [ abs(mean_Ey_mat[i,j]) < maxx + 0.3*ranx ? mean_Ey_mat[i,j] : NaN for i = 1:npred, j=1:npred ]
# plot([scatter(x=X_pred[:,1], y=mean_Ey, mode="scatter")])

q025_Ey_mat = Matrix(reshape(q025_Ey, fill(npred,Kuse)...))
q025_Ey_mat_plot = [ abs(q025_Ey_mat[i,j]) < maxx + 0.1*ranx ? q025_Ey_mat[i,j] : NaN for i = 1:npred, j=1:npred ]

q975_Ey_mat = Matrix(reshape(q975_Ey, fill(npred,Kuse)...))
q975_Ey_mat_plot = [ abs(q975_Ey_mat[i,j]) < maxx + 0.1*ranx ? q975_Ey_mat[i,j] : NaN for i = 1:npred, j=1:npred ]


## quantile
mean_Q_mat = Matrix(reshape(mean_Q, fill(npred,Kuse)...))
maximum(mean_Q_mat)
mean_Q_mat_plot = [ abs(mean_Q_mat[i,j]) < maxx + 0.3*ranx ? mean_Q_mat[i,j] : NaN for i = 1:npred, j=1:npred ]

q025_Q_mat = Matrix(reshape(q025_Q, fill(npred,Kuse)...))
q025_Q_mat_plot = [ abs(q025_Q_mat[i,j]) < maxx + 0.1*ranx ? q025_Q_mat[i,j] : NaN for i = 1:npred, j=1:npred ]
q975_Q_mat = Matrix(reshape(q975_Q, fill(npred,Kuse)...))
q975_Q_mat_plot = [ abs(q975_Q_mat[i,j]) < maxx + 0.1*ranx ? q975_Q_mat[i,j] : NaN for i = 1:npred, j=1:npred ]



trace1 = surface(Dict(
  :x => X_pred_vals,
  :y => X_pred_vals,
  :z => q025_Ey_mat_plot,
  # :z => q025_Q_mat_plot, # uncomment to view quantile surfaces
  :colorscale => "Viridis",
  :opacity => 0.7,
  :showscale => false,
  :type => "surface"
))

trace2 = surface(Dict(
  :x => X_pred_vals,
  :y => X_pred_vals,
  :z => mean_Ey_mat_plot,
  # :z => mean_Q_mat_plot, # uncomment to view quantile surfaces
  :colorscale => "Viridis",
  :opacity => 0.9,
  :showscale => false,
  :type => "surface"
))

trace3 = surface(Dict(
  :x => X_pred_vals,
  :y => X_pred_vals,
  :z => q975_Ey_mat_plot,
  # :z => q975_Q_mat_plot, # uncomment to view quantile surfaces
  :colorscale => "Viridis",
  :opacity => 0.7,
  :showscale => false,
  :type => "surface"
))

trace0 = scatter3d(Dict(
    :x => model.X[:,γindx_use[1]],
    :y => model.X[:,γindx_use[2]],
    :z => model.y,
    :marker => Dict(
        :color => "rgba(0,0,0,0.70)",
        :size => 2
    ),
    :mode => "markers"
))

data = [trace2, trace0] # mean surface and data only
data = [trace1, trace2, trace3, trace0] # include intervals

layout = Layout(Dict(
  :hovermode => "closest",
  :margin => Dict(
    :r => 10,
    :t => 25,
    :b => 40,
    :l => 60
  ),
  :scene => Dict( :aspectmode => "manual", # doesn't seem to be working
    :xaxis => Dict(:title => "y[t-$(γindx_use[1])]",
                   # :range => [minx - 0.1*ranx, maxx + 0.1*ranx]),
                   :autorange => true),
    :yaxis => Dict(:title => "y[t-$(γindx_use[2])]",
                   # :range => [minx - 0.1*ranx, maxx + 0.1*ranx]),
                   :autorange => true),
    :zaxis => Dict(:title => "y[t]",
                   # :range => [minx - 0.1*ranx, maxx + 0.1*ranx])
                   :autorange => true),
  ),
  :showlegend => false
))

p = plot(data, layout)





### one dimension
size(Ey)
plotR(X_pred_vals, mean_Ey, ylab="", xlab="",
    ylim=[minimum(q025_Ey), maximum(q975_Ey)], lwd=2, typ="l")
plotR(X_pred_vals, mean_Ey, ylab="", xlab="",
    ylim=[minimum(y), maximum(y)], lwd=2, typ="l")

linesR(X_pred_vals, q025_Ey, lty=2)
linesR(X_pred_vals, q975_Ey, lty=2)
pointsR(model.X[:,γindx_use[1]], model.y)

R"title(main='Pink Salmon')"
R"title(ylab='log y[t]')"
R"title(xlab='log y[t-2] ')"
# linesR(X_pred_vals, mean_Q, lwd=2, typ="l", col="blue")
# linesR(X_pred_vals, p05_Q, lty=2, col="blue")
# linesR(X_pred_vals, p95_Q, lty=2, col="blue")

## individual iteration mean functions
for i in 1:5
    linesR(X_pred_vals, Ey[i,:], col=i, lwd=1)
end


R"par(mar=c(4,4,1,1)+0.1)"
plotR(X_pred_vals, mean_Ey, ylab="", xlab="",
    # ylim=[minimum(q025_Ey), maximum(q975_Ey)],
    ylim = [-1.0, 13.0], # for pink salmon
    lwd=2, typ="l")
polygonR(vcat(X_pred_vals, reverse(X_pred_vals)), vcat(q025_Ey, reverse(q975_Ey)), col="gray90", border="NA")
pointsR(model.X[:,γindx_use], model.y)
linesR(X_pred_vals, mean_Ey, lwd=2)

R"title(main='Pink Salmon');
    title(ylab='log y[t]');
    title(xlab='log y[t-2]')"
R"abline(0,1, lty=3)"

R"dev.off()"


plotR(model.y, typ="o", ylab="") # time series
plotR(model.X[:,2], model.y, ylab="", xlab="") # lag scatter plot






### Density estimation
## This section requires several variables assigned at the beginning of the 'Functionals' section.

using Distributions

n_y = 25
y_grid = collect(range(minx - 0.3*ranx, length=n_y,
                        stop=(maxx + 0.3*ranx)))
# y_grid = collect(range(minimum(model.y)-2.0, length=n_y, stop=maximum(model.y)+2.0))


n_star_Py = 4 # number of x locations at which to estimate the transition density

## single lag
n_star_Py = 3
X_star_Py = deepcopy(X_pred[1:n_star_Py,:])
X_star_Py[1, [γindx_use[1]] ] = [ 1.0 ]
X_star_Py[2, [γindx_use[1]] ] = [ 3.0 ]
X_star_Py[3, [γindx_use[1]] ] = [ 5.0 ]

## two lag
X_star_Py = deepcopy(X_pred[1:n_star_Py,:])
X_star_Py[1, [γindx_use[1], γindx_use[2]] ] = [ 1.0, 2.0 ]
X_star_Py[2, [γindx_use[1], γindx_use[2]] ] = [ 2.0, 1.0 ]
X_star_Py[3, [γindx_use[1], γindx_use[2]] ] = [ 2.5, 2.5 ]
X_star_Py[4, [γindx_use[1], γindx_use[2]] ] = [ 3.0, 2.5 ]

# iter_use = 1:length(sims)
iter_use = range(1, stop=length(sims), step=5)

X_star_Py

## construct weight functions and transition densities from posterior samples
ldwPy = ldensweight_mat(X_star_Py, sims[iter_use])
lPy, Ey_0 = getlogdens_EY(X_star_Py, y_grid, ldwPy, sims[iter_use])

Py = exp.(lPy)
size(Py)

mean_Py = permutedims(mean(Py, dims=1)[1,:,:])
q025_Py = [quantile(Py[:,j,i], 0.025) for i = 1:n_y, j = 1:n_star_Py]
# p50_Py = [quantile(Py[:,j,i], 0.5) for i = 1:n_y, j = 1:n_star_Py]
q975_Py = [quantile(Py[:,j,i], 0.975) for i = 1:n_y, j = 1:n_star_Py]

@rput mean_Py
R"library(lattice); levelplot(mean_Py)" # quick, discretized estimate of transtition density

xind = 1 # select row of n_star_Py for which to evaluate density
println(X_star_Py[xind,:])

## Transition densities for individual MCMC iterations
ii = 1 # posterior sample
j = 1 # row of n_star_Py
plotR(y_grid, Py[ii,j,:], typ="l", ylab="", xlab="")
for ii in range(1, stop=1000, step=3)
    linesR(y_grid, Py[ii,j,:], typ="l", ylab="", xlab="", col=ii)
end
R"dev.off()"

R"par(mfrow=c(2,2))"
for xind in 1:n_star_Py
    plotR(y_grid, q025_Py[:,xind],
        main="y[t-1] = $(round(X_star_Py[xind,γindx_use[1]], digits=2))",
        ylab="density", xlab="y[t] (minutes)", typ="l", ylim=[0.0, maximum(q975_Py[:,xind])], lty=2)
    linesR(y_grid, q975_Py[:,xind], lty=2)
    linesR(y_grid, mean_Py[:,xind], lty=1, lwd=2)
end
R"dev.off()"

