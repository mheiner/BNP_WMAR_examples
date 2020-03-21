using BSON
using RCall

R"load('Ricker_sigle_lognormal.rda')"
R"ls()"

@rget y X y_valid X_valid Y_valid tprime

bson("Ricker_sigle_lognormal.bson", y = deepcopy(y),
    X = deepcopy(X), y_valid = deepcopy(y_valid),
    X_valid = deepcopy(X_valid), Y_valid = deepcopy(Y_valid),
    tprime = deepcopy(tprime))

R"rm(list=ls())"

R"load('Ricker_twolag_lognormal.rda')"
R"ls()"

@rget y X y_valid X_valid Y_valid tprime

bson("Ricker_twolag_lognormal.bson", y = deepcopy(y),
    X = deepcopy(X), y_valid = deepcopy(y_valid),
    X_valid = deepcopy(X_valid), Y_valid = deepcopy(Y_valid),
    tprime = deepcopy(tprime))

R"rm(list=ls())"

R"load('Ricker_single.rda')"
R"ls()"

@rget y X y_valid X_valid Y_valid tprime

bson("Ricker_single.bson", y = deepcopy(y),
    X = deepcopy(X), y_valid = deepcopy(y_valid),
    X_valid = deepcopy(X_valid), Y_valid = deepcopy(Y_valid),
    tprime = deepcopy(tprime))



### plotting

# using Plotly
# nplot = 500

# trace0 = scatter3d(Dict(
#     :x => X[1:nplot,1],
#     :y => X[1:nplot,2],
#     :z => y[1:nplot],
#     :marker => Dict(
#         :color => "rgba(0,0,0,0.70)",
#         :size => 2
#     ),
#     :mode => "markers"
# ))
# data = [trace0]
# maximum(X[1:nplot,1:2])

# layout = Layout(Dict(
#   :hovermode => "closest",
#   :margin => Dict(
#     :r => 10,
#     :t => 25,
#     :b => 40,
#     :l => 60
#   ),
#   :scene => Dict(
#     :xaxis => Dict(:title => "y[t-1]",
#                    :range => [0.0, 12.0]),
#     :yaxis => Dict(:title => "y[t-2]",
#                     :range => [0.0, 12.0]),
#     :zaxis => Dict(:title => "y[t]",
#                     :range => [0.0, 12.0]),
#   :showlegend => false
# )))

# p = plot(data, layout)

# nplot = 250
# Rplot = R"plot"
#
# R"pdf('appdata2_plots/appdata2_lag2only_hsced.pdf', width=4.5, height=3.5);
#   par(mar=c(4,4,1,1)+0.1)"
# Rplot(X[1:nplot,2], y[1:nplot], ylab="y[t]", xlab="y[t-2]", main="")
# R"dev.off()"
#
# R"pdf('appdata2_plots/appdata2_lag2only.pdf', width=4.5, height=3.5);
#   par(mar=c(4,4,1,1)+0.1)"
# Rplot(X[1:nplot,2], y[1:nplot], ylab="y[t]", xlab="y[t-2]", main="")
# R"dev.off()"
