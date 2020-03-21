julia postsim_to_R.jl $1 &
wait
Rscript --vanilla postProcess_loop.R $1 &
wait
