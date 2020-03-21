## ARGS
# 1: srnd; seed for rng
# 2: n; sample size (T - L)
# 3: K; number of lags (L)
# 4: H; DP druncation
# 5: Σx_type; full or diag
# 6: γ_type; fixed, global, or local (lag selection)
# 7: γ_init; all or none (initialization of lag selection indicators)
# 8: SNR; Signal-to-noise parameter for the prior
# 9: whichdat; name of the dataset

nice julia runMCMC.jl 1121191 70 2 40 full fixed all 5 Ricker_single &> logs/Ricker_single_n70_K2_full_fixed_all_5_1.txt &
nice julia runMCMC.jl 1121192 70 2 40 full fixed all 5 Ricker_single &> logs/Ricker_single_n70_K2_full_fixed_all_5_2.txt &
nice julia runMCMC.jl 1121193 70 2 40 full fixed all 5 Ricker_single &> logs/Ricker_single_n70_K2_full_fixed_all_5_3.txt &

# nice julia runMCMC.jl 1121191 300 2 40 full fixed all 5 Ricker_single &> logs/Ricker_single_n300_K2_full_fixed_all_5_1.txt &
# nice julia runMCMC.jl 1121192 300 2 40 full fixed all 5 Ricker_single &> logs/Ricker_single_n300_K2_full_fixed_all_5_2.txt &
# nice julia runMCMC.jl 1121193 300 2 40 full fixed all 5 Ricker_single &> logs/Ricker_single_n300_K2_full_fixed_all_5_3.txt &


nice julia runMCMC.jl 1121191 70 5 40 diag global all 5 Ricker_single &> logs/Ricker_single_n70_K5_diag_global_all_5_1.txt &
nice julia runMCMC.jl 1121192 70 5 40 diag global all 5 Ricker_single &> logs/Ricker_single_n70_K5_diag_global_all_5_2.txt &
nice julia runMCMC.jl 1121193 70 5 40 diag global none 5 Ricker_single &> logs/Ricker_single_n70_K5_diag_global_none_5_1.txt &
nice julia runMCMC.jl 1121194 70 5 40 diag global none 5 Ricker_single &> logs/Ricker_single_n70_K5_diag_global_none_5_2.txt &

# nice julia runMCMC.jl 1121191 300 5 40 diag global all 5 Ricker_single &> logs/Ricker_single_n300_K5_diag_global_all_5_1.txt &
# nice julia runMCMC.jl 1121192 300 5 40 diag global all 5 Ricker_single &> logs/Ricker_single_n300_K5_diag_global_all_5_2.txt &
# nice julia runMCMC.jl 1121193 300 5 40 diag global none 5 Ricker_single &> logs/Ricker_single_n300_K5_diag_global_none_5_1.txt &
# nice julia runMCMC.jl 1121194 300 5 40 diag global none 5 Ricker_single &> logs/Ricker_single_n300_K5_diag_global_none_5_2.txt &


nice julia runMCMC.jl 1121191 70 5 40 diag local all 5 Ricker_single &> logs/Ricker_single_n70_K5_diag_local_all_5_1.txt &
nice julia runMCMC.jl 1121192 70 5 40 diag local all 5 Ricker_single &> logs/Ricker_single_n70_K5_diag_local_all_5_2.txt &
nice julia runMCMC.jl 1121193 70 5 40 diag local none 5 Ricker_single &> logs/Ricker_single_n70_K5_diag_local_none_5_1.txt &
nice julia runMCMC.jl 1121194 70 5 40 diag local none 5 Ricker_single &> logs/Ricker_single_n70_K5_diag_local_none_5_2.txt &

# nice julia runMCMC.jl 1121191 300 5 40 diag local all 5 Ricker_single &> logs/Ricker_single_n300_K5_diag_local_all_5_1.txt &
# nice julia runMCMC.jl 1121192 300 5 40 diag local all 5 Ricker_single &> logs/Ricker_single_n300_K5_diag_local_all_5_2.txt &
# nice julia runMCMC.jl 1121193 300 5 40 diag local none 5 Ricker_single &> logs/Ricker_single_n300_K5_diag_local_none_5_1.txt &
# nice julia runMCMC.jl 1121194 300 5 40 diag local none 5 Ricker_single &> logs/Ricker_single_n300_K5_diag_local_none_5_2.txt &

## 22 jobs

wait
