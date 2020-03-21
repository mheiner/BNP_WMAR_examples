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

nice julia runMCMC.jl 1121191 25 5 25 diag global all 5 pinkSalmon_first30_log &> logs/pinkSalmon_first30_log_n25_K5_diag_global_all_5_1.txt &
nice julia runMCMC.jl 1121192 25 5 25 diag global all 5 pinkSalmon_first30_log &> logs/pinkSalmon_first30_log_n25_K5_diag_global_all_5_2.txt &
nice julia runMCMC.jl 1121193 25 5 25 diag global none 5 pinkSalmon_first30_log &> logs/pinkSalmon_first30_log_n25_K5_diag_global_none_5_1.txt &
nice julia runMCMC.jl 1121194 25 5 25 diag global none 5 pinkSalmon_first30_log &> logs/pinkSalmon_first30_log_n25_K5_diag_global_none_5_2.txt &

nice julia runMCMC.jl 1121191 25 5 25 diag local all 5 pinkSalmon_first30_log &> logs/pinkSalmon_first30_log_n25_K5_diag_local_all_5_1.txt &
nice julia runMCMC.jl 1121192 25 5 25 diag local all 5 pinkSalmon_first30_log &> logs/pinkSalmon_first30_log_n25_K5_diag_local_all_5_2.txt &
nice julia runMCMC.jl 1121193 25 5 25 diag local none 5 pinkSalmon_first30_log &> logs/pinkSalmon_first30_log_n25_K5_diag_local_none_5_1.txt &
nice julia runMCMC.jl 1121194 25 5 25 diag local none 5 pinkSalmon_first30_log &> logs/pinkSalmon_first30_log_n25_K5_diag_local_none_5_2.txt &

## 8 jobs

wait
