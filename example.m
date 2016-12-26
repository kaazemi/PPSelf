clear; clc; close all;
p = 50; s = 3; n = 1000;
sys = gen_synthetic_data( p, s, n );
sys_est =  Estimate_pp(sys);
plot_hawkes_synthetic
