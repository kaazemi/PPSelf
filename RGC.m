clear; clc; close all
% fs = 50e-6;
binsize = 500;

load RGC_data;
rgc_obs = n13r10;

for t = 1:length(rgc_obs)/binsize;
data(t) = sum(rgc_obs(1+(t-1)*binsize:t*binsize))>0;
end
% n = floor(length(data)/3);
n = 60;

sys.pi_min = 0;
sys.pi_max = 0.49;
sys.p = 50;
sys.lambda = .25;
sys.s_star = 3;

sys.y = data(1:n+sys.p)';
sys.test_data = data(n+sys.p+1:end)';


sys_est =  Estimate_pp(sys);
plot_hawkes
