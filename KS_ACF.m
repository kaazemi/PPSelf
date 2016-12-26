function KS_ACF(test_data,theta,mu)

delta = 1;
%% Calculate ISI distribution Empirically
pd = gen_data(theta,mu,delta);

n_test = min(length(test_data),10000);
p = length(theta);

test_data = test_data(1:n_test);
rate = zeros(size(test_data));
for t = p:length(test_data)
rate(t+1) = delta*(mu+ test_data(t:-1:t-p+1)'*theta);
end

spike_times = find(test_data ==1);
spike_times = spike_times(spike_times>p);
num_spikes = length(spike_times);

for sp = 1:num_spikes-1
int_rate(sp) = sum(rate(spike_times(sp)+1:spike_times(sp+1)));
end
unif_rate = cdf(pd,int_rate);

%%
unif_sorted = sort(unif_rate);
unifs = ((1:num_spikes-1)-1/2)/(num_spikes-1);

% 1.36 for 95% (harder to pass)
% 1.63 for 99% (easier to pass)
upper_bound = unifs + 1.36/sqrt(num_spikes-1);
lower_bound = unifs - 1.36/sqrt(num_spikes-1);


figure
subplot(1,2,1);
plot([0 unif_sorted 1],[0 unifs 1],'LineWidth',2); hold on
plot([0 unifs 1],[0  unifs 1],'k--','LineWidth',2);
% plot([0 exp_sorted 1],[0 unifs 1]);
plot(upper_bound,unifs,'k--','LineWidth',2)
plot(lower_bound,unifs,'k--','LineWidth',2)
xlim([0 1])
ylim([0 1])
title('Kolmogorov Smirnov Test')
xlabel('Uniform Quantiles')
ylabel('Empirical CDF')
axis square

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ACF %%%%%%%%%%%%
equal1 = find(unif_rate == 1);
unif_rate(equal1) = unif_rate(max(equal1-1,1));
equal0 = find(unif_rate == 0);
unif_rate(equal0) = unif_rate(max(equal0-1,1));

v = norminv(unif_rate,0,1);
% 1.96 for 95% (harder to pass)
% 2.575 for 99% (easier to pass);
subplot(1,2,2); autocorr_test(v,[],[],2.575); axis square
set(gcf,'units','normalized','outerposition',[0 0 .5 0.6],'defaulttextinterpreter','latex')

end

%% Calculate ISI distribution Empirically
function pd = gen_data(theta,mu,delta)
    p = length(theta);
    y = zeros(1000*p,1);
    rate = zeros(1000*p,1);

    for t = p:1000*p
        rate(t+1) = mu+ y(t:-1:t-p+1)'*theta;
        y(t+1) = rand <= rate(t);
    end

    spike_times = find(y ==1);
    spike_times = spike_times(spike_times>p);
    num_spikes = length(spike_times);
    int_rate = zeros(num_spikes-1,1);

    for sp = 1:num_spikes-1
        int_rate(sp) = sum(rate(spike_times(sp)+1:spike_times(sp+1)));
    end
pd = learn_cdf(int_rate);

end

function pd = learn_cdf(data)
[f, a] = ecdf(data);
a(1) = a(1)-.000001;
pd = makedist('PieceWiselinear','x',a,'Fx',f);
end
