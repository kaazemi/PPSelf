close all;
disp('Plotting the Estimates')
%%
if 1 
    m = min([sys_est.theta_ML;sys_est.theta_sp;sys_est.theta_OMP])-.01;
    M = max([sys_est.theta_ML;sys_est.theta_sp;sys_est.theta_OMP]) + 0.01;
    figure;
    subplot(4,1,1); stem(sys.theta_gt,'linewidth',2); title('ML');
    ylim([m M]);
    subplot(4,1,2); stem(sys_est.theta_ML,'linewidth',2); title('ML');
    ylim([m M]);
    subplot(4,1,3); stem(sys_est.theta_sp,'linewidth',2); title('ell_1');
    ylim([m M]);
    subplot(4,1,4); stem(sys_est.theta_OMP,'linewidth',2); title('OMP');
    ylim([m M]);
    set(gcf,'units','normalized','outerposition',[0 0 .55 0.95],'defaulttextinterpreter','latex')

    KS_ACF(sys.test_data, sys_est.theta_ML,sys_est.mu_ML); title('ML estimate')
    KS_ACF(sys.test_data, sys_est.theta_sp,sys_est.mu_sp); title('ell_regularized estimate')
    KS_ACF(sys.test_data, sys_est.theta_OMP,sys_est.mu_OMP); title('OMP estimate')
end
