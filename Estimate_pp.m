function  sys_est = Estimate_pp(sys)
y = sys.y;
p = sys.p; n = length(y)-p;
lambda = sys.lambda;
A = hankel(y(1:n),y(n:end-1));  % Matrix of covariates
pi_min = sys.pi_min;                  % from theory
pi_max = sys.pi_max;                  % from theory
regressors = y(p+1:end);
if ~isfield(sys,'mu'); mu0 = mean(y)/2; else mu0 = sys.mu; end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ML_CVX %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Maximum Likelihood Estimator
    [mu_ML,theta_ML] = solve_ML_cvx(pi_min,pi_max,A,regressors,mu0);
disp('Maximum Likelihood estimate calculated')     
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ell_1 regularization Estimates %%%%%%%
% This is an implementation of ell_1 regularization in cvx
    cvx_quiet(true);
      cvx_begin 
        variable theta_sp(p)
        variable mu_sp
        minimize( (-(1-regressors)'*(1-mu_sp-A*theta_sp)-regressors'*log(mu_sp+A*theta_sp))/n  + lambda*norm(theta_sp,1))
        subject to
%         sum(theta_sp) <= 0.99;                                    % for stationarity
%         mu_sp >= pi_min;                                          % Avoids all zero processes
        norm(theta_sp,1)-sum(theta_sp) - 2*mu_sp <=  -2* pi_min ; % for theoretical guarantees
        norm(theta_sp,1)+sum(theta_sp) + 2*mu_sp <= 2*(pi_max);   % for theoretical guarantees
        mu_sp == mu0;
      cvx_end
      cvx_quiet(false);
disp('ell_1 regularized estimate calculated')     

%%%%%%%%%%%%%%%%%%%%%%%%%%%% OMP Estimates %%%%%%%%%%%%%%%%%   
% This is an implementation of the greedy method in cvx

    theta_OMP = zeros(p,1);
       s_star = sys.s_star;
       support = [];
       B = A;
       mu_OMP = mu0; % Initialization
       for s = 1:s_star
           l_grad = A'*(regressors./(mu_OMP+A*theta_OMP))-A'*((1-regressors)./(1-mu_OMP-A*theta_OMP));
           l_grad(support) = 0;
           [m,I] = max(abs( l_grad ));
           support = union(support,I);
           B = A(:,support);
           [mu_OMP, theta_OMP] = solve_ML_cvx(pi_min,pi_max,B,regressors,mu0);
           fprintf('Iteration %d of OMP completed \n',s);
%            mu0 = mu_OMP*(1-sum(theta_OMP));
       theta = zeros(p,1);
       theta(support) = theta_OMP;
       theta_OMP = theta;
%        stem(flipud(theta_OMP)); hold on; drawnow;
       end

disp('OMP estimate calculated')       
       
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Writing Outputs %%%%%%%%%%%%%%%%%%%%%%
n_test = length(sys.test_data)-p;
A_test = hankel(sys.test_data(1:n_test),sys.test_data(n_test:end-1));  % Matrix of covariates

sys_est.mu_ML = mu_ML;
sys_est.rates_ML = mu_ML+A*theta_ML;
sys_est.rates_ML_test = mu_ML+A_test*theta_ML;
sys_est.theta_ML = flipud(theta_ML);

sys_est.mu_sp = mu_sp;
sys_est.rates_sp = mu_sp+A*theta_sp;
sys_est.rates_sp_test = mu_sp+A_test*theta_sp;
sys_est.theta_sp = flipud(theta_sp);

sys_est.mu_OMP = mu_OMP;
sys_est.rates_OMP = mu_OMP+A*theta_OMP;
sys_est.rates_OMP_test = mu_OMP+A_test*theta_OMP;
sys_est.theta_OMP = flipud(theta_OMP);
%%%%%%%%%%%%%%%%%%%%%%%%%%% MSE for Synthetic %%%%%%%%%%%%%%%%%%%%%%
if isfield(sys,'theta_gt');
    sys_est.MSE_ML = norm(sys.theta_gt-sys_est.theta_ML)/norm(sys.theta_gt);
    sys_est.MSE_sp = norm(sys.theta_gt-sys_est.theta_sp)/norm(sys.theta_gt);
    sys_est.MSE_OMP = norm(sys.theta_gt-sys_est.theta_OMP)/norm(sys.theta_gt);
end

end



function [mu, theta] = solve_ML_cvx(pi_min,pi_max,A,regressors,varargin)
if nargin >4 
    mu0 = varargin{1};
end

[n,p] = size(A);
cvx_quiet(true);
  cvx_begin 
    variable theta(p)
    variable mu
    minimize( (-(1-regressors)'*(1-mu-A*theta)-regressors'*log(mu+A*theta))/n )
    subject to
%     sum(theta) <= 0.99;                              % for stationarity
%     mu >= pi_min;                                    % Avoids all zero processes
    norm(theta,1)-sum(theta) - 2*mu <=  -2* pi_min ; % for theoretical guarantees
    norm(theta,1)+sum(theta) + 2*mu <=  2*(pi_max);  % for theoretical guarantees
    if nargin>4
    mu == mu0;                                        % Only if mu is assumed to be known
    end
  cvx_end
  cvx_quiet(false);
  
end

