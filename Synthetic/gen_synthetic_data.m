function sys = gen_synthetic_data( p, s, n, varargin )
    if nargin > 3 
        slack = varargin{1};
    else
        slack = n;
    end

theta  = zeros(p,1);
support = [sort(randsample(p,s))];
theta(support) = [0.1 0.1 0.1] ;
mu = .18;
y = zeros(1000*(n+p),1);
rates = zeros(1000*(n+p),1);

for t = p:1000*(n+p)
    rates(t+1) = mu+ y(t:-1:t-p+1)'*theta;
    y(t+1) = rand <= rates(t+1);
end

%%
sys.y = y(slack+1:slack+n+p);
sys.p = p;
sys.mu = mu;
sys.lambda = .1;
sys.theta_gt = theta;
sys.test_data = y(slack+n+p+1:end);
sys.pi_min = 0;
sys.pi_max = 0.49;
sys.s_star = s;

end

