function [J,gamma,alpha,beta] = ObjFuncCalc(sigma,data)
%This function takes in guesses for sigma, gamma, alpha and beta and then
%(1) simulates deltas via contraction mapping
%(2) uses 2SLS to instrument for p in constructing residuals
%(3) Evaluates the GMM Objective Function (e'W)(e'W)'

rng(123);
v = random('Normal',0,1,1000,1);
len = length(data);
sigma_guess = sigma;
delta_guess = random('Normal',1,.01,len,1);
x = data(:,6);
w = data(:,7);
W = [ones([len,1]),x,w,data(:,10)];
X_0 = data(:,5:6);
X_ones = ones([len,1]);
X = [X_ones,X_0];


%Contraction Mapping for Delta
t=.0001;
tol=t*ones([len,1]);
niter = 1;
itermax = 100;

while niter <= itermax
   delta_next = delta_guess + (log(data(:,4)) - log(simShareCalc(delta_guess,data,sigma_guess,v)));
   diff = log(data(:,4)) - log(simShareCalc(delta_guess,data,sigma_guess,v));
   if abs(diff) < tol
       break
   else
   delta_guess = delta_next;
   niter=niter+1;
   end
end
delta_final=delta_next;

%Price is correlated with Xi, so use IV(2SLS)
coef = (X'*W*W'*X)\(X'*W*W'*delta_final);
resid_gmm = delta_final - X*coef;
gamma = coef(1);
alpha = coef(2);
beta = coef(3);

%GMM criterion
J = (resid_gmm'*W)*(resid_gmm'*W)';
