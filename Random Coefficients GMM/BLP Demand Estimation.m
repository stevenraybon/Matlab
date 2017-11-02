clear variables 
close all 

rng(123);
%{  1  ,   2  ,  3   ,   4   ,   5   , 6 , 7,       8,             9,    10}
%{'mkt','time','firm','share','price','x','w','outside share', 'delta', 'x-j'}
data = csvread('hw4_data.csv',1,0);
data(:,8)=0;
data(:,9)=0;
data(:,10)=0;
len = length(data);

%For convenience, create individual vectors for each column used below
shares = data(:,4);
price = data(:,5);
x = data(:,6);
w = data(:,7);

%% OLS %%%%%%%%%%%%%%%%
%%Generate the outside shares for creation of y
for i=1:len
    
    if data(i,2)==1 && data(i,3)==1
        out_share = 1 - data(i,4) - data(i+1,4);
        data(i,8) = out_share;
        data(i+1,8) = out_share;
    elseif data(i,2)==2 && data(i,3)==1
        out_share2 = 1 - data(i,4) - data(i+1,4) - data(i+2,4);
        data(i,8) = out_share2;
        data(i+1,8) = out_share2;
        data(i+2,8) = out_share2;       
    else
    end
    
end

%Create y: ln(share_i/outside share) per each market
%%These are deltas
data(:,9) = log(data(:,4))-log(data(:,8));
delta = data(:,9);
%%Create X vector and add column of 1's for intercept (gamma) estimation
X_0 = data(:,5:6);
X_ones = ones([len,1]);
X = [X_ones,X_0];
y = delta;

%%OLS Estimators
beta_ols = X\y;
%or betas = inv(X'*X)*(X'*y)

%%OLS Standard Errors%%
%%Homoskedstic Std Errors

%Get the variance of the residuals
resid_ols = y - X*beta_ols;
resid2_ols = resid_ols.^2;
res_var_ols = mean(resid2_ols);

%Now calculate Var-Cov Matrix
VarMat_OLS = res_var_ols*inv(X'*X);
stderr_ols = [VarMat_OLS(1,1)^.5,VarMat_OLS(2,2)^.5,VarMat_OLS(3,3)^.5];

%% 2SLS %%
%%Generate X-j and add it to the data matrix
for i=1:len
    
    if data(i,2)==1 && data(i,3)==1
        data(i,10) = data(i+1,6);
        data(i+1,10) = data(i,6);
    elseif data(i,2)==2 && data(i,3)==1
        data(i,10) = sum([data(i+1,6), data(i+2,6)]);
        data(i+1,10) = sum([data(i,6), data(i+2,6)]);
        data(i+2,10) = sum([data(i,6), data(i+1,6)]);
    else
    end
    
end

W = [ones([len,1]),x,w,data(:,10)];
X_hat = W*inv(W'*W)*W'*X;
beta_2sls = X_hat\y;

% 2SLS Std Errors %
%%Homoskedastic Standard Errors

%Obtain the variance of the residuals
resid_2sls = y - X_hat*beta_2sls;
resid2_2sls = resid_2sls.^2;
res_var_2sls = mean(resid2_2sls);

%Now calculate Var-Cov Matrix
VarMat_2sls = res_var_2sls*inv(X_hat'*X_hat);
stderr_2sls = [VarMat_2sls(1,1)^.5,VarMat_2sls(2,2)^.5,VarMat_2sls(3,3)^.5];

%% Random Coefficient GMM Estimation %%%%%%%%

%Initial Guesses
sigma_0 = 0.5;

options_unc = optimoptions('fminunc','Algorithm','quasi-newton', ...
'SpecifyObjectiveGradient',false,'Display','iter','MaxIterations',5000,...
'MaxFunEvals',5000);

%Use fminunc to minimize the GMM Objective function using only sigma as
%input
sigma_star = fminunc(@(sigma)ObjFuncCalc(sigma,data),sigma_0,options_unc);

%Now use Sigma_hat to find other parameters
[J_star, gamma_star, alpha_star, beta_star] = ObjFuncCalc(sigma_star,data);

