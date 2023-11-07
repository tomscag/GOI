% Causal diagrams
clear; clc
N       = 10*10^3;
nvartot = 5;

addpath("~/gdrive/Dottorato_Matlab/ToolBox/robince-gcmi/matlab/")

%% COLLIDER : 2 --> 1 && 3 --> 1   (SYNERGY)
X = randn(N,nvartot);
X(:,1) = 0.8*X(:,2) + 0.3*X(:,3) + 0.15*randn(N,1);

out = goi_gradients(X,nvartot);  

fprintf("COLLIDER: 2 --> 1 && 3 --> 1  \n")
fprintf("First-order gradients: %2.2f %2.2f %2.2f %2.2f %2.2f \n\n",out(1).O_val)
% out(2).O_val
% 0.5*log((2*pi*exp(1))^5*det(cov(X)))/log(2)


%% COUNFOUNDING 1 --> 2 && 1 --> 3 (REDUNDANCY)
X = randn(N,nvartot);
X(:,2) = 0.8*X(:,1) + 0.15*randn(N,1);
X(:,3) = 0.8*X(:,1) + 0.15*randn(N,1);

out = goi_gradients(X,nvartot);  

fprintf("COUNFOUNDERS: 1 --> 2 && 1 --> 3\n")
fprintf("First-order gradients: %2.2f %2.2f %2.2f %2.2f %2.2f \n\n",out(1).O_val)
% out(2)

%% MEDIATOR : 3 --> 2 --> 1    (REDUNDANCY)
X      = randn(N,nvartot);
X(:,2) = 0.8*X(:,3) + 0.15*randn(N,1);
X(:,1) = 0.8*X(:,2) + 0.15*randn(N,1);

out = goi_gradients(X,nvartot); 

fprintf("MEDIATOR: X3 --> X2 --> X1 \n")
fprintf("First-order gradients: %2.2f %2.2f %2.2f %2.2f %2.2f \n\n",out(1).O_val)


%% COLLIDER + MEDIATOR  (Synergy + Redundancy)
% 2 --> 1 && 3 --> 1 
% 1 --> 4 && 4 --> 3

X = randn(N,nvartot);

X(:,1) = 0.8*X(:,2) + 0.3*X(:,3) + 0.15*randn(N,1);

X(:,4) = 0.8*X(:,1) + 0.15*randn(N,1);
% X(:,3) = 0.8*X(:,4) + 0.15*randn(N,1);


X(:,5) = 0.8*X(:,4) + 0.15*randn(N,1);

out_grad  = goi_gradients(X,nvartot); 
out_oinfo = goi_oinfo_all(X,nvartot);

% [out_grad(2).index_var out_grad(2).O_val']
[out_grad(3).index_var out_grad(3).O_val']
[out_oinfo(3).index_var out_oinfo(3).O_val']
%% OTHERS
X = randn(N,nvartot);
X(:,1) = 0.8*X(:,2) + 0.3*X(:,3) + 0.4*X(:,4) + 0.8*X(:,5) + 0.15*randn(N,1);
X(:,4) = 0.2*X(:,5) + 0.15*randn(N,1);
out = goi_gradients(X,nvartot); 

fprintf("OTHERS \n")
fprintf("First-order gradients: %2.2f %2.2f %2.2f %2.2f %2.2f \n\n",out(1).O_val)

figure; heatmap(table(out(2).index_var(:,1),out(2).index_var(:,2),out(2).O_val'),'Var1','Var2','ColorVariable','Var3')
title('Second-order gradients')
