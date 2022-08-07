% TEST
clear
N       = 10^3;
nvartot = 5;
X = randn(N,nvartot);


%% COLLIDER : X2 --> X1 && X3 --> X1   (SYNERGY)
X(:,1) = 0.8*X(:,2) + 0.3*X(:,3) + 0.15*randn(N,1);

tic; out = goi_gradients(X,nvartot);  toc

out(1)
out(2)
% 0.5*log((2*pi*exp(1))^5*det(cov(X)))/log(2)


%% COUNFOUNDING X1 --> X2 && X1 --> X3 (REDUNDANCY)
X(:,2) = 0.8*X(:,1) + 0.15*randn(N,1);
X(:,3) = 0.8*X(:,1) + 0.15*randn(N,1);

tic, out = goi_gradients(X,nvartot);  toc

out(1)
out(2)

%% MEDIATOR : X3 --> X2 --> X1    (REDUNDANCY)

X(:,2) = 0.8*X(:,3) + 0.15*randn(N,1);
X(:,1) = 0.8*X(:,2) + 0.15*randn(N,1);



tic; out = goi_gradients(X,nvartot);  toc

out(1)