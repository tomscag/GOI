function oinf = goi_oinfo(X,biascorrect)
% Compute the O-information for the matrix X (nobs x nvartot)
% 
if nargin < 2; biascorrect = 'true'; end

X = goi_copnorm(X); % Copula normalization

% matrix size
[~, nvartot] = size(X);


H = goi_ent_g(X, biascorrect);

oinf = (nvartot - 2)*H;

for i = 1:nvartot
    oinf = oinf + goi_ent_g(X(:,i), biascorrect) - ...
           goi_ent_g(X(:,setdiff(1:nvartot,i)), biascorrect);
end


end