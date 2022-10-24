function oinfo = goi_oinfo_all(ts,maxorder,biascorrect)
% Compute O-information for all the multiplets until a max order
% INPUT
%       ts          input (observations x variables), time series or static/behavioral data
%       maxorder    maximum order for the gradients  
%       biascorrect apply or not bias correction for entropy calculation

oinfo   = struct;
if nargin < 3; biascorrect = 'true'; end
    
% matrix size
[N, nvartot] = size(ts);

% Copula normalization
for i = 1:nvartot;  X(:,i) = goi_copnorm(ts(:,i)); end  


for isize = 3:maxorder  % Loop over gradient orders
    clear O_val
    C     = nchoosek(1:nvartot,isize); % List all the combinations
    ncomb = size(C,1);
    
    for icomb = 1:ncomb   % parfor        
        O_val(icomb) = goi_oinfo(X(:,C(icomb,:)));
    end

    oinfo(isize).O_val  = O_val;
    oinfo(isize).index_var   = C;
end