function goi = goi_gradients(ts,maxorder,biascorrect)
% Compute low order O-infomation descriptor of data
% INPUT
%       ts          input (observations x variables), time series or static/behavioral data
%       maxorder    maximum order for the gradients  
%       biascorrect apply or not bias correction for entropy calculation

goi   = struct;
if nargin < 3; biascorrect = 'true'; end
    
% matrix size
[N, nvartot] = size(ts);

% Copula normalization
for i = 1:nvartot  X(:,i) = goi_copnorm(ts(:,i)); end  


for isize = 1:maxorder  % Loop over gradient orders
    clear O_val_univ
    C     = nchoosek(1:nvartot,isize); % List all the combinations
    ncomb = size(C,1);
    subset_index = dec2bin(0:2^isize-1) - '0';
    
    for icomb = 1:ncomb   % parfor
        O_val_univ(icomb) = 0;
        c = C(icomb,:);
        for n = 1:height(subset_index) % TO DO: this loop is quite inefficient as 'goi_oinfo' computes many times the same quantities
            idx = c(find(subset_index(n,:)));
            O_val_univ(icomb) = O_val_univ(icomb) + (-1)^numel(idx)*goi_oinfo(X(:,setdiff(1:nvartot,idx )));
        end
    end

    goi(isize).O_val_univ  = O_val_univ;
    goi(isize).index_var   = C;
end

