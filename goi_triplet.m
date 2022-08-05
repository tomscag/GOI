% DEPRECATED (use goi_oinfo)%

function out = goi_triplet(ts)

% Compute low order O-infomation descriptor of data

% ts = input (observations x variables), time series or static/behavioral data

out   = struct;


% matrix size
[N, nvartot] = size(ts);

% Copula normalization
for i = 1:nvartot  X(:,i) = goi_copnorm(ts(:,i)); end  

% TRIPLET
isize = 3; 

C     = nchoosek(1:nvartot,isize);
ncomb = size(C,1);


for icomb = 1:ncomb   % parfor
    O_val_trip(icomb) = goi_oinfo(X) - goi_oinfo(X(:,setdiff(1:nvartot,C(icomb,1)))) ...
        - goi_oinfo(X(:,setdiff(1:nvartot,C(icomb,2)))) -  goi_oinfo(X(:,setdiff(1:nvartot,C(icomb,3)))) ...
        + goi_oinfo(X(:,setdiff(1:nvartot,[C(icomb,1) C(icomb,2)]))) + goi_oinfo(X(:,setdiff(1:nvartot,[C(icomb,1) C(icomb,3)]))) ...
        + goi_oinfo(X(:,setdiff(1:nvartot,[C(icomb,2) C(icomb,3)]))) ...
        - goi_oinfo(X(:,setdiff(1:nvartot,[C(icomb,1) C(icomb,2) C(icomb,3)])));
end

out.O_val_trip = O_val_trip;


% PAIR
isize = 2; 
C     = nchoosek(1:nvartot,isize);
ncomb = size(C,1);
for icomb = 1:ncomb   % 
    O_val_pair(icomb) = goi_oinfo(X) ...
        - goi_oinfo(X(:,setdiff(1:nvartot,C(icomb,1)))) - goi_oinfo(X(:,setdiff(1:nvartot,C(icomb,2)))) + ...
        goi_oinfo(X(:,setdiff(1:nvartot,[C(icomb,1) C(icomb,2)])));
end
out.O_val_pair = O_val_pair;


% UNIVAR
isize = 1; 

C     = nchoosek(1:nvartot,isize);
ncomb = size(C,1);
for icomb = 1:ncomb   % 
    O_val_univ(icomb) = goi_oinfo(X) - goi_oinfo(X(:,setdiff(1:nvartot,C(icomb,1))));
end
out.O_val_univ = O_val_univ;



end