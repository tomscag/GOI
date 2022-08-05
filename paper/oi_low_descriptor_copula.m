function out = oi_low_descriptor_copula(data)

% Compute low order O-infomation descriptor of data


out   = struct;
[~, N] = size(data);


% Copula normalization
for i = 1:N  data(:,i) = goi_copnorm(data(:,i)); end  


% Univariate Descriptors
dT = nan(N,1);
dD = nan(N,1);
dO = nan(N,1);

% Pairwise Descriptors
ddT = nan(N,N);
ddD = nan(N,N);
ddO = nan(N,N);
ddI = nan(N,N);   % Interaction information

for i = 1:N

    X       = data(:,i);
    Z_not_i = data(:,setdiff(1:N,i));

    dT(i)   = goi_mi_gg(X, Z_not_i);  
    dD(i)   = (N-1)*dT(i);

    for k = 1:(N-1)
        dD(i) = dD(i) - goi_mi_gg(X, Z_not_i(:,setdiff(1:N-1,k)))  ;
    end
    dO(i) = dT(i) - dD(i);


    for j = i+1:N

        Y = data(:,j);
        Z = data(:,setdiff(1:N,[i j]));

        % Pairwise Total correlation  

        ddT(i,j) = goi_mi_gg(X, Z_not_i) - goi_mi_gg(X, Z);

        % Pairwise Dual Total correlation

        ddD(i,j) = (N-1)*ddT(i,j);
        
        for k = 1:(N-2)
            ddD(i,j) = ddD(i,j) - ( goi_mi_gg(X,[Z(:, setdiff(1:N-2,k)) Y])  - goi_mi_gg(X,Z(:, setdiff(1:N-2,k)) ) )  ;
        end


        ddO(i,j) = ddT(i,j) - ddD(i,j);


        % Interaction information w_ij = I(X;Y) + I(X;Z) -I(X;[Z Y])

        ddI(i,j) = goi_mi_gg(X,Y) + goi_mi_gg(X,Z) - goi_mi_gg(X,[Z Y]);
            


        ddT(j,i) = ddT(i,j);
        ddD(j,i) = ddD(i,j);
        ddO(j,i) = ddO(i,j);
        ddI(j,i) = ddI(i,j);

    end
end

out.dT  = dT;
out.dD  = dD;
out.dO  = dO;
out.ddT = ddT;
out.ddD = ddD;
out.ddO = ddO;
out.ddI = ddI;


end