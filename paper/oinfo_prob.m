function out = oinfo_prob(prob)


%   INPUT 
%       prob        -   Probability mass function of X = (X_1, X_2, ... , X_N)

%   OUTPUT              out (struct)
%       OI          -   O-information    O = (N-2)H(X) + sum(  H(X_j) - H(X^N_{-j}))
%       TC          -   Total Correlation
%       DTC         -   Dual Total Correlation
%       d_oi        -   O-information delta due to the j-th variable  
%                       d_oi(j) = (2-N)I(X_j;X) + sum_i I(X_j;X_{-ij})
%       d_tc        -   Total Correlation delta due to the j-th
%                       variable   d_tc(j) = I(X_j;X)
%       d_dtc       -   Dual Total Correlation delta due to the j-th
%                       variable   d_dtc(j) = I(X_j;X) -  sum_i I(X_j;X_{-ij})
%       p_oi        -   Pairwise O-information
%       p_tc        -   Pairwise Total correlation
%       p_dtc       -   Pairwise Dual Total correlation
%       w_oi        -   O-information Local (Rosas) w_oi(i,j) = I(X_i;X_j; X_{-ij})

N       = length(size(prob));  % Number of spin

%%% 1) O-information %%%%

H      = sum(-prob(prob>0).*log2(prob(prob>0)));   % Joint entropy for N spin
h      = zeros(N,1);                               % Marginal entropies vector
H_not  = zeros(N,1);


for j = 1:N
    p_j        = squeeze(sum(prob,setdiff(1:N,j)));
    p_not_j    = squeeze(sum(prob,j));
    h(j)     = sum( -p_j(p_j>0).*log2(p_j(p_j>0)));
    H_not(j) = sum( -p_not_j(p_not_j>0).*log2(p_not_j(p_not_j>0)));
    clear p_j p_not_j
end

OI  = (N-2)*H  + sum(h - H_not); 
TC  = sum(h) - H;                   % Total Corelation
DTC = H - sum(H - H_not);




%%% 2) delta-Oinfo of each variable %%%

% d_O_jth = (1-N)I(X_j;X) + sum_i I(X_j;X_{-i})

d_oi   =  zeros(N,1);
d_tc   =  zeros(N,1);   % Variazione di Total Correlation
d_dtc  =  zeros(N,1);

for j = 1:N        
    mi         =  h(j) + H_not(j) - H;  % I(X_j; X_not_j)
    d_oi(j)    =  (1- (N-1) )*mi;
    d_tc(j)    =  mi;
    d_dtc(j)   =  (N-1)*mi;
    var_not_j  =  setdiff(1:N,j);   
    
    for i = 1:N-1
        
        p_not_ij    = squeeze(sum(prob,[var_not_j(i) j] ));
        h_not_ij    = sum( -p_not_ij(p_not_ij>0).*log2(p_not_ij(p_not_ij>0)));
        
        p_not_i    = squeeze(sum(prob,var_not_j(i))); 
        h_not_i =        sum( -p_not_i(p_not_i>0).*log2(p_not_i(p_not_i>0)));
        
        cmi =  h(j) + h_not_ij - h_not_i;
        d_oi(j) = d_oi(j) + cmi;
        d_dtc(j)   = d_dtc(j) - cmi;
        clear p_not_j_i
    end
end

%%% 3) Pairwise O-information 
p_oi   = zeros(N,N);
w_oi   = zeros(N,N);
p_tc   = zeros(N,N);   % Pairwise total correlation
p_dtc  = zeros(N,N);

for i = 1:N
    for j = i+1:N

        % Il primo pezzo: la variazione di dOj_i = O(A+i+j) - O(A+j) e quindi è
        % la variazione di O-info grazie alla variabile i e già lo abbiamo
        dO1 = d_oi(i);

        % Calcoliamo il secondo pezzo: la variazione di dO_i = O(A+i) - O(A)

        % I(Xi;X^n_{-ij})
        p_not_ij    = squeeze(sum(prob,[i j] ));
        h_not_ij    = sum( -p_not_ij(p_not_ij>0).*log2(p_not_ij(p_not_ij>0)));
        dO2 = (3-N)*(h(i)  +  h_not_ij - H_not(j) );

        var_not_ij   =  setdiff(1:N,[i j]);
        for k = 1:N-2
            p_not_ijk    = squeeze(sum(prob, [i j var_not_ij(k)]  ));
            h_not_ijk    = sum( -p_not_ijk(p_not_ijk>0).*log2(p_not_ijk(p_not_ijk>0)));

            p_not_jk    = squeeze(sum(prob,[j var_not_ij(k)] ));
            h_not_jk    = sum( -p_not_jk(p_not_jk>0).*log2(p_not_jk(p_not_jk>0)));

            dO2    = dO2 + ( h(i) + h_not_ijk - h_not_jk );
        end


        p_oi(i,j) = dO1 - dO2;  % Pairwise O-information
        p_oi(j,i) = p_oi(i,j);


        % Local O-information Rosas
        %         w_oi(i,j) = mi_gg(i,j) + mi_gg(i,A) - mi_gg(i,[A j]);
        p_ij    = squeeze(sum(prob, setdiff(1:N,[i j ])  ));
        h_ij    = sum( -p_ij(p_ij>0).*log2(p_ij(p_ij>0)));
        w1 = h(i) + h(j) - h_ij;   % I(X_i; X_j)


        p_not_ij    = squeeze(sum(prob,[i j]   ));
        h_not_ij    = sum( -p_not_ij(p_not_ij>0).*log2(p_not_ij(p_not_ij>0)));
        w2 = h(i) + h_not_ij - H_not(j);       % I(X_i; X^n_{-ij})

        w3 = h(i) + H_not(i) - H;     % I(X_i; X^n_{-i})
        w_oi(i,j) = w1 + w2 - w3;

        % Prova
        w_oi(i,j) = h(i) + h(j) + h_not_ij - H_not(j) - H_not(i) - h_ij + H;

        w_oi(j,i) = w_oi(i,j);


        % Ora calcoliamo il contributo della total correlation 
        p_tc(i,j) = H_not(j) - h_not_ij - H + H_not(i);
        p_tc(j,i) = p_tc(i,j);

        % E infine quello della dual total correlation
        p_dtc(i,j) = (N-1)*p_tc(i,j);

        var_not_ij   =  setdiff(1:N,[i j]);
        for k = 1:N-2
            p_not_ijk    = squeeze(sum(prob, [i j var_not_ij(k)]  ));
            h_not_ijk    = sum( -p_not_ijk(p_not_ijk>0).*log2(p_not_ijk(p_not_ijk>0)));

            p_not_jk    = squeeze(sum(prob,[j var_not_ij(k)] ));
            h_not_jk    = sum( -p_not_jk(p_not_jk>0).*log2(p_not_jk(p_not_jk>0)));

            p_not_ik    = squeeze(sum(prob,[i var_not_ij(k)] ));
            h_not_ik    = sum( -p_not_ik(p_not_ik>0).*log2(p_not_ik(p_not_ik>0)));

            p_dtc(i,j)  = p_dtc(i,j) - ( h_not_jk - h_not_ijk - H_not(var_not_ij(k)) + h_not_ik   );
            p_dtc(j,i)  = p_dtc(i,j);
        end

    end
end


out.OI     = OI;
out.TC     = TC;
out.DTC    = DTC;
out.d_oi   = d_oi;
out.d_tc   = d_tc;
out.d_dtc  = d_dtc;
out.p_oi   = p_oi;
out.p_tc   = p_tc;
out.p_dtc  = p_dtc;
out.w_oi   = w_oi;
out.H      = H;
out.h      = h;
out.H_not  = H_not;



end