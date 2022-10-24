function [P,spin,E,EE] = ising_exact_probability(N,H,J,K,beta)
% FUNZIONE ACCESSORIA: serve solo per calcolare le probabilità, per tutto
% il resto usare ising_exact


% ising_analytic  Analytic solution for the 2-D Ising model in squared lattice
%                 with high-order interactions 
%   Input
%       N           -   Number of spin
%       J           -   Coupling matrix pairwise
%       H           -   External field
%       beta        -   Temperature
%   Output
%       P           -   Vettore di probabilità delle configurazioni di spin
%       spin        -   Lista le possibili configurazioni
%       E           -   2*2*...*2 probability mass function (2^N entries)

%   Per ottenere la probabilità basta fare exp(-beta*E) e normalizzare


% if nargin < 3
%     L    = sqrt(N);  % Per reticoli quadrati
%     J    = diag(ones(N-1,1),-1)   + diag(ones(N-1,1),1) + ...    % Coupling matrix
%         diag(ones(N-L,1),-L)   + diag(ones(N-L,1),L) + ...    % First neibor. interactions
%         diag(ones(L,1),N-L)    + diag(ones(L,1),-(N-L)) + ... % Boundary conditions;
%         diag(ones(1,1),N-1)    + diag(ones(1,1),-(N-1));
% end
if nargin < 2
    H    = zeros(1,N);      % External magnetic field
end

E    = zeros(2*ones(1,N));        % Store the energy of every configuration
spin = 2*de2bi(0:(2^N-1),N) - 1;  % Each row is a possible spin configuration

idx  = de2bi(0:(2^N-1),N) + 1;    % Map each spin configuration in a {1,2} array
% idx  = idx(1:height(idx)/2,:);  % Delete symmetry Z2 (flip all spins)

Z = 0;
% Iterate over all the spin configurations
for conf_k = 1:height(idx)

    e = - spin(conf_k,:)*J*spin(conf_k,:)'/2  - H*spin(conf_k,:)';
    
    % Somma il contributo di tutte le triplette
    for ii = 1:N
        for jj = 1:N
            for kk = 1:N
                e = e - K(ii,jj,kk)*spin(conf_k,ii)*spin(conf_k,jj)*spin(conf_k,kk)/6;
            end
        end
    end


    % Find the linear indices for prob
    sub      = num2cell(idx(conf_k,:));       % Subscript [row, col, page, ecc]
    idx_prob = sub2ind(size(E),sub{:});
%     E(idx_prob) = exp(-beta * E);
    E(idx_prob) = e;
    P(conf_k)  = exp(-beta*e);
    EE(conf_k) = e;
end
P = P/sum(P);
% Z    = sum(prob(:));
% % prob = prob/Z;
% Z    = log(Z);
end