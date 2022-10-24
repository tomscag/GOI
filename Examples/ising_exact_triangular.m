function [P,spin,Ematr,E] = ising_exact_triangular(N,T,beta)
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
%       E           -   Vettore di energie di ogni configurazione di spin
%       spin        -   Lista le possibili configurazioni
%       Ematr       -   2*2*...*2 probability mass function (2^N entries)


%   Per ottenere la probabilità basta fare exp(-beta*E) e normalizzare




Ematr = zeros(2*ones(1,N));      % Store the energy of every configuration
spin  = 2*de2bi(0:(2^N-1),N) - 1;  % Each row is a possible spin configuration
idx   = de2bi(0:(2^N-1),N) + 1;    % Map each spin configuration in a {1,2} array


% Iterate over all the spin configurations
for conf_k = 1:height(idx)

%     e = - spin(conf_k,:)*J*spin(conf_k,:)'/2  - H*spin(conf_k,:)';
    e = 0;
    for i = 1:N
        e = e - spin(conf_k,T(i,1))*spin(conf_k,T(i,2))*spin(conf_k,T(i,3));
    end

    % Find the linear indices for prob
    sub      = num2cell(idx(conf_k,:));       % Subscript [row, col, page, ecc]
    idx_prob = sub2ind(size(Ematr),sub{:});
    Ematr(idx_prob) = e;
    E(conf_k)  = e;
    P(conf_k)   = exp(-beta*e);    
end
P = P/sum(P);
end