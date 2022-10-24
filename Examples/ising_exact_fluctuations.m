function [M,E] = ising_exact_fluctuations(P,E,spin)

% Compute fluctuations of magnetization 

m = 0;   m1 = 0;
for i = 1:size(spin,1)
    m = m + P(i)*(mean(spin(i,:)).^2);
    m1 = m1 + P(i)*(mean(spin(i,:)));
end
M = m - m1^2;

% Calcoliamo le fluttuazioni dell' energia
e = 0;    e0 = 0;
for i = 1:size(spin,1)
    e = e + P(i)*E(i)^2;
    e0 = e0 + P(i)*E(i);
end

E = e -e0^2;


end