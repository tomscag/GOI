% Gradients of O-information in Ising system with three-spins interactions 

clear; clc

N = 8;  % Number of spins
% Create a triangular geometry for the three-spins interactions
T   = zeros(N,3);
T(1,:) = [1 2 3];
T(2,:) = [2 3 4];
T(3,:) = [3 4 5];
T(4,:) = [4 5 6];
T(5,:) = [5 6 7];
T(6,:) = [6 7 8];
T(7,:) = [7 8 1];
T(8,:) = [8 1 2];



%% Calcoliamo i gradienti
% Cerchiamo anche il picco delle fluttuazioni al variare della temperatura beta

nbeta  = 100;
beta   = linspace(0.05,3,nbeta);

grad_univ   = zeros(nbeta,N);    % Gradients O-information univariato
grad_pair   = zeros(nbeta,N*(N-1)/2);  % Gradients O-information bivariato
grad_trip   = zeros(nbeta,N*(N-1)*(N-2)/6);  % Gradients O-information trivariate
oinfo_loc   = zeros(nbeta,N*(N-1)/2);  % Local O-information Rosas

for w = 1:nbeta
    disp(beta(w))
    [P,spin,H,E] = ising_exact_triangular(N,T,beta(w));
    [dm,de] = ising_exact_fluctuations(P,E,spin);
    M(w) = dm;

    % Calcoliamo i gradienti
    prob = exp(-beta(w)*H);
    prob = prob/sum(prob(:));
 
    out =  gradients_exact(prob,3); % oinfo_Ising(prob);
    grad_univ(w,:) = out(1).O_grad;
    grad_pair(w,:) = out(2).O_grad;
    grad_trip(w,:) = out(3).O_grad;
    oinfo_loc(w,:) = out(2).O_loc;
end

[~, idx] = min(abs(diff(M)));
fprintf('Beta critico: %.4f\n',beta(idx))
% figure; plot(beta,M)


%%
clear pairs oiloc trip
pairs(:,1) = grad_pair(:,find(sum(out(2).index_var == [2,3],2) == 2));
pairs(:,2) = grad_pair(:,find(sum(out(2).index_var == [2,4],2) == 2));
pairs(:,3) = grad_pair(:,find(sum(out(2).index_var == [2,6],2) == 2));
pairs(:,4) = grad_pair(:,find(sum(out(2).index_var == [2,7],2) == 2));

oiloc(:,1) = oinfo_loc(:,find(sum(out(2).index_var == [2,3],2) == 2));
oiloc(:,2) = oinfo_loc(:,find(sum(out(2).index_var == [2,4],2) == 2));
oiloc(:,3) = oinfo_loc(:,find(sum(out(2).index_var == [2,6],2) == 2));
oiloc(:,4) = oinfo_loc(:,find(sum(out(2).index_var == [2,7],2) == 2));

trip(:,1) = grad_trip(:,find(sum(out(3).index_var == [1,2,3],2) == 3));
trip(:,2) = grad_trip(:,find(sum(out(3).index_var == [1,2,4],2) == 3));
trip(:,3) = grad_trip(:,find(sum(out(3).index_var == [1,2,5],2) == 3));
trip(:,4) = grad_trip(:,find(sum(out(3).index_var == [1,3,6],2) == 3));
trip(:,5) = grad_trip(:,find(sum(out(3).index_var == [1,3,5],2) == 3));
% trip(:,6) = grad_trip(:,find(sum(out(3).index_var == [2,5,7],2) == 3));


spin_pair = {'2-3','2-4','2-6','2-7'};
spin_trip = {'1-2-3','1-2-4','1-2-5','1-3-6','1-3-5'};



%% Figure
f_name   = 'Helvetica';
ax_fsize = 23;

% figure('Position',[488.0000   40.2000  873.8000  721.8000])

figure('Position',[295.4000  324.2000  445.6000  340.8000])
% subplot 322
plot(beta, M,'Color','k');   axis square; 
title('susceptibility','FontName',f_name)
xlabel('\beta','interpreter','tex','FontSize',ax_fsize); ylabel('$\chi$','Interpreter','latex','FontSize',ax_fsize)

figure('Position',[295.4000  324.2000  445.6000  340.8000]); 
% subplot 323
labels = {'spin 1'};
plot(beta, grad_univ(:,1),'k-'); axis square ;  legend(labels,'Location','best','FontSize',16,'Box','off');
title('first order gradients','FontName',f_name)
xlabel('\beta','interpreter','tex','FontSize',ax_fsize); ylabel('$\partial_i\Omega ({\bf s}^8)$','Interpreter','latex','FontSize',ax_fsize)



figure('Position',[295.4000  324.2000  445.6000  340.8000])
% subplot 324
plot(beta, pairs);   axis square; leg = legend(spin_pair,'Location','northeast','FontSize',13,'Box','off'); set(leg,'ItemTokenSize',[20 18]);
title('second order gradients','FontName',f_name)
xlabel('\beta','interpreter','tex','FontSize',ax_fsize); ylabel('$\partial_{ij}^2\Omega({\bf s}^8)$','Interpreter','latex','FontSize',ax_fsize)

figure('Position',[295.4000  324.2000  445.6000  340.8000])
% subplot 325
plot(beta, trip);   axis square; leg = legend(spin_trip,'Location','northeast','FontSize',13,'Box','off'); set(leg,'ItemTokenSize',[20 18]);
title('third order gradients','FontName',f_name)
xlabel('\beta','interpreter','tex','FontSize',ax_fsize); ylabel('$\partial_{ijk}^3\Omega({\bf s}^8)$','Interpreter','latex','FontSize',ax_fsize)


figure('Position',[295.4000  324.2000  445.6000  340.8000])
% subplot 326
plot(beta,oiloc);  legend(spin_pair,'Location','best','FontSize',13,'Box','off'); axis square
title('local O-information','FontName',f_name)
xlabel('\beta','interpreter','tex','FontSize',ax_fsize); ylabel('$I(s_i;s_j;{\bf s}^8_{-ij})$','Interpreter','latex','FontSize',ax_fsize-1)
% set(gca, 'Units', 'normalized','Position', [0.13, .18, 0.7750, 0.7195]); % [left bottom width height]
% set(get(gca,'title'),'Units', 'normalized','Position',[0.4 1.0 1 ])   % [left bottom width ]




