%%% Toy Model Ising Frustrato
%   Costruiamo una geometria esagonale con link -1 e un settimo spin al
%   centro collegato agli altri con spin +1
clear; clc
addpath('C:\Users\tomsc\Dropbox\Matlab\Dottorato Matlab\Projects\my_Code\GHOI')
N = 7;
B      = [0,0,0,0,0,0,0];        % External magnetic field

% Create hexagonal geometry
J = zeros(N,N);   

J(1,:) = [0,1,1,1,1,1,1];
J(2,:) = [1,0,-1,0,0,0,-1];
J(3,:) = [1,-1,0,-1,0,0,0];
J(4,:) = [1,0,-1,0,-1,0,0];
J(5,:) = [1,0,0,-1,0,-1,0];
J(6,:) = [1,0,0,0,-1,0,-1];
J(7,:) = [1,-1,0,0,0,-1,0];



% Uncomment to flip all the couplings
 J      = - J;

K         = zeros(N,N,N);
% Uncomment to add three-spins interactions with coupling k
% k         = 0;             % Coupling 
% K(1,2,3)  = k;
% K(1,3,4)  = k;
% K(1,4,5)  = k;
% K(1,5,6)  = k;
% K(1,6,7)  = k;
% K(1,7,2)  = k;
% K         = K + permute(K,[1 3 2]) + permute(K,[2 1 3]) + ...  
%     permute(K,[2 3 1]) + permute(K,[3 1 2]) + permute(K,[3 2 1]); % Simmetrizzo


nbeta       = 60;
beta        = linspace(0.05,3,nbeta);   % Inverse temperatures list

grad_univ   = zeros(nbeta,N);    % Gradients O-information univariato
grad_pair   = zeros(nbeta,N*(N-1)/2);  % Gradients O-information bivariato
grad_trip   = zeros(nbeta,N*(N-1)*(N-2)/6);  % Gradients O-information trivariate
oinfo_loc   = zeros(nbeta,N*(N-1)/2);  



% Energie delle 2^N configurazioni (K parametro interazioni a tre spins)
[Ematr,E,spin] = ising_exact_hoi(N,B,J,K);  
%% Compute gradients


for w = 1:nbeta
    disp(beta(w))

    P = exp(-beta(w)*E); P = P/sum(P);
    [dm,de] = ising_exact_fluctuations(P,E,spin);
    M(w) = dm;

    
    prob = exp(-beta(w)*Ematr);
    prob = prob/sum(prob(:));
 
    out =  gradients_exact(prob,3); % oinfo_Ising(prob);
    grad_univ(w,:) = out(1).O_grad;
    grad_pair(w,:) = out(2).O_grad;
    grad_trip(w,:) = out(3).O_grad;
    oinfo_loc(w,:) = out(2).O_loc;
end


%% 
clear pairs trip oiloc
set(0,'defaultAxesFontSize', 14);
set(0,'DefaultLineLineWidth',1.25)


pairs(:,1) = grad_pair(:,find(sum(out(2).index_var == [2,3],2) == 2));
pairs(:,2) = grad_pair(:,find(sum(out(2).index_var == [2,4],2) == 2));
pairs(:,3) = grad_pair(:,find(sum(out(2).index_var == [2,5],2) == 2));
pairs(:,4) = grad_pair(:,find(sum(out(2).index_var == [1,2],2) == 2));

oiloc(:,1) = oinfo_loc(:,find(sum(out(2).index_var == [2,3],2) == 2));
oiloc(:,2) = oinfo_loc(:,find(sum(out(2).index_var == [2,4],2) == 2));
oiloc(:,3) = oinfo_loc(:,find(sum(out(2).index_var == [2,6],2) == 2));
oiloc(:,4) = oinfo_loc(:,find(sum(out(2).index_var == [1,2],2) == 2));

trip(:,1) = grad_trip(:,find(sum(out(3).index_var == [1,2,3],2) == 3));
trip(:,2) = grad_trip(:,find(sum(out(3).index_var == [1,2,4],2) == 3));
trip(:,3) = grad_trip(:,find(sum(out(3).index_var == [1,2,5],2) == 3));
trip(:,4) = grad_trip(:,find(sum(out(3).index_var == [2,3,4],2) == 3));
trip(:,5) = grad_trip(:,find(sum(out(3).index_var == [2,3,5],2) == 3));
trip(:,6) = grad_trip(:,find(sum(out(3).index_var == [2,4,6],2) == 3));


spin_pair = {'2-3','2-4','2-5','1-2'};
spin_trip = {'1-2-3','1-2-4','1-2-5','2-3-4','2-3-5','2-4-6'};




%% Figure
f_name   = 'Helvetica';
ax_fsize = 23;

figure('Position',[488.0000   40.2000  873.8000  721.8000])

% figure('Position',[295.4000  324.2000  445.6000  340.8000])
subplot 322
plot(beta, M,'Color','k');   axis square; 
title('susceptibility','FontName',f_name)
xlabel('\beta','interpreter','tex','FontSize',ax_fsize); ylabel('$\chi$','Interpreter','latex','FontSize',ax_fsize)

% figure('Position',[295.4000  324.2000  445.6000  340.8000]); 
subplot 323
labels = {'spin 1','spin 7'}; 
plot(beta, grad_univ(:,1),'k-', beta, grad_univ(:,7),'k--'); axis square ;  legend(labels,'Location','best','FontSize',16,'Box','off');
title('first order gradients','FontName',f_name)
xlabel('\beta','interpreter','tex','FontSize',ax_fsize); ylabel('$\partial_i\Omega ({\bf s}^7)$','Interpreter','latex','FontSize',ax_fsize)



% figure('Position',[295.4000  324.2000  445.6000  340.8000])
subplot 324
plot(beta, pairs);   axis square; leg = legend(spin_pair,'Location','northeast','FontSize',13,'Box','off'); set(leg,'ItemTokenSize',[20 18]);
title('second order gradients','FontName',f_name)
xlabel('\beta','interpreter','tex','FontSize',ax_fsize); ylabel('$\partial_{ij}^2\Omega({\bf s}^7)$','Interpreter','latex','FontSize',ax_fsize)

% figure('Position',[295.4000  324.2000  445.6000  340.8000])
subplot 325
plot(beta, trip);   axis square; leg = legend(spin_trip,'Location','northeast','FontSize',13,'Box','off'); set(leg,'ItemTokenSize',[20 18]);
title('third order gradients','FontName',f_name)
xlabel('\beta','interpreter','tex','FontSize',ax_fsize); ylabel('$\partial_{ijk}^3\Omega({\bf s}^7)$','Interpreter','latex','FontSize',ax_fsize)


% figure('Position',[295.4000  324.2000  445.6000  340.8000])
subplot 326
plot(beta,oiloc);  legend(spin_pair,'Location','best','FontSize',16,'Box','off'); axis square
title('local O-information','FontName',f_name)
xlabel('\beta','interpreter','tex','FontSize',ax_fsize); ylabel('$I(s_i;s_j;{\bf s}^7_{-ij})$','Interpreter','latex','FontSize',ax_fsize-1)
% set(gca, 'Units', 'normalized','Position', [0.13, .18, 0.7750, 0.7195]); % [left bottom width height]
% set(get(gca,'title'),'Units', 'normalized','Position',[0.4 1.0 1 ])   % [left bottom width ]








%% Figura Physical Review Letter
% f_name   = 'Helvetica';
% ax_fsize = 23;
% 
% figure('Position',[295.4000  324.2000  445.6000  340.8000])
% plot(beta, grad_univ(:,7),'k-.', beta, grad_univ(:,spin_c),'k-'); axis square ;  legend(labels,'Location','best','FontSize',16,'Box','off');
% title('first order gradients','FontName',f_name)
% xlabel('\beta','interpreter','tex','FontSize',ax_fsize); ylabel('$\partial_i\Omega ({\bf s}^7)$','Interpreter','latex','FontSize',ax_fsize)
% 
% figure('Position',[295.4000  324.2000  445.6000  340.8000])
% plot(beta, pairs);   axis square; leg = legend(spin_pair,'Location','northeast','FontSize',13,'Box','off'); set(leg,'ItemTokenSize',[20 18]);
% title('second order gradients','FontName',f_name)
% xlabel('\beta','interpreter','tex','FontSize',ax_fsize); ylabel('$\partial_{ij}^2\Omega({\bf s}^7)$','Interpreter','latex','FontSize',ax_fsize)
% 
% figure('Position',[295.4000  324.2000  445.6000  340.8000])
% plot(beta, trip);   axis square; leg = legend(spin_trip,'Location','northeast','FontSize',13,'Box','off'); set(leg,'ItemTokenSize',[20 18]);
% title('third order gradients','FontName',f_name)
% xlabel('\beta','interpreter','tex','FontSize',ax_fsize); ylabel('$\partial_{ijk}^3\Omega({\bf s}^7)$','Interpreter','latex','FontSize',ax_fsize)
% 
% figure('Position',[295.4000  324.2000  445.6000  340.8000])
% plot(beta, oiloc); axis square;  leg= legend(spin_pair,'Location','northwest','FontSize',13,'Box','off'); set(leg,'ItemTokenSize',[20 18]);
% title('local O-information','FontName',f_name)
% xlabel('\beta','interpreter','tex','FontSize',ax_fsize); ylabel('$I(s_i;s_j;{\bf s}^7_{-ij})$','Interpreter','latex','FontSize',ax_fsize-1)
% % set(gca, 'Units', 'normalized','Position', [0.13, .18, 0.7750, 0.7195]); % [left bottom width height]
% % set(get(gca,'title'),'Units', 'normalized','Position',[0.4 1.0 1 ])   % [left bottom width ]
% 






