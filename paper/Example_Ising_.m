%%% Toy Model Ising Frustrato
%   Costruiamo una geometria esagonale con link -1 e un settimo spin al
%   centro collegato agli altri con spin +1
clear; clc

geometria = 'esagonale';   % esagonale  pentagonale solo_periferici

switch geometria
    case 'esagonale'
        N = 7;
        B      = [0,0,0,0,0,0,0];        % External magnetic field
        
        J = zeros(N,N);                  % Geometria esagonale
        J(1,:) = [0,-1,0,0,0,-1,+1]; 
        J(2,:) = [-1,0,-1,0,0,0,+1]; 
        J(3,:) = [0,-1,0,-1,0,0,+1]; 
        J(4,:) = [0,0,-1,0,-1,0,+1]; 
        J(5,:) = [0,0,0,-1,0,-1,+1]; 
        J(6,:) = [-1,0,0,0,-1,0,+1]; 
        J(7,:) = [1,1,1,1,1,1,0]; 

    case 'pentagonale'
        N      = 6;
        J      = zeros(N,N);           % Geometria pentagonale
        B      = [0,0,0,0,0,0];        % External magnetic field
        J(1,:) = [0,-1,0,0,-1,+1]; 
        J(2,:) = [-1,0,-1,0,0,+1]; 
        J(3,:) = [0,-1,0,-1,0,+1]; 
        J(4,:) = [0,0,-1,0,-1,+1]; 
        J(5,:) = [-1,0,0,-1,0,+1]; 
        J(6,:) = [1,1,1,1,1,0]; 
    case 'solo_periferici'
        N = 6;
        B      = [0,0,0,0,0,0];        % External magnetic field
        
        J = zeros(N,N);                  % Geometria esagonale
        J(1,:) = [0,-1,0,0,0,-1]; 
        J(2,:) = [-1,0,-1,0,0,0]; 
        J(3,:) = [0,-1,0,-1,0,0]; 
        J(4,:) = [0,0,-1,0,-1,0]; 
        J(5,:) = [0,0,0,-1,0,-1]; 
        J(6,:) = [-1,0,0,0,-1,0];       
end



nbeta       = 60;
beta        = linspace(0.05,3,nbeta);
O           = zeros(nbeta,1);    % O-information
TC          = zeros(nbeta,1);    % Total Correlation
DTC         = zeros(nbeta,1);    % Dual Total Correlation
delta_O     = zeros(nbeta,N);    % Delta O-information
H           = zeros(nbeta,1);    % Joint Entropy
H_j         = zeros(nbeta,N);    % Entropy of each spin
H_not_j     = zeros(nbeta,N);    % Conditional Entropy
p_oi        = zeros(nbeta,N,N);  % Pairwise O-information
w_oi        = zeros(nbeta,N,N);  % Local O-information Rosas
p_tc        = zeros(nbeta,N,N);  % Pairwise O-information
p_dtc       = zeros(nbeta,N,N);  % Pairwise O-information

E = ising_exact(N,B,J);      % Energie delle 2^N configurazioni

%%
for w = 1:nbeta
    disp(beta(w))
    prob = exp(-beta(w)*E);
    prob = prob/sum(prob(:));
    out =  oinfo_prob(prob); % oinfo_Ising(prob); 

    O(w)   = out.OI;       delta_O(w,:) = out.d_oi;    w_oi(w,:,:) = out.w_oi;
    TC(w)  = out.TC;      d_tc(w,:)    = out.d_tc;
    DTC(w) = out.DTC;     d_dtc(w,:)   = out.d_dtc; 
    p_oi(w,:,:)  = out.p_oi;
    p_tc(w,:,:)  = out.p_tc; 
    p_dtc(w,:,:) = out.p_dtc;
        
    H(w)         = out.H; 
    H_j(w,:)     = out.h; 
    H_not_j(w,:) = out.H_not;

end

clear p_oi
p_oi = p_tc - p_dtc;

% Plot
set(0,'defaultAxesFontSize', 14); 
set(0,'DefaultLineLineWidth',1.25)

clear A B
switch geometria
    case 'esagonale'
        A(:,1) = squeeze(p_oi(:,1,2));  B(:,1) = squeeze(w_oi(:,1,2)); 
        A(:,2) = squeeze(p_oi(:,1,3));  B(:,2) = squeeze(w_oi(:,1,3));
        A(:,3) = squeeze(p_oi(:,1,4));  B(:,3) = squeeze(w_oi(:,1,4));
        A(:,4) = squeeze(p_oi(:,1,7));  B(:,4) = squeeze(w_oi(:,1,7));
        T(:,1) = squeeze(p_tc(:,1,2));  D(:,1) = squeeze(p_dtc(:,1,2)); 
        T(:,2) = squeeze(p_tc(:,1,3));  D(:,2) = squeeze(p_dtc(:,1,3));
        T(:,3) = squeeze(p_tc(:,1,4));  D(:,3) = squeeze(p_dtc(:,1,4));
        T(:,4) = squeeze(p_tc(:,1,7));  D(:,4) = squeeze(p_dtc(:,1,7));
%         spin_pair = {'1-2','1-3','1-4','1-7'};
        spin_pair = {'2-3','2-4','2-5','1-2'};
    case 'pentagonale'
        A(:,1) = squeeze(p_oi(:,1,2));  B(:,1) = squeeze(w_oi(:,1,2));
        A(:,2) = squeeze(p_oi(:,1,3));  B(:,2) = squeeze(w_oi(:,1,3));
        A(:,3) = squeeze(p_oi(:,1,6));  B(:,3) = squeeze(w_oi(:,1,6));
        T(:,1) = squeeze(p_tc(:,1,2));  D(:,1) = squeeze(p_dtc(:,1,2)); 
        T(:,2) = squeeze(p_tc(:,1,3));  D(:,2) = squeeze(p_dtc(:,1,3));
        T(:,3) = squeeze(p_tc(:,1,6));  D(:,3) = squeeze(p_dtc(:,1,6));
        spin_pair = {'1-2','1-3','1-6'};
    case 'solo_periferici'
        A(:,1) = squeeze(p_oi(:,1,2));  B(:,1) = squeeze(w_oi(:,1,2));
        A(:,2) = squeeze(p_oi(:,1,3));  B(:,2) = squeeze(w_oi(:,1,3));
        A(:,3) = squeeze(p_oi(:,1,4));  B(:,3) = squeeze(w_oi(:,1,4));
        T(:,1) = squeeze(p_tc(:,1,2));  D(:,1) = squeeze(p_dtc(:,1,2)); 
        T(:,2) = squeeze(p_tc(:,1,3));  D(:,2) = squeeze(p_dtc(:,1,3));
        T(:,3) = squeeze(p_tc(:,1,4));  D(:,3) = squeeze(p_dtc(:,1,4));
        spin_pair = {'1-2','1-3','1-4'};
end


%
if strcmp(geometria,'pentagonale'); spin_c = 6; else; spin_c = 7; end
labels = {'spin 7','spin 1'};

%% Figura Physical Review Letter
f_name   = 'Helvetica';
ax_fsize = 23;
figure('Position',[295.4000  324.2000  445.6000  340.8000])
plot(beta, delta_O(:,1),'k-.', beta, delta_O(:,spin_c),'k-'); axis square ;  legend(labels,'Location','best','FontSize',16,'Box','off');
title('first order gradients','FontName',f_name)
xlabel('\beta','interpreter','tex','FontSize',ax_fsize); ylabel('$\partial_i\Omega ({\bf s}^7)$','Interpreter','latex','FontSize',ax_fsize)

figure('Position',[295.4000  324.2000  445.6000  340.8000])
% plot(beta, A(:,1),'k:',beta, A(:,2),'k--',beta, A(:,3),'k.-',beta, A(:,4),'k-');   axis square
plot(beta, A);   axis square; leg = legend(spin_pair,'Location','northeast','FontSize',13,'Box','off'); set(leg,'ItemTokenSize',[20 18]);
title('second order gradients','FontName',f_name)
xlabel('\beta','interpreter','tex','FontSize',ax_fsize); ylabel('$\partial_{ij}^2\Omega({\bf s}^7)$','Interpreter','latex','FontSize',ax_fsize)


figure('Position',[295.4000  324.2000  445.6000  340.8000])
% plot(beta, B(:,1),'k:',beta, B(:,2),'k--',beta, B(:,3),'k.-',beta, B(:,4),'k-');  legend(spin_pair,'Location','best'); axis square
clf
plot(beta, B); axis square;  leg= legend(spin_pair,'Location','northwest','FontSize',13,'Box','off'); set(leg,'ItemTokenSize',[20 18]);
title('local O-information','FontName',f_name)
xlabel('\beta','interpreter','tex','FontSize',ax_fsize); ylabel('$I(s_i;s_j;{\bf s}^7_{-ij})$','Interpreter','latex','FontSize',ax_fsize-1)
% set(gca, 'Units', 'normalized','Position', [0.13, .18, 0.7750, 0.7195]); % [left bottom width height]
% set(get(gca,'title'),'Units', 'normalized','Position',[0.4 1.0 1 ])   % [left bottom width ]







