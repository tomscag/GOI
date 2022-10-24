%% Gradients O-information analysis of US econometric time-series
clear; clc

%% Load dataset
load Serie_MacroEconomiche.mat
DataTable = DataQuaterly;  % DataQuaterly  DataMonthly
Description([17 21 24 28],:) = [];  
Data      = rmmissing(tick2ret(DataTable,'Method','continuous'));
Data      = Data(1:244,:);  % Escludi punti dopo 1 Jan, 2020 (Covid)
[nobs,nvars]  = size(Data);
labels    = DataTable.Properties.VariableNames;

%% Compute gradients
nsurr     = 250;  % Number of surrogates
maxorder  = 2;
goi = goi_gradients(Data{:,:},maxorder);


%% First-order gradients
fprintf('%s\t',Description(11,:))
fprintf('First-order gradients\n')
for i = 1:nvars
    fprintf('%s\t',Description(i+12,:))
    fprintf('%2.2f\n',goi(1).O_val(i))
end


%% Figure second-order gradients

% load results_FRED_500_surr.mat
node_size = 20;    % Node scaling factor
edge_size = 20;    % Edge scaling factor
f_size = 13;

h1 = figure('Position',[324.4000  324.2000  445.6000  340.8000]);
% X = O; X(isnan(X))=0; X(flagO == -1) = 0;

G = graph(goi(2).index_var(:,1), goi(2).index_var(:,2), goi(2).O_val,labels,'omitselfloops');
plot([NaN NaN], [NaN NaN], '-');  hold on; % Just to insert legend, this plot are invisible
plot([NaN NaN], [NaN NaN], 'r-'); hold on 
P = plot(G,'-ok','layout','circle','MarkerSize',7,... 
         'LineWidth',edge_size*abs(G.Edges.Weight)); axis off, axis square
highlight(P,'Edges',find( G.Edges.Weight < 0),'EdgeColor','blue') ;
highlight(P,'Edges',find( G.Edges.Weight > 0),'EdgeColor','red') 
title('gradients of O-information','FontSize',f_size)
legend({'negative','positive'},'Orientation','horizontal','Location','southoutside','Box','off')

set(gca, 'Units', 'normalized','Position', [0.30, .3, 0.5,0.5]); % [left bottom width height]
set(get(gca,'title'),'Units', 'normalized','Position',[0.4 1.2 1 ])   % [left bottom width ]


