
%% 
clear
load Serie_MacroEconomiche.mat
DataTable = DataQuaterly;  % DataQuaterly  DataMonthly
Description([17 21 24 28],:) = [];  
% DataTable = removevars(DataTable,{'GS10','M1SL','M2SL','TB3MS'});

nsurr = 250;

Data      = rmmissing(tick2ret(DataTable,'Method','continuous'));
Data      = Data(1:244,:);  % Escludi punti dopo 1 Jan, 2020 (Covid)
[nobs,N]  = size(Data);
labels    = DataTable.Properties.VariableNames;
% for i=1:14; names{i} = extractBetween(lab(i,:),'(',')'); names{i} = names{i}(1); end


out = oi_low_descriptor_copula(Data{:,:});
O   = out.ddO; O(isnan(O)) = 0;
T   = out.ddT; T(isnan(T)) = 0;
D   = out.ddD; D(isnan(D)) = 0;
I   = out.ddI; I(isnan(I)) = 0;
dO  = out.dO;

C = corr(Data{:,:});
[PC, p] = partialcorr(Data{:,:}); 

for s = 1:nsurr
    disp(s)
    p = datasample(1:nobs,nobs,'Replace',true); %  Sampling with replacemente of Rows
    XS = Data{p,:};
    outS = oi_low_descriptor_copula(XS);
    Osur(s,:,:)    = outS.ddO; 
    Tsur(s,:,:)    = outS.ddT; 
    Dsur(s,:,:)    = outS.ddD;
    Isur(s,:,:)    = outS.ddI; 
    dOsur(s,:)     = outS.dO;
    PCsur(s,:,:)   = partialcorr(XS); 
    Csurr(s,:,:)   = corrcoef(XS);
    tsurr(s,:)     = outS.dT;
    dsurr(s,:)     = outS.dD;
    osurr(s,:)     = outS.dO;
end

% Disp Result
% fprintf('Method: %s\n',method)
fprintf('%s\t',Description(11,:))
fprintf('Delta TC \tDelta DTC\t Delta O\n')
for i = 1:N
    fprintf('%s\t',Description(i+12,:))
    fprintf('%2.2f\t\t',out.dT(i))
    fprintf('%2.2f\t\t',out.dD(i))
    fprintf('%2.2f\n',out.dO(i))
end

% Test significatività
alpha= 0.05;
for i = 1:N
    flago5(i) = prod(sign(prctile(osurr(:,i),[100*0.05/2,100*(1-0.05/2)])));
    for j = 1:N
        flagO(i,j) = prod(sign(prctile(Osur(:,i,j),[100*alpha/2,100*(1-alpha/2)])));
        flagI(i,j) = prod(sign(prctile(Isur(:,i,j),[100*alpha/2,100*(1-alpha/2)])));
        flagP(i,j) = prod(sign(prctile(PCsur(:,i,j),[100*alpha/2,100*(1-alpha/2)])));
        flagC(i,j) = prod(sign(prctile(Csurr(:,i,j),[100*alpha/2,100*(1-alpha/2)])));
    end
end
fprintf('Correlation between dO and w: \t%1.3f\n',corr(O(:),I(:)))
% figure; scatter(O(:),W(:)), axis square

for i = 1:14
    dObounds(i,:) = prctile(osurr(:,i),[100*0.05/2,100*(1-0.05/2)]);
end



%% Figura 

% load results_FRED_500_surr.mat
node_size = 20;      % For node marker size
edge_size = 20;
f_size = 13;

% O-information
h1 = figure('Position',[324.4000  324.2000  445.6000  340.8000]);
X = O; X(isnan(X))=0; X(flagO == -1) = 0;
dX = dO;
G    = graph(X,labels,'omitselfloops');
plot([NaN NaN], [NaN NaN], '-');hold on; plot([NaN NaN], [NaN NaN], 'r-'); hold on % Just to put legend, this plot are invisible
P = plot(G,'-ok','layout','circle','MarkerSize',7,... %node_size*abs(dX),...
         'LineWidth',edge_size*abs(G.Edges.Weight)); axis off, axis square
highlight(P,'Edges',find( G.Edges.Weight < 0),'EdgeColor','blue') ;
highlight(P,'Edges',find( G.Edges.Weight > 0),'EdgeColor','red') 
title('gradients of O-information','FontSize',f_size)
legend({'negative','positive'},'Orientation','horizontal','Location','southoutside','Box','off')

set(gca, 'Units', 'normalized','Position', [0.30, .3, 0.5,0.5]); % [left bottom width height]
set(get(gca,'title'),'Units', 'normalized','Position',[0.4 1.2 1 ])   % [left bottom width ]

h1 = figure('Position',[324.4000  324.2000  445.6000  340.8000]);
node_size = 20;      % For node marker size
edge_size = 20;
% Interaction Information
X = I; X(isnan(X))=0; X(flagI == -1) = 0;
G    = graph(X,labels,'omitselfloops');
plot([NaN NaN], [NaN NaN], '-');hold on; plot([NaN NaN], [NaN NaN], 'r-'); hold on % Just to put legend, this plot are invisible
P = plot(G,'-ok','layout','circle','MarkerSize',7,...
         'LineWidth',edge_size*abs(G.Edges.Weight)); axis off, axis square
highlight(P,'Edges',find( G.Edges.Weight < 0),'EdgeColor','blue') ;
highlight(P,'Edges',find( G.Edges.Weight > 0),'EdgeColor','red') 

title('local O-information','FontSize',f_size) % I(X_i;X_j;\bf{X}^n_{-ij})
%legend({'negative','positive'},'Orientation','horizontal','Location','southoutside','Box','off')
c = colorbar('southoutside'); c.TickLabels  = {'-','0','+'};
set(gca, 'Units', 'normalized','Position', [0.30, .3, 0.5,0.5]); % [left bottom width height]
set(get(gca,'title'),'Units', 'normalized','Position',[0.4 1.2 1 ])   % [left bottom width ]



col = colormap(brewermap(2,'PrGn'));
node_size = 10;      % For node marker size
edge_size = 10;
h1 = figure('Position',[324.4000  324.2000  445.6000  340.8000]);
%  Correlation
X = PC; X(isnan(X))=0; X(flagP == -1) = 0; %  X(p>alpha) = 0;
G    = graph(X,labels,'omitselfloops');
plot([NaN NaN], [NaN NaN], '-','Color',col(1,:));hold on; plot([NaN NaN], [NaN NaN], 'r-','Color',col(2,:)); hold on % Just to put legend, this plot are invisible
P = plot(G,'-ok','layout','circle','MarkerSize',7,...
         'LineWidth',edge_size*abs(G.Edges.Weight)); axis off, axis square
highlight(P,'Edges',find( G.Edges.Weight < 0),'EdgeColor',col(1,:)) ; % Green
highlight(P,'Edges',find( G.Edges.Weight > 0),'EdgeColor',col(2,:)) % Orange

title('partial correlation','FontSize',f_size)
legend({'negative','positive'},'Orientation','horizontal','Location','southoutside','Box','off')

set(gca, 'Units', 'normalized','Position', [0.30, .3, 0.5,0.5]); % [left bottom width height]
set(get(gca,'title'),'Units', 'normalized','Position',[0.4 1.2 1 ])   % [left bottom width ]

