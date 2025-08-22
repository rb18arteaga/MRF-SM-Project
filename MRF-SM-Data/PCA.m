clear all
clc
%% Read average results
dataREAD = readtable('average_results.csv');

name = dataREAD.OriginalVariableNames;  
values = dataREAD.Var1;  

% Initialize structures for each MRF configuration
MRF1 = struct();
MRF2 = struct();
MRF3 = struct();
MRF4 = struct();
SM = struct();

% Organize data by MRF configuration
for i = 1:length(name)
    current_name = name{i};
    current_value = values(i);
    
    if startsWith(current_name, 'MRF1_')
        field = strrep(current_name, 'MRF1_', ''); 
        MRF1.(field) = current_value;
    elseif startsWith(current_name, 'MRF2_')
        field = strrep(current_name, 'MRF2_', '');
        MRF2.(field) = current_value;
    elseif startsWith(current_name, 'MRF3_')
        field = strrep(current_name, 'MRF3_', '');
        MRF3.(field) = current_value;
    elseif startsWith(current_name, 'MRF4_')
        field = strrep(current_name, 'MRF4_', '');
        MRF4.(field) = current_value;
    elseif startsWith(current_name, 'SM_')
        field = strrep(current_name, 'SM_', '');
        SM.(field) = current_value;
    else
        warning('Unclassified variable: %s', current_name);
    end
end

% Convert structures to tables
table_MRF1 = struct2table(MRF1);
table_MRF2 = struct2table(MRF2);
table_MRF3 = struct2table(MRF3);
table_MRF4 = struct2table(MRF4);
table_SM = struct2table(SM);

disp('MRF1 data:');
disp(table_MRF1);

%% Speed comparison: MRF2, MRF3, MRF1, MRF4, SM
Method = ["MRF"; "MRF"; "MRF"; "MRF";"SM"];
w_omega = [31.7 , 44.1, 69.7, 99.5, 69.7]';
h = [
    double(table_MRF2{1, 11}),  
    double(table_MRF3{1, 11}),  
    double(table_MRF1{1, 11}),  
    double(table_MRF4{1, 11}), 
    double(table_SM{1, 11})
];

T_max = [
    double(table_MRF2{1, 5}),  
    double(table_MRF3{1, 5}),  
    double(table_MRF1{1, 5}),  
    double(table_MRF4{1, 5}), 
    double(table_SM{1, 5})
];

Q = [
    double(table_MRF2{1, 6}),  
    double(table_MRF3{1, 6}),  
    double(table_MRF1{1, 6}),  
    double(table_MRF4{1, 6}), 
    double(table_SM{1, 6})
];

Tavg = [
    double(table_MRF2{1, 8}),  
    double(table_MRF3{1, 8}),  
    double(table_MRF1{1, 8}),  
    double(table_MRF4{1, 8}), 
    double(table_SM{1, 8})
];

Cd_t = [
    double(table_MRF2{1, 1}),  
    double(table_MRF3{1, 1}),  
    double(table_MRF1{1, 1}),  
    double(table_MRF4{1, 1}), 
    double(table_SM{1, 1})
];

V_vanes = [
    double(table_MRF2{1, 9}),  
    double(table_MRF3{1, 9}),  
    double(table_MRF1{1, 9}),  
    double(table_MRF4{1, 9}), 
    double(table_SM{1, 9})
];

% Create final comparison table
T = table(Method, w_omega, h, T_max, Q, Cd_t, Tavg,V_vanes, ...
    'VariableNames', {'Method', 'w_omega', 'h', 'T_max', 'Q','Cd_t', 'Tavg','V_vanes' });

%% Variables for PCA
vars = {'h','T_max','Q','Tavg'};    % variables for PCA
outPrefix = 'PCA_output';           % file prefix for saved results

% Check columns
miss = setdiff(vars, T.Properties.VariableNames);
if ~isempty(miss)
    error('Missing required columns: %s', strjoin(miss, ', '));
end

% Remove rows with NaN in these variables
T = rmmissing(T, 'DataVariables', vars);
n = height(T);
if n < 2
    error('Too few observations after removing NaN: n = %d. PCA not possible.', n);
end

% Group for coloring (Method) if exists
if ismember('Method', T.Properties.VariableNames)
    grp = categorical(T.Method);
else
    grp = categorical(repmat("All", n, 1));
end

% w_omega for optional labeling
hasw_omega = ismember('w_omega', T.Properties.VariableNames);

%% Enable LaTeX font for all figures
set(groot, 'defaultAxesTickLabelInterpreter','latex'); 
set(groot, 'defaultLegendInterpreter','latex');
set(groot, 'defaultTextInterpreter','latex');

%% Build matrix and normalize (z-score)
X = T{:, vars};           
mu = mean(X, 1);
sigma = std(X, [], 1);
sigma(sigma==0) = 1;
Xz = (X - mu) ./ sigma;   

%% Run PCA
[coeff, score, latent, tsquared, explained] = pca(Xz); 

%% Scree / Pareto plot
figure('Units','centimeters','Position',[2 2 18 10]);
subplot(1,2,1);
pareto(explained);
xlabel('Principal Components');
ylabel('Explained Variance (\%)');
title('Pareto - Explained Variance by Component');
subplot(1,2,2);
plot(cumsum(explained), '-o', 'LineWidth',1.4); grid on;
xlabel('Number of Components');
ylabel('Cumulative Explained Variance (\%)');
title('Cumulative Variance');
saveas(gcf, [outPrefix '_scree.png']);

%% PC1 vs PC2 scatter (colored by Method)
pcx = 1; pcy = 2;
figure('Units','centimeters','Position',[2 2 16 12]); hold on; box on;
uCats = categories(grp);
colors = lines(numel(uCats));
markers = {'o','s','^','d','v','>','<','p','h'};

for k = 1:numel(uCats)
    idx = grp == uCats{k};
    scatter(score(idx,pcx), score(idx,pcy), 80, 'Marker', markers{mod(k-1,numel(markers))+1}, ...
        'MarkerEdgeColor','k', 'MarkerFaceColor', colors(k,:), 'DisplayName', char(uCats{k}));
end

% Annotate w_omega if exists
if hasw_omega
    for i = 1:n
        text(score(i,pcx)+0.02, score(i,pcy), sprintf('$\\omega$=%g', T.w_omega(i)), 'FontSize',8);
    end
end

xlabel(sprintf('PC%d (%.1f\\%%)', pcx, explained(pcx)));
ylabel(sprintf('PC%d (%.1f\\%%)', pcy, explained(pcy)));
title(sprintf('PCA: PC%d vs PC%d', pcx, pcy));
legend('Location','best'); grid on;

%% Draw loadings (arrows) over PC1-PC2
scale = max(abs(score(:,[pcx pcy])), [], 'all') * 0.8;
L = coeff(:, [pcx pcy]) * scale;
for j = 1:size(L,1)
    quiver(0, 0, L(j,1), L(j,2), 0, 'k', 'LineWidth', 1.2, 'MaxHeadSize', 0.6);
    text(L(j,1)*1.08, L(j,2)*1.08, vars{j}, 'FontWeight','bold', 'FontSize',10);
end

saveas(gcf, [outPrefix '_PC1_vs_PC2_biplot.png']);

%% Biplot
figure('Units','centimeters','Position',[2 2 14 10]);
biplot(coeff(:,[pcx pcy]), 'Scores', score(:,[pcx pcy]), 'VarLabels', vars);
title('Biplot (PC1 vs PC2)');

%% Export results
pcNames = arrayfun(@(k) sprintf('PC%d', k), 1:size(score,2), 'UniformOutput', false);
scoreTbl = array2table(score, 'VariableNames', pcNames);
outTbl = [T(:, intersect({'Method','w_omega'}, T.Properties.VariableNames)), scoreTbl];
writetable(outTbl, [outPrefix '_scores.csv']);

nPC = size(coeff,2);
pcNamesDyn = arrayfun(@(k) sprintf('PC%d', k), 1:nPC, 'UniformOutput', false);
writetable([table(vars','VariableNames',{'Var'}) array2table(coeff, 'VariableNames', pcNamesDyn)], ...
    [outPrefix '_loadings.csv']);

explTbl = table((1:numel(explained))', explained, cumsum(explained), ...
    'VariableNames', {'PC','ExplainedPct','CumExplainedPct'});
writetable(explTbl, [outPrefix '_explained.csv']);

fprintf('PCA completed. Results saved with prefix: %s_\n', outPrefix);
