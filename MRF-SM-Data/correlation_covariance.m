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
data = table(Method, w_omega, h, T_max, Q, Cd_t, Tavg,V_vanes, ...
    'VariableNames', {'Method', 'w_omega', 'h', 'T_max', 'Q','Cd_t', 'Tavg','V_vanes' });

vars = {'h', 'T_max', 'Q','Cd_t', 'Tavg'}; % variables to analyze
subsetOption = 'all'; % 'all' | 'MRF' | 'SM' | 'paired'
nBootstrap = 2000;     % bootstrap for Spearman CI

% Convert Method if it exists
if ismember('Method', data.Properties.VariableNames)
    data.Method = categorical(data.Method);
end

%% Select subset
switch lower(subsetOption)
    case 'all'
        dataSub = data;
    case 'mrf'
        if ~ismember('Method', data.Properties.VariableNames), error('No Method column to filter MRF'); end
        dataSub = data(data.Method=='MRF', :);
    case 'sm'
        if ~ismember('Method', data.Properties.VariableNames), error('No Method column to filter SM'); end
        dataSub = data(data.Method=='SM', :);
    case 'paired'
        if ~ismember('Method', data.Properties.VariableNames), error('No Method column to pair'); end
        w_omega_mrf = unique(data.w_omega(data.Method=='MRF'));
        w_omega_sm  = unique(data.w_omega(data.Method=='SM'));
        commonw_omega = intersect(w_omega_mrf, w_omega_sm);
        dataSub = data(ismember(data.w_omega, commonw_omega), :);
    otherwise
        error('Unrecognized subsetOption');
end

if isempty(dataSub)
    error('No data in selected subset.');
end

%% Prepare matrices and counters
nTotal = height(dataSub);
fprintf('Using %d rows for analysis (subset=%s)\n', nTotal, subsetOption);

w_omega = dataSub.w_omega;
X = zeros(nTotal, numel(vars));
for i=1:numel(vars)
    X(:,i) = dataSub{:, vars{i}};
end

%% Correlations (Pearson and Spearman) with w_omega
results = table('Size',[numel(vars), 10], ...
    'VariableTypes',{'string','double','double','double','double','double','double','double','double','double'}, ...
    'VariableNames',{'Var','Pearson_r','Pearson_p','Pearson_CI_low','Pearson_CI_high', ...
                     'Spearman_rho','Spearman_p','Spearman_CI_low','Spearman_CI_high','Spearman_CI_width'});
results.Var = vars';

for i=1:numel(vars)
    y = X(:,i);
    ok = ~isnan(w_omega) & ~isnan(y);
    x_ok = w_omega(ok); y_ok = y(ok);
    n = numel(x_ok);

    if n < 3
        warning('Too few data points for variable %s (n=%d). Skipping.', vars{i}, n);
        results.Pearson_r(i) = NaN; results.Pearson_p(i)=NaN;
        results.Pearson_CI_low(i) = NaN; results.Pearson_CI_high(i) = NaN;
        results.Spearman_rho(i)=NaN; results.Spearman_p(i)=NaN; results.Spearman_CI_width(i)=NaN;
        continue;
    end

    % Pearson
    [r_p, p_p] = corr(x_ok, y_ok, 'Type', 'Pearson', 'Rows', 'complete');
    results.Pearson_r(i) = r_p;
    results.Pearson_p(i) = p_p;
    z = atanh(r_p);
    se = 1/sqrt(n-3);
    zcrit = norminv(0.975);
    r_lo = tanh(z - zcrit*se); 
    r_hi = tanh(z + zcrit*se);
    results.Pearson_CI_low(i) = r_lo;
    results.Pearson_CI_high(i) = r_hi;

    % Spearman
    [rho_s, p_s] = corr(x_ok, y_ok, 'Type', 'Spearman', 'Rows', 'complete');
    results.Spearman_rho(i) = rho_s;
    results.Spearman_p(i) = p_s;

    bootStats = zeros(nBootstrap,1);
    for b = 1:nBootstrap
        idx = randsample(n, n, true);
        xb = x_ok(idx); yb = y_ok(idx);
        bootStats(b) = corr(xb, yb, 'Type', 'Spearman', 'Rows', 'complete');
    end
    ci_s = prctile(bootStats, [2.5 97.5]);
    results.Spearman_CI_low(i) = ci_s(1);
    results.Spearman_CI_high(i) = ci_s(2);
    results.Spearman_CI_width(i) = ci_s(2) - ci_s(1);
end

%% Covariance and correlation matrices
cov_matrix = cov(X, 'omitrows');
[Rpearson, Ppearson] = corr(X, 'Type', 'Pearson', 'Rows', 'complete');
[Rspearman, Pspearman] = corr(X, 'Type', 'Spearman', 'Rows', 'complete');

%% Display results
disp('--- Summary correlations vs w_omega ---');
disp(results);

disp('--- Covariance matrix (variables in this order) ---');
disp(vars);
disp(cov_matrix);

disp('--- Correlation matrix (Pearson) ---');
disp(Rpearson);
disp('p-values (Pearson):'); disp(Ppearson);

%% Visualizations (LaTeX font)
set(groot, 'defaultAxesTickLabelInterpreter','latex'); 
set(groot, 'defaultLegendInterpreter','latex');
set(groot, 'defaultTextInterpreter','latex');

% Heatmap Pearson
figure('Name','Pearson correlation heatmap','Units','centimeters','Position',[2 2 12 10]);
imagesc(Rpearson); axis square; colorbar; colormap(parula);
set(gca,'XTick',1:numel(vars),'XTickLabel',vars,'YTick',1:numel(vars),'YTickLabel',vars);
title('Pearson correlation (variables)', 'Interpreter', 'latex');

for i=1:numel(vars)
    for j=1:numel(vars)
        text(j,i,sprintf('%.2f',Rpearson(i,j)),'HorizontalAlignment','center','Color','w','FontSize',10, 'Interpreter', 'latex');
    end
end

% Scatter plots
for i=1:numel(vars)
    y = X(:,i);
    ok = ~isnan(w_omega) & ~isnan(y);
    x_ok = w_omega(ok); y_ok = y(ok);
    if numel(x_ok) < 2, continue; end
    figure('Name',['w\_omega vs ', vars{i}],'Units','centimeters','Position',[2 2 14 8]);
    scatter(x_ok, y_ok, 50, 'filled'); hold on;
    p = polyfit(x_ok, y_ok, 1);
    xg = linspace(min(x_ok), max(x_ok), 200);
    plot(xg, polyval(p,xg), 'r-', 'LineWidth',1.5);
    r_p = results.Pearson_r(i); p_p = results.Pearson_p(i);
    r_s = results.Spearman_rho(i); p_s = results.Spearman_p(i);
    txt = sprintf('Pearson $r=%.3f$ (p=%.3g)\\\\Spearman $\\rho=%.3f$ (p=%.3g)', r_p, p_p, r_s, p_s);
    xlabel('$\omega$', 'Interpreter','latex'); 
    ylabel(vars{i}, 'Interpreter','latex');
    title(['$\omega$ vs ', vars{i}], 'Interpreter','latex');
    text(0.02, 0.95, txt, 'Units','normalized','VerticalAlignment','top','BackgroundColor','w','EdgeColor','k', 'Interpreter','latex');
    grid on; hold off;
end

%% Save results
writetable(results, 'correlation_summary_vs_w_omega.csv');
writematrix(cov_matrix, 'covariance_matrix.csv');
writematrix(Rpearson, 'correlation_matrix_pearson.csv');

fprintf('Results saved: correlation_summary_vs_w_omega.csv, covariance_matrix.csv, correlation_matrix_pearson.csv\n');