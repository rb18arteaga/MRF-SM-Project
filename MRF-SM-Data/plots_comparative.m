clear all 
clc

%% Read average
dataREAD = readtable('average_results.csv');

name = dataREAD.OriginalVariableNames;  
values = dataREAD.Var1;  

MRF1 = struct();
MRF2 = struct();
MRF3 = struct();
MRF4 = struct();
SM = struct();

for i = 1:length(name)
    name_actual = name{i};
    value_actual = values(i);
    
    if startsWith(name_actual, 'MRF1_')
        field = strrep(name_actual, 'MRF1_', ''); 
        MRF1.(field) = value_actual;
    elseif startsWith(name_actual, 'MRF2_')
        field = strrep(name_actual, 'MRF2_', '');
        MRF2.(field) = value_actual;
    elseif startsWith(name_actual, 'MRF3_')
        field = strrep(name_actual, 'MRF3_', '');
        MRF3.(field) = value_actual;
    elseif startsWith(name_actual, 'MRF4_')
        field = strrep(name_actual, 'MRF4_', '');
        MRF4.(field) = value_actual;
    elseif startsWith(name_actual, 'SM_')
        field = strrep(name_actual, 'SM_', '');
        SM.(field) = value_actual; 
    else
        warning('Unclassified variable %s', name_actual);
    end
end

table_MRF1 = struct2table(MRF1);
table_MRF2 = struct2table(MRF2);
table_MRF3 = struct2table(MRF3);
table_MRF4 = struct2table(MRF4);
table_SM = struct2table(SM);

%% Main table
Method = ["MRF2"; "MRF3"; "MRF1"; "MRF4";"SM"];
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

T = table(Method, w_omega, h, T_max, Q, Cd_t, Tavg, ...
    'VariableNames', {'Method', 'w_omega', 'h', 'T_max', 'Q','Cd_t', 'Tavg' });

%% ==== Comparison A ====
T_max_exp1 = 500; % K
T_avg_mrf2 = double(table_MRF2{1, 12}); % K https://www.sciencedirect.com/science/article/pii/S1290072917316149?via%3Dihub#bib2
err_rel_A = abs(T_avg_mrf2 - T_max_exp1) / T_max_exp1 * 100;
T_max_sim_min = T_avg_mrf2 * (1 - err_rel_A/100);
T_max_sim_max = T_avg_mrf2 * (1 + err_rel_A/100);

figure('Color','w','Position',[100 100 600 400]);
b = bar([T_max_exp1, T_avg_mrf2], 'FaceColor','flat');
b.CData(1,:) = [0, 0.45, 0.74];
b.CData(2,:) = [0.85, 0.33, 0.1];
set(gca,'XTickLabel',{'Experimental','MRF2 (Simulation)'}, ...
    'TickLabelInterpreter','latex','FontSize',12);
ylabel('$T_{\mathrm{avg}}$ [K]', 'Interpreter','latex', 'FontSize',14);
title('Comparison A: $T_{\mathrm{avg}}$', 'Interpreter','latex', 'FontSize',14);
grid on; hold on;

errorbar(2, T_avg_mrf2, T_avg_mrf2 - T_max_sim_min, ...
    T_max_sim_max - T_avg_mrf2, 'k','LineWidth',1.5,'CapSize',12);

text(1, T_max_exp1 + 8, sprintf('%.1f', T_max_exp1), ...
    'HorizontalAlignment','center','FontSize',12,'Interpreter','latex');
text(2.2, T_avg_mrf2, sprintf('Err=%.1f\\%%', err_rel_A), ...
    'HorizontalAlignment','left','FontSize',12, 'VerticalAlignment','top','Interpreter','latex');
text(2, T_avg_mrf2 + 8, sprintf('%.1f', T_avg_mrf2), ...
    'HorizontalAlignment','left','FontSize',12,'VerticalAlignment','middle','Interpreter','latex');

%% ==== Comparison A2 ====

% Experimental & simulation data
htc_exp1 = 33; % W/m2K https://www.sciencedirect.com/science/article/pii/S1290072917316149?via%3Dihub#bib2
htc_mrf2 = double(table_MRF2{1, 11}); % W/m2K 

% Relative error (%)
err_rel_A2 = abs(htc_mrf2 - htc_exp1) / htc_exp1 * 100;

% Error bars: ± relative error
htc_sim_min = htc_mrf2 * (1 - err_rel_A2/100);
htc_sim_max = htc_mrf2 * (1 + err_rel_A2/100);

% Create figure
figure('Color','w','Position',[100 100 600 400]);
b = bar([htc_exp1, htc_mrf2], 'FaceColor','flat');
b.CData(1,:) = [0, 0.45, 0.74]; % Experimental color
b.CData(2,:) = [0.85, 0.33, 0.1]; % Simulation color

% Axis formatting
set(gca,'XTickLabel',{'Experimental','MRF2 (Simulation)'}, ...
    'TickLabelInterpreter','latex','FontSize',12);
ylabel('$HTC$ [W/m$^2$K]', 'Interpreter','latex', 'FontSize',14);
title('Comparison A2: $HTC$', 'Interpreter','latex', 'FontSize',14);
grid on; hold on;

% Error bars
errorbar(2, htc_mrf2, htc_mrf2 - htc_sim_min, ...
    htc_sim_max - htc_mrf2, 'k','LineWidth',1.5,'CapSize',12);

% Text annotations
text(1, htc_exp1 + 1, sprintf('%.1f', htc_exp1), ...
    'HorizontalAlignment','center','FontSize',12,'Interpreter','latex');
text(2.2, htc_mrf2, sprintf('Err=%.1f\\%%', err_rel_A2), ...
    'HorizontalAlignment','left','FontSize',12, ...
    'VerticalAlignment','top','Interpreter','latex');
text(2, htc_mrf2 + 1, sprintf('%.1f', htc_mrf2), ...
    'HorizontalAlignment','left','FontSize',12, ...
    'VerticalAlignment','middle','Interpreter','latex');




%% ==== Comparison B ====
data_htcRpm_exp = readmatrix("2htc-rpm.csv"); %https://publications.lib.chalmers.se/records/fulltext/255591/255591.pdf
rad_s_exp2 = data_htcRpm_exp(:,1)*((2*pi)/60);  
htcW_exp2 = data_htcRpm_exp(:,2); 

sim_radps = T{:,2};
sim_HTC = T{:,3};

colores_mrf = lines(4);
color_sm = [0.85, 0.33, 0.1];

figure('Color','w','Position',[100 100 600 400]);
scatter(rad_s_exp2, htcW_exp2, 70, '^', ...
    'MarkerFaceColor',[0.5,0.5,0.5], 'MarkerEdgeColor','k', ...
    'DisplayName','Experimental');
hold on;

for i = 1:height(T)
    if startsWith(T.Method{i}, 'MRF')
        c = colores_mrf(strcmp(T.Method{i}, {'MRF2','MRF3','MRF1','MRF4'}), :);
        scatter(sim_radps(i), sim_HTC(i), 70, 'o', ...
            'MarkerFaceColor',c, 'MarkerEdgeColor','k', ...
            'DisplayName', sprintf('%s ($\\omega$=%.1f)', T.Method{i}, T.w_omega(i)));
        text(sim_radps(i), sim_HTC(i)+2, T.Method{i}, ...
            'FontSize',10, 'HorizontalAlignment','center','Interpreter','latex');
    else
        scatter(sim_radps(i), sim_HTC(i), 70, 's', ...
            'MarkerFaceColor',color_sm, 'MarkerEdgeColor','k', ...
            'DisplayName', sprintf('%s ($\\omega$=%.1f)', T.Method{i}, T.w_omega(i)));
        text(sim_radps(i), sim_HTC(i)+2, T.Method{i}, ...
            'FontSize',10, 'HorizontalAlignment','center','Interpreter','latex');
    end
end

idx_mrf = startsWith(T.Method,'MRF');
p_mrf = polyfit(sim_radps(idx_mrf), sim_HTC(idx_mrf), 1);
x_fit = linspace(min(sim_radps(idx_mrf)), max(sim_radps(idx_mrf)), 100);
y_fit = polyval(p_mrf, x_fit);
plot(x_fit, y_fit, '--', 'Color',[0 0 0], 'LineWidth', 1.5, ...
    'DisplayName','Linear fit (MRF)');

xlabel('$\omega$ [rad/s]', 'Interpreter','latex','FontSize',14);
ylabel('HTC [W/m$^2$K]', 'Interpreter','latex','FontSize',14);
title('Comparison B: HTC vs $\omega$', 'Interpreter','latex','FontSize',14);
legend('Location','best','Interpreter','latex');
grid on;

%% ==== Comparison C ====
data_htcTemp_exp = readmatrix("htc-temp-exp-chalmers.csv"); %https://publications.lib.chalmers.se/records/fulltext/255591/255591.pdf
temp_K_exp2 = data_htcTemp_exp(:,1) + 273; 
htcT_exp2 = data_htcTemp_exp(:,2);

sim_Temp = T{:,4}; % K

figure('Color','w','Position',[100 100 600 400]);
scatter(temp_K_exp2, htcT_exp2, 70, '^', ...
    'MarkerFaceColor',[0.5,0.5,0.5], 'MarkerEdgeColor','k', ...
    'DisplayName','Experimental');
hold on;

for i = 1:height(T)
    if startsWith(T.Method{i}, 'MRF')
        c = colores_mrf(strcmp(T.Method{i}, {'MRF2','MRF3','MRF1','MRF4'}), :);
        scatter(sim_Temp(i), sim_HTC(i), 70, 'o', ...
            'MarkerFaceColor',c, 'MarkerEdgeColor','k', ...
            'DisplayName', sprintf('%s (T=%.1f)', T.Method{i}, T.T_max(i)));
        text(sim_Temp(i), sim_HTC(i)+2, T.Method{i}, ...
            'FontSize',10, 'HorizontalAlignment','center','Interpreter','latex');
    else
        scatter(sim_Temp(i), sim_HTC(i), 70, 's', ...
            'MarkerFaceColor',color_sm, 'MarkerEdgeColor','k', ...
            'DisplayName', sprintf('%s (T=%.1f)', T.Method{i}, T.T_max(i)));
        text(sim_Temp(i), sim_HTC(i)+2, T.Method{i}, ...
            'FontSize',10, 'HorizontalAlignment','center','Interpreter','latex');
    end
end

idx_mrf = startsWith(T.Method,'MRF');
p_mrf_T = polyfit(sim_Temp(idx_mrf), sim_HTC(idx_mrf), 1);
x_fit_T = linspace(min(sim_Temp(idx_mrf)), max(sim_Temp(idx_mrf)), 100);
y_fit_T = polyval(p_mrf_T, x_fit_T);
plot(x_fit_T, y_fit_T, '--', 'Color',[0 0 0], 'LineWidth', 1.5, ...
    'DisplayName','Linear fit (MRF)');

xlabel('Temperature [K]', 'Interpreter','latex','FontSize',14);
ylabel('HTC [W/m$^2$K]', 'Interpreter','latex','FontSize',14);
title('Comparison C: HTC vs Temperature', 'Interpreter','latex','FontSize',14);
legend('Location','best','Interpreter','latex');
%grid off;