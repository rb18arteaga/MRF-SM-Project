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

%% Define data: Method, w_omega, htc
Method = ["MRF"; "MRF"; "MRF"; "MRF"; "SM"];
w_omega = [31.7, 44.1, 69.7, 99.5, 69.7]';
htc = [
    double(table_MRF2{1, 11}),
    double(table_MRF3{1, 11}),
    double(table_MRF1{1, 11}),
    double(table_MRF4{1, 11}),
    double(table_SM{1, 11})
];

data = table(Method, w_omega, htc, 'VariableNames', {'Method', 'w_omega', 'h'});



%% Regression preparation
data.Method = categorical(data.Method);
data = rmmissing(data);
data.w_omega2 = data.w_omega.^2;

%% Separate MRF and SM datasets
mrf = data(data.Method == 'MRF', :);
sm  = data(data.Method == 'SM', :);

if height(mrf) < 2
    error('At least 2 MRF points are required for linear fit (>=3 recommended for quadratic).');
end

%% Fit linear and quadratic models
lm_lin = fitlm(mrf, 'h ~ w_omega');
lm_quad = fitlm(mrf, 'h ~ w_omega + w_omega2');

disp('--- Linear Model (MRF) ---'); disp(lm_lin);
disp('--- Quadratic Model (MRF) ---'); disp(lm_quad);

AIC_lin = lm_lin.ModelCriterion.AIC;
AIC_quad = lm_quad.ModelCriterion.AIC;
fprintf('AIC linear: %.3f, AIC quadratic: %.3f\n', AIC_lin, AIC_quad);

%% Centered models
mrf.w_omega_c = mrf.w_omega - mean(mrf.w_omega);
mrf.w_omega_c2 = mrf.w_omega_c.^2;

lm_lin_c = fitlm(mrf, 'h ~ w_omega_c');
lm_quad_c = fitlm(mrf, 'h ~ w_omega_c + w_omega_c2');

disp('--- Centered Linear Model ---'); disp(lm_lin_c);
disp('--- Centered Quadratic Model ---'); disp(lm_quad_c);

%% Corrected AIC (AICc)
n_mrf = height(mrf);
k_lin = lm_lin.NumCoefficients;
k_quad = lm_quad.NumCoefficients;

AICc_lin = AIC_lin + (2*k_lin*(k_lin+1)) / (n_mrf - k_lin - 1);
AICc_quad = AIC_quad + (2*k_quad*(k_quad+1)) / (n_mrf - k_quad - 1);

fprintf('AICc linear: %.3f, AICc quadratic: %.3f\n', AICc_lin, AICc_quad);

%% Prediction grid
w_omega_min = min(data.w_omega);
w_omega_max = max(data.w_omega);
w_omegaGrid = linspace(w_omega_min, w_omega_max, 200)';
predictTbl = table(w_omegaGrid, w_omegaGrid.^2, 'VariableNames', {'w_omega','w_omega2'});

[y_lin, ylin_CI] = predict(lm_lin, predictTbl);
[y_quad, yquad_CI] = predict(lm_quad, predictTbl);

%% Predict for SM point (if exists)
if ~isempty(sm)
    sm_w_omega = sm.w_omega(1);
    sm_h = sm.h(1);
    smTbl = table(sm_w_omega, sm_w_omega.^2, 'VariableNames', {'w_omega','w_omega2'});
    [h_sm_pred_quad, hsm_ci] = predict(lm_quad, smTbl);
    fprintf('SM: w_omega=%.1f, h_obs=%.3f, h_pred_quad=%.3f, CI=[%.3f, %.3f]\n', ...
            sm_w_omega, sm_h, h_sm_pred_quad, hsm_ci(1), hsm_ci(2));
end

%% Plot 1: Confidence bands
figure('Units','centimeters','Position',[2 2 18 10]); hold on; box on;
set(gca,'TickLabelInterpreter','latex');
set(findall(gca,'Type','text'),'Interpreter','latex');

scatter(mrf.w_omega, mrf.h, 60, 'o', 'filled', 'DisplayName','MRF (data)');
if ~isempty(sm)
    scatter(sm_w_omega, sm_h, 100, 's', 'filled', ...
        'MarkerEdgeColor','k', 'MarkerFaceColor',[0.85 0.3 0.3], ...
        'DisplayName','SM (observed)');
end
plot(w_omegaGrid, y_lin, 'b-', 'LineWidth',1.6, 'DisplayName','Linear model (MRF)');
plot(w_omegaGrid, y_quad, 'g-', 'LineWidth',1.6, 'DisplayName','Quadratic model (MRF)');

fill([w_omegaGrid; flipud(w_omegaGrid)], ...
     [ylin_CI(:,1); flipud(ylin_CI(:,2))], ...
     [0.8 0.9 1], 'EdgeColor','none', 'FaceAlpha', 0.35);
fill([w_omegaGrid; flipud(w_omegaGrid)], ...
     [yquad_CI(:,1); flipud(yquad_CI(:,2))], ...
     [0.85 1 0.85], 'EdgeColor','none', 'FaceAlpha', 0.25);

xlabel('$\omega$', 'Interpreter','latex');
ylabel('$h \; (W/m^2K)$', 'Interpreter','latex');
title('$h$ vs $\omega$: MRF models and SM point', 'Interpreter','latex');
legend('Location','best','Interpreter','latex');
grid on; hold off;

%% Plot 2: Prediction bands
figure('Units','centimeters','Position',[2 2 18 10]); hold on; box on;
set(gca,'TickLabelInterpreter','latex');
set(findall(gca,'Type','text'),'Interpreter','latex');

scatter(mrf.w_omega, mrf.h, 60, 'o', 'filled', 'DisplayName','MRF (data)');
if ~isempty(sm)
    scatter(sm_w_omega, sm_h, 100, 's', 'filled', ...
        'MarkerEdgeColor','k', 'MarkerFaceColor',[0.85 0.3 0.3], ...
        'DisplayName','SM (observed)');
end

[y_lin_pred, ylin_PI] = predict(lm_lin, predictTbl, 'Prediction', 'observation');
[y_quad_pred, yquad_PI] = predict(lm_quad, predictTbl, 'Prediction', 'observation');

plot(w_omegaGrid, y_lin_pred, 'b-', 'LineWidth', 1.6, 'DisplayName','Linear (MRF)');
plot(w_omegaGrid, y_quad_pred, 'g-', 'LineWidth', 1.6, 'DisplayName','Quadratic (MRF)');

fill([w_omegaGrid; flipud(w_omegaGrid)], ...
     [ylin_PI(:,1); flipud(ylin_PI(:,2))], ...
     [0.6 0.8 1], 'EdgeColor','none', 'FaceAlpha', 0.25, ...
     'DisplayName','Prediction (linear)');

fill([w_omegaGrid; flipud(w_omegaGrid)], ...
     [yquad_PI(:,1); flipud(yquad_PI(:,2))], ...
     [0.7 1 0.7], 'EdgeColor','none', 'FaceAlpha', 0.25, ...
     'DisplayName','Prediction (quadratic)');

xlabel('$\omega$', 'Interpreter','latex');
ylabel('$h \; (W/m^2K)$', 'Interpreter','latex');
title('$h$ vs $\omega$: Prediction bands', 'Interpreter','latex');
legend('Location','best','Interpreter','latex');
grid on; hold off;

%% Linear Mixed-Effects Model
try
    lme = fitlme(data, 'h ~ w_omega + (1|Method)');
    disp('--- Linear Mixed-Effects model ---'); disp(lme);
    disp('Random effects estimates:'); disp(randomEffects(lme));
catch ME
    warning('Could not fit LME: %s', ME.message);
    lme = [];
end

%% Fixed-effects models
mdl_fixed = fitlm(data, 'h ~ w_omega + Method');
disp('--- Fixed-effects Model ---'); disp(mdl_fixed);

mdl_inter = fitlm(data, 'h ~ w_omega*Method');
disp('--- Interaction Model (w_omega*Method) ---'); disp(mdl_inter);

%% Bootstrap for Quadratic Model Coefficients (MRF)
nboot = 2000;
n = height(mrf);
if n < 2
    warning('Bootstrap not performed: too few MRF points.');
else
    coefs_boot = zeros(nboot, 3);
    for b = 1:nboot
        idx = randsample(n, n, true);
        mrf_b = mrf(idx, :);
        lm_b = fitlm(mrf_b, 'h ~ w_omega + w_omega2');
        coefs_boot(b, :) = lm_b.Coefficients.Estimate';
    end
    ci_lower = prctile(coefs_boot, 2.5);
    ci_upper = prctile(coefs_boot, 97.5);
    coef_names = lm_quad.Coefficients.Properties.RowNames;
    est = lm_quad.Coefficients.Estimate;
    T = table(coef_names, est, ci_lower', ci_upper', ...
              'VariableNames', {'Coef','Estimate','CI_2p5','CI_97p5'});
    disp('--- Bootstrap CIs (quadratic, MRF) ---'); disp(T);
end

%% Residual diagnostics
figure;
plotResiduals(lm_quad,'fitted');
title('Residuals vs Fitted (Quadratic MRF)', 'Interpreter','latex');
set(gca,'TickLabelInterpreter','latex');

%% Save processed data
writetable(data, 'CFD_results_processed.csv');
%% 
%% 
%% 