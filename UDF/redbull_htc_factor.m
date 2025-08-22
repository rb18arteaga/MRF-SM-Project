% redbull_htc_factor.m
% MATLAB script focused only on computing HTC correction factor f to adjust MRF results to SM.
% Minimal statistics — intended for the case: 4 MRF runs + 1 SM run.
%
% Expected CSV: results_summary.csv
% Required columns (case-insensitive): Method, RPM_rad_s, RPM_rpm, h_mean_Wm2K, T_max_K, Q_W, Cd, case_id
%
% Output: results_factor.csv in ./results/ with original and adjusted h (and Q if available)
%
% Author: Daniel Arteaga
% Date: 2025
clear; close all; clc;


% ---------- SETTINGS ----------
interp_method = 'pchip';    % interpolation method if needed ('linear','pchip','spline')
allow_interp = true;        % allow interpolation if SM omega not present in MRF omegas
outdir = 'results';
if ~exist(outdir,'dir'), mkdir(outdir); end
verbose = true;

% Regression settings (only linear added)
use_regression = true;
regression_R2_threshold = 0.8; % if R2 >= threshold, use linear prediction; otherwise fallback to interp
% ------------------------------

%% Read data
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
        warning('Variable no clasificada: %s', name_actual);
    end
end


table_MRF1 = struct2table(MRF1);
table_MRF2 = struct2table(MRF2);
table_MRF3 = struct2table(MRF3);
table_MRF4 = struct2table(MRF4);
table_SM = struct2table(SM);


disp('Datos de MRF1:');
disp(table_MRF1);


%% 
% %speed MRF2, MRF3, MRF1, MRF4, SM

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


data = table(Method, w_omega, h, T_max, Q, Cd_t, Tavg, 'VariableNames', {'Method', 'w_omega', 'h', 'T_max', 'Q','Cd_t', 'Tavg' });

% 1) Load input table 'data' or try CSVs
if ~exist('data','var')
    % try common filenames
    if exist('results_summary.csv','file')
        data = readtable('results_summary.csv');
    elseif exist('average_results.csv','file')
        data = readtable('average_results.csv');
    else
        error('No variable "data" found and no results_summary.csv or average_results.csv file available.');
    end
end

% ensure column names consistent
% Expected columns: Method, w_omega, h, T_max, Q, Cd_t, Tavg
vars = data.Properties.VariableNames;
lowervars = lower(vars);

% find required columns robustly
findcol = @(names) find(ismember(lowervars, lower(names)),1);
ci_method = findcol({'Method','method'});
ci_omega  = findcol({'w_omega','omega','rpm_rad_s','rpm_rad_s'});
ci_h      = findcol({'h','h_mean_wm2k','h_mean_wm2k'});
ci_Q      = findcol({'q','q_w','q_w'});
ci_tmax   = findcol({'t_max','t_max_k','T_max_K'});

if isempty(ci_method) || isempty(ci_omega) || isempty(ci_h)
    error('Input table must contain columns: Method, w_omega (omega), h (h mean). Found: %s', strjoin(vars,','));
end

% normalize table columns
T = table();
T.Method  = string(data{:,ci_method});
T.w_omega = double(data{:,ci_omega});
T.h       = double(data{:,ci_h});
if ~isempty(ci_Q), T.Q = double(data{:,ci_Q}); end
if ~isempty(ci_tmax), T.T_max = double(data{:,ci_tmax}); end

% ensure Method strings trimmed and uppercase/lower for detection
mstr = lower(strtrim(string(T.Method)));
isSM = (mstr == "sm") | contains(mstr,'sliding');
isMRF = (mstr == "mrf") | contains(mstr,'mrf');

if ~any(isSM)
    % try exact "SM" case-insensitive
    isSM = strcmpi(mstr,'SM');
end
if ~any(isMRF)
    isMRF = strcmpi(mstr,'MRF');
end

if ~any(isMRF)
    error('No MRF rows detected in your input table (Method column).');
end
if ~any(isSM)
    error('No SM row detected in your input table (Method column). Need at least one SM run.');
end

T_MRF = T(isMRF,:);
T_SM  = T(isSM,:);

% If multiple SM rows, use first (warn)
if height(T_SM) > 1
    warning('Multiple SM rows detected. The script will use the FIRST SM row as anchor.');
end
SM_row = T_SM(1,:);

% Extract arrays
omega_MRF = T_MRF.w_omega;
h_MRF     = T_MRF.h;

omega_SM = SM_row.w_omega;
h_SM     = SM_row.h;

if verbose
    fprintf('Detected %d MRF rows and %d SM row(s). Using SM at omega = %.4f rad/s, h_SM = %.6f\n', ...
        height(T_MRF), height(T_SM), omega_SM, h_SM);
end

% 2) obtain h_MRF at SM omega (exact match or interpolation or linear regression)
tol = 1e-9;
idx_exact = find(abs(omega_MRF - omega_SM) <= tol, 1);
matched_note = '';

if ~isempty(idx_exact)
    h_MRF_at_SM = h_MRF(idx_exact);
    matched_note = 'exact';
else
    % No exact match: try linear regression first (if enabled)
    used_regression = false;
    if use_regression && numel(omega_MRF) >= 2
        % linear fit
        p = polyfit(omega_MRF, h_MRF, 1); % p(1)*omega + p(2)
        h_fit = polyval(p, omega_MRF);
        % compute R^2
        SS_res = sum((h_MRF - h_fit).^2);
        SS_tot = sum((h_MRF - mean(h_MRF)).^2);
        if SS_tot == 0
            R2 = 1.0; % all points equal
        else
            R2 = 1 - SS_res / SS_tot;
        end
        % predict at SM omega
        h_pred_lin = polyval(p, omega_SM);
        if verbose
            fprintf('Linear regression: h = a + b*omega -> a=%g, b=%g, R2=%g\n', p(2), p(1), R2);
            fprintf('Linear prediction at omega_SM=%.4f -> h_pred_lin=%.6f\n', omega_SM, h_pred_lin);
        end
        % decide whether to use regression
        if R2 >= regression_R2_threshold
            h_MRF_at_SM = h_pred_lin;
            matched_note = sprintf('linear_regression (R2=%.3f)', R2);
            used_regression = true;
        else
            if verbose
                fprintf('Linear R2 = %.3f < threshold %.3f -> will fallback to interpolation (if allowed)\n', R2, regression_R2_threshold);
            end
        end
    end

    % If regression not used or disabled, use interpolation (if allowed)
    if ~used_regression
        if allow_interp
            [omega_sorted, sidx] = sort(omega_MRF);
            h_sorted = h_MRF(sidx);
            if numel(omega_sorted) < 2
                error('Not enough MRF points for interpolation.');
            end
            h_MRF_at_SM = interp1(omega_sorted, h_sorted, omega_SM, interp_method, 'extrap');
            matched_note = sprintf('interpolated (%s)', interp_method);
            if verbose
                fprintf('Interpolation at omega_SM=%.4f -> h_interp=%.6f\n', omega_SM, h_MRF_at_SM);
            end
        else
            error('No exact MRF at SM omega and interpolation disabled.');
        end
    end
end

% 3) compute correction factor f_ref (h-based)
f_ref = h_SM / h_MRF_at_SM;

% 4) apply factor globally to all MRF h values
h_MRF_adj = h_MRF * f_ref;

% 5) build output table and save
outTab = table();
outTab.Method = T_MRF.Method;
outTab.w_omega = T_MRF.w_omega;
outTab.h_MRF = h_MRF;
outTab.h_MRF_adj = h_MRF_adj;
outTab.f_applied = repmat(f_ref, height(T_MRF), 1);
outTab.isMatchedToSM = false(height(T_MRF),1);
if ~isempty(idx_exact)
    outTab.isMatchedToSM(idx_exact) = true;
end

% Add SM reference column only on matched row (for readability)
outTab.h_SM_reference = NaN(height(outTab),1);
if ~isempty(idx_exact)
    outTab.h_SM_reference(idx_exact) = h_SM;
end

writetable(outTab, fullfile(outdir,'results_factor.csv'));
if verbose, fprintf('Wrote results to %s\n', fullfile(outdir,'results_factor.csv')); end

% 6) summary & plot
fprintf('--- SUMMARY ---\n');
fprintf('h_SM (omega=%.4f) = %.6f\n', omega_SM, h_SM);
fprintf('h_MRF used (at SM omega) = %.6f (%s)\n', h_MRF_at_SM, matched_note);
fprintf('Global correction factor f_ref = %.6f\n', f_ref);
fprintf('Applied to %d MRF cases. Adjusted h range: [%.6f, %.6f]\n', numel(h_MRF_adj), min(h_MRF_adj), max(h_MRF_adj));

% plot
figure('Visible','off');
plot(T_MRF.w_omega, T_MRF.h, 'o-', 'LineWidth',1.5); hold on;
plot(T_MRF.w_omega, h_MRF_adj, 's--', 'LineWidth',1.5);
plot(omega_SM, h_SM, 'kp', 'MarkerFaceColor','y', 'MarkerSize',12);
xlabel('\omega (rad/s)');
ylabel('Area-weighted mean h (W/m^2K)');
legend('h_{MRF}','h_{MRF} adjusted (global f)','h_{SM}','Location','best');
grid on;
saveas(gcf, fullfile(outdir,'h_before_after.png'));
if verbose, fprintf('Saved figure: %s\n', fullfile(outdir,'h_before_after.png')); end

% Save f_ref to a small text file for UDF reading
fid = fopen(fullfile(outdir,'f_ref.txt'),'w');
fprintf(fid, 'f_ref %.12g\n', f_ref);
fprintf(fid, 'omega_SM %.12g\n', omega_SM);
fclose(fid);
if verbose, fprintf('Saved f_ref to %s\n', fullfile(outdir,'f_ref.txt')); end

% Done