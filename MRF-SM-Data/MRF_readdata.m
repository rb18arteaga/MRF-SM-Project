% file: process_MRF_SM_averages.m
clc
clear all
clearvars

%% ---------------------------
% Read data MRF (steady-state) and SM (transient)
% ----------------------------
% MRF1 - 69.7 rad/s
MRF1_CD_ty =  readmatrix("MRF1-69.7\report-coef-drag-tyre-rfile.out", FileType="text");
MRF1_CD_tyrim =  readmatrix("MRF1-69.7\report-coef-drag-tyre-rim-rfile.out", FileType="text");
MRF1_CL_ty =  readmatrix("MRF1-69.7\report-coef-lift-tyre-rfile.out", FileType="text");
MRF1_CL_tyrim =  readmatrix("MRF1-69.7\report-coef-lift-tyre-rim-rfile.out", FileType="text");

MRF1_maxT_soliddisc =  readmatrix("MRF1-69.7\report-exp-max-sttemp-discsolid-rfile.out", FileType="text");

MRF1_heatflux_disczones =  readmatrix("MRF1-69.7\report-sum-heatflux-disczones-rfile.out", FileType="text");
MRF1_heatflux_soliddisc =  readmatrix("MRF1-69.7\report-sum-heatflux-soliddisc-rfile.out", FileType="text");

MRF1_savgT_soliddisc =  readmatrix("MRF1-69.7\report-surf-arweavg-sttemp-soliddisc-rfile.out", FileType="text");

MRF1_savgV_vanes =  readmatrix("MRF1-69.7\report-surf-arweavg-vanes-velocity-rfile.out", FileType="text");

MRF1_savg_HTC_soliddisc =  readmatrix("MRF1-69.7\report-surf-arweavg-wallflux-shtc-disc-solid-rfile.out", FileType="text");
MRF1_savg_HTC_disczones = readmatrix("MRF1-69.7\report-surfarweavg-wallflux-shtc-mid-disc-rfile.out", FileType="text");

MRF1_vavgT_soliddisc =  readmatrix("MRF1-69.7\report-volavg-statemp-soliddisc-rfile.out", FileType="text");
MRF1_vmaxT_pads12 =  readmatrix("MRF1-69.7\report-volmax-statemp-pads12-rfile.out", FileType="text");

%% MRF2 - 31.7 rad/s
MRF2_CD_ty =  readmatrix("MRF2-31.7\report-coef-drag-tyre-rfile.out", FileType="text");
MRF2_CD_tyrim =  readmatrix("MRF2-31.7\report-coef-drag-tyre-rim-rfile.out", FileType="text");
MRF2_CL_ty =  readmatrix("MRF2-31.7\report-coef-lift-tyre-rfile.out", FileType="text");
MRF2_CL_tyrim =  readmatrix("MRF2-31.7\report-coef-lift-tyre-rim-rfile.out", FileType="text");

MRF2_maxT_soliddisc =  readmatrix("MRF2-31.7\report-exp-max-sttemp-discsolid-rfile.out", FileType="text");

MRF2_heatflux_disczones =  readmatrix("MRF2-31.7\report-sum-heatflux-disczones-rfile.out", FileType="text");
MRF2_heatflux_soliddisc =  readmatrix("MRF2-31.7\report-sum-heatflux-soliddisc-rfile.out", FileType="text");

MRF2_savgT_soliddisc =  readmatrix("MRF2-31.7\report-surf-arweavg-sttemp-soliddisc-rfile.out", FileType="text");

MRF2_savgV_vanes =  readmatrix("MRF2-31.7\report-surf-arweavg-vanes-velocity-rfile.out", FileType="text");

MRF2_savg_HTC_soliddisc =  readmatrix("MRF2-31.7\report-surf-arweavg-wallflux-shtc-disc-solid-rfile.out", FileType="text");
MRF2_savg_HTC_disczones = readmatrix("MRF2-31.7\report-surfarweavg-wallflux-shtc-mid-disc-rfile.out", FileType="text");

MRF2_vavgT_soliddisc =  readmatrix("MRF2-31.7\report-volavg-statemp-soliddisc-rfile.out", FileType="text");
MRF2_vmaxT_pads12 =  readmatrix("MRF2-31.7\report-volmax-statemp-pads12-rfile.out", FileType="text");

%% MRF3 - 44.1 rad/s
MRF3_CD_ty =  readmatrix("MRF3-44.1\report-coef-drag-tyre-rfile.out", FileType="text");
MRF3_CD_tyrim =  readmatrix("MRF3-44.1\report-coef-drag-tyre-rim-rfile.out", FileType="text");
MRF3_CL_ty =  readmatrix("MRF3-44.1\report-coef-lift-tyre-rfile.out", FileType="text");
MRF3_CL_tyrim =  readmatrix("MRF3-44.1\report-coef-lift-tyre-rim-rfile.out", FileType="text");

MRF3_maxT_soliddisc =  readmatrix("MRF3-44.1\report-exp-max-sttemp-discsolid-rfile.out", FileType="text");

MRF3_heatflux_disczones =  readmatrix("MRF3-44.1\report-sum-heatflux-disczones-rfile.out", FileType="text");
MRF3_heatflux_soliddisc =  readmatrix("MRF3-44.1\report-sum-heatflux-soliddisc-rfile.out", FileType="text");

MRF3_savgT_soliddisc =  readmatrix("MRF3-44.1\report-surf-arweavg-sttemp-soliddisc-rfile.out", FileType="text");

MRF3_savgV_vanes =  readmatrix("MRF3-44.1\report-surf-arweavg-vanes-velocity-rfile.out", FileType="text");

MRF3_savg_HTC_soliddisc =  readmatrix("MRF3-44.1\report-surf-arweavg-wallflux-shtc-disc-solid-rfile.out", FileType="text");
MRF3_savg_HTC_disczones = readmatrix("MRF3-44.1\report-surfarweavg-wallflux-shtc-mid-disc-rfile.out", FileType="text");

MRF3_vavgT_soliddisc =  readmatrix("MRF3-44.1\report-volavg-statemp-soliddisc-rfile.out", FileType="text");
MRF3_vmaxT_pads12 =  readmatrix("MRF3-44.1\report-volmax-statemp-pads12-rfile.out", FileType="text");

%% MRF4 - 95.1 rad/s
MRF4_CD_ty =  readmatrix("MRF4-95.1\report-coef-drag-tyre-rfile.out", FileType="text");
MRF4_CD_tyrim =  readmatrix("MRF4-95.1\report-coef-drag-tyre-rim-rfile.out", FileType="text");
MRF4_CL_ty =  readmatrix("MRF4-95.1\report-coef-lift-tyre-rfile.out", FileType="text");
MRF4_CL_tyrim =  readmatrix("MRF4-95.1\report-coef-lift-tyre-rim-rfile.out", FileType="text");

MRF4_maxT_soliddisc =  readmatrix("MRF4-95.1\report-exp-max-sttemp-discsolid-rfile.out", FileType="text");

MRF4_heatflux_disczones =  readmatrix("MRF4-95.1\report-sum-heatflux-disczones-rfile.out", FileType="text");
MRF4_heatflux_soliddisc =  readmatrix("MRF4-95.1\report-sum-heatflux-soliddisc-rfile.out", FileType="text");

MRF4_savgT_soliddisc =  readmatrix("MRF4-95.1\report-surf-arweavg-sttemp-soliddisc-rfile.out", FileType="text");

MRF4_savgV_vanes =  readmatrix("MRF4-95.1\report-surf-arweavg-vanes-velocity-rfile.out", FileType="text");

MRF4_savg_HTC_soliddisc =  readmatrix("MRF4-95.1\report-surf-arweavg-wallflux-shtc-disc-solid-rfile.out", FileType="text");
MRF4_savg_HTC_disczones = readmatrix("MRF4-95.1\report-surfarweavg-wallflux-shtc-mid-disc-rfile.out", FileType="text");

MRF4_vavgT_soliddisc =  readmatrix("MRF4-95.1\report-volavg-statemp-soliddisc-rfile.out", FileType="text");
MRF4_vmaxT_pads12 =  readmatrix("MRF4-95.1\report-volmax-statemp-pads12-rfile.out", FileType="text");

%% SM - Sliding Mesh (transient) - 75 rad/s
SM_CD_ty =  readmatrix("SM\report-coef-drag-tyre-rfile.out", FileType="text");
SM_CD_tyrim =  readmatrix("SM\report-coef-drag-tyre-rim-rfile.out", FileType="text");
SM_CL_ty =  readmatrix("SM\report-coef-lift-tyre-rfile.out", FileType="text");
SM_CL_tyrim =  readmatrix("SM\report-coef-lift-tyre-rim-rfile.out", FileType="text");

SM_maxT_soliddisc =  readmatrix("SM\report-exp-max-sttemp-discsolid-rfile.out", FileType="text");

SM_heatflux_disczones =  readmatrix("SM\report-sum-heatflux-disczones-rfile.out", FileType="text");
SM_heatflux_soliddisc =  readmatrix("SM\report-sum-heatflux-soliddisc-rfile.out", FileType="text");

SM_savgT_soliddisc =  readmatrix("SM\report-surf-arweavg-sttemp-soliddisc-rfile.out", FileType="text");

SM_savgV_vanes =  readmatrix("SM\report-surf-arweavg-vanes-velocity-rfile.out", FileType="text");

SM_savg_HTC_soliddisc =  readmatrix("SM\report-surf-arweavg-wallflux-shtc-disc-solid-rfile.out", FileType="text");
SM_savg_HTC_disczones = readmatrix("SM\report-surfarweavg-wallflux-shtc-mid-disc-rfile.out", FileType="text");

SM_vavgT_soliddisc =  readmatrix("SM\report-volavg-statemp-soliddisc-rfile.out", FileType="text");
SM_vmaxT_pads12 =  readmatrix("SM\report-volmax-statemp-pads12-rfile.out", FileType="text");

%% ---------------------------
% Variables list to average
% ----------------------------
variables = {
    'MRF1_CD_ty', 'MRF1_CD_tyrim', 'MRF1_CL_ty', 'MRF1_CL_tyrim', ...
    'MRF1_maxT_soliddisc', 'MRF1_heatflux_disczones', 'MRF1_heatflux_soliddisc', ...
    'MRF1_savgT_soliddisc', 'MRF1_savgV_vanes', 'MRF1_savg_HTC_soliddisc', ...
    'MRF1_savg_HTC_disczones', 'MRF1_vavgT_soliddisc', 'MRF1_vmaxT_pads12', ...
    'MRF2_CD_ty', 'MRF2_CD_tyrim', 'MRF2_CL_ty', 'MRF2_CL_tyrim', ...
    'MRF2_maxT_soliddisc', 'MRF2_heatflux_disczones', 'MRF2_heatflux_soliddisc', ...
    'MRF2_savgT_soliddisc', 'MRF2_savgV_vanes', 'MRF2_savg_HTC_soliddisc', ...
    'MRF2_savg_HTC_disczones', 'MRF2_vavgT_soliddisc', 'MRF2_vmaxT_pads12', ...
    'MRF3_CD_ty', 'MRF3_CD_tyrim', 'MRF3_CL_ty', 'MRF3_CL_tyrim', ...
    'MRF3_maxT_soliddisc', 'MRF3_heatflux_disczones', 'MRF3_heatflux_soliddisc', ...
    'MRF3_savgT_soliddisc', 'MRF3_savgV_vanes', 'MRF3_savg_HTC_soliddisc', ...
    'MRF3_savg_HTC_disczones', 'MRF3_vavgT_soliddisc', 'MRF3_vmaxT_pads12', ...
    'MRF4_CD_ty', 'MRF4_CD_tyrim', 'MRF4_CL_ty', 'MRF4_CL_tyrim', ...
    'MRF4_maxT_soliddisc', 'MRF4_heatflux_disczones', 'MRF4_heatflux_soliddisc', ...
    'MRF4_savgT_soliddisc', 'MRF4_savgV_vanes', 'MRF4_savg_HTC_soliddisc', ...
    'MRF4_savg_HTC_disczones', 'MRF4_vavgT_soliddisc', 'MRF4_vmaxT_pads12',...
    'SM_CD_ty', 'SM_CD_tyrim', 'SM_CL_ty', 'SM_CL_tyrim', ...
    'SM_maxT_soliddisc', 'SM_heatflux_disczones', 'SM_heatflux_soliddisc', ...
    'SM_savgT_soliddisc', 'SM_savgV_vanes', 'SM_savg_HTC_soliddisc', ...
    'SM_savg_HTC_disczones', 'SM_vavgT_soliddisc', 'SM_vmaxT_pads12'
};

%% ---------------------------
% Calculate averages
% - For MRF variables (steady-state) average last 2000 rows (column 2)
% - For SM variables (transient) average last 30 time-steps (column 2)
% ----------------------------
averages = struct();

for i = 1:numel(variables)
    var_name = variables{i};
    % get data variable from workspace
    data = eval(var_name);
    
    % Determine required number of rows depending on MRF or SM
    if startsWith(var_name, 'SM_')
        required_n = 30; % last 30 time-steps for SM
    else
        required_n = 2000; % last 2000 iterations for MRF
    end
    
    % Validate data size and choose column to average
    if isempty(data)
        warning('%s is empty. Average not calculated.', var_name);
        averages.(var_name) = NaN;
        continue;
    end
    
    nrows = size(data,1);
    ncols = size(data,2);
    
    % select which column to average: prefer column 2 if it exists, else last column
    colToUse = min(2, ncols);
    
    if nrows >= required_n
        start_idx = nrows - required_n + 1;
        values_to_avg = data(start_idx:nrows, colToUse);
        avg_val = mean(values_to_avg, 'omitnan'); % ignore NaNs
        averages.(var_name) = avg_val;
        fprintf('Average of last %d data points for %s: %.6g\n', required_n, var_name, avg_val);
    else
        % If there are fewer rows than required, warn and set NaN
        warning('%s has %d rows (< %d). Average not calculated.', var_name, nrows, required_n);
        averages.(var_name) = NaN;
    end
end

%% ---------------------------
% Convert averages struct to table and save results
% ----------------------------
averages_table = struct2table(averages); % one-row table with variable names as column headers
% transpose to get variable names in rows for easier saving/inspection
averages_table_vertical = rows2vars(averages_table); 

% Save to CSV
writetable(averages_table_vertical, 'average_results.csv');

%% ---------------------------
% Read average results back and organize by configuration
% ----------------------------
dataREAD = readtable('average_results.csv');

% variables: OriginalVariableNames (column) and Var1 (values)
name = dataREAD.OriginalVariableNames;
values = dataREAD.Var1;

% initialize structs
MRF1 = struct(); MRF2 = struct(); MRF3 = struct(); MRF4 = struct(); SM = struct();

for i = 1:length(name)
    current_name = name{i};
    current_value = values(i);
    
    if startsWith(current_name, 'MRF1_')
        field = strip(strrep(current_name, 'MRF1_', ''));
        MRF1.(field) = current_value;
    elseif startsWith(current_name, 'MRF2_')
        field = strip(strrep(current_name, 'MRF2_', ''));
        MRF2.(field) = current_value;
    elseif startsWith(current_name, 'MRF3_')
        field = strip(strrep(current_name, 'MRF3_', ''));
        MRF3.(field) = current_value;
    elseif startsWith(current_name, 'MRF4_')
        field = strip(strrep(current_name, 'MRF4_', ''));
        MRF4.(field) = current_value;
    elseif startsWith(current_name, 'SM_')
        field = strip(strrep(current_name, 'SM_', ''));
        SM.(field) = current_value;
    else
        warning('Unclassified variable: %s', current_name);
    end
end

% Convert structures to tables (will throw if no fields; guard with isempty)
if ~isempty(fieldnames(MRF1)), table_MRF1 = struct2table(MRF1); else table_MRF1 = table(); end
if ~isempty(fieldnames(MRF2)), table_MRF2 = struct2table(MRF2); else table_MRF2 = table(); end
if ~isempty(fieldnames(MRF3)), table_MRF3 = struct2table(MRF3); else table_MRF3 = table(); end
if ~isempty(fieldnames(MRF4)), table_MRF4 = struct2table(MRF4); else table_MRF4 = table(); end
if ~isempty(fieldnames(SM)),   table_SM   = struct2table(SM);   else table_SM   = table(); end

disp('MRF1 data:');
disp(table_MRF1);

%% ---------------------------
% Speed comparison: MRF2, MRF3, MRF1, MRF4, SM
% Note: ordering kept same as original comparison block
% ----------------------------
Method = ["MRF"; "MRF"; "MRF"; "MRF"; "SM"];

% Angular speeds in rad/s for each case (match the folders / headers)
w_omega = [31.7 , 44.1, 69.7, 95.1, 69.7]';  % MRF2, MRF3, MRF1, MRF4, SM

% The following indices (column numbers) assume that when the per-case tables
% were created the ordering of fields matches the previous script. If you find
% indexing errors here, print the tables and adapt indices to the correct columns.
try
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
catch ME
    warning('Indexing error while building comparison arrays: %s\nCheck the per-case tables to update column indices accordingly.', ME.message);
    % create empty arrays with NaN to allow the script to continue
    h = NaN(5,1); T_max = NaN(5,1); Q = NaN(5,1); Tavg = NaN(5,1); Cd_t = NaN(5,1); V_vanes = NaN(5,1);
end

% Create final comparison table
data = table(Method, w_omega, h, T_max, Q, Cd_t, Tavg, V_vanes, ...
    'VariableNames', {'Method', 'w_omega', 'h', 'T_max', 'Q','Cd_t', 'Tavg','V_vanes' });

% Save final data to CSV
writetable(data, 'data_avg.csv');

disp('Comparison table saved to data_avg.csv');