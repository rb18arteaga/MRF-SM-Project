clc
clear all

% ===== Load experimental data =====
data5d8 = readmatrix("5d8.csv");
data6d8 = readmatrix("6d8.csv");
data7d8 = readmatrix("7d8.csv");

x5d8 = data5d8(:,1);  v5d8 = data5d8(:,2);
x6d8 = data6d8(:,1);  v6d8 = data6d8(:,2);
x7d8 = data7d8(:,1);  v7d8 = data7d8(:,2);

% ===== Moving Wall - Medium Mesh =====
medium5d8 = readmatrix("medium5d8-data.csv");
medium6d8 = readmatrix("medium6d8-data.csv");
medium7d8 = readmatrix("medium7d8-data.csv");

xmed = medium5d8(:,5);

v_medium5d8 = medium5d8(:,1);
indexM5 = v_medium5d8 ~= 0;
x_M_5 = xmed(indexM5);  v_M_5 = v_medium5d8(indexM5);

v_medium6d8 = medium6d8(:,1);
indexM6 = v_medium6d8 ~= 0;
x_M_6 = xmed(indexM6);  v_M_6 = v_medium6d8(indexM6);

v_medium7d8 = medium7d8(:,1);
indexM7 = v_medium7d8 ~= 0;
x_M_7 = xmed(indexM7);  v_M_7 = v_medium7d8(indexM7);

% ===== MRF Mesh =====
MRF5d8 = readmatrix("MRF5d8.csv");
MRF6d8 = readmatrix("MRF6d8.csv");
MRF7d8 = readmatrix("MRF7d8.csv");

xMRF = MRF5d8(:,5);

v_MRF5d8 = MRF5d8(:,1);
indexMRF5 = v_MRF5d8 ~= 0;
x_MRF_5 = xMRF(indexMRF5);  v_MRF_5 = v_MRF5d8(indexMRF5);

v_MRF6d8 = MRF6d8(:,1);
indexMRF6 = v_MRF6d8 ~= 0;
x_MRF_6 = xMRF(indexMRF6);  v_MRF_6 = v_MRF6d8(indexMRF6);

v_MRF7d8 = MRF7d8(:,1);
indexMRF7 = v_MRF7d8 ~= 0;
x_MRF_7 = xMRF(indexMRF7);  v_MRF_7 = v_MRF7d8(indexMRF7);

% ===== Sliding Mesh =====
SM5d8 = readmatrix("SM5d8.csv");
SM6d8 = readmatrix("SM6d8.csv");
SM7d8 = readmatrix("SM7d8.csv");

xSM = SM5d8(:,5);

v_SM5d8 = SM5d8(:,1);
indexSM5 = v_SM5d8 ~= 0;
x_SM_5 = xSM(indexSM5);  v_SM_5 = v_SM5d8(indexSM5);

v_SM6d8 = SM6d8(:,1);
indexSM6 = v_SM6d8 ~= 0;
x_SM_6 = xSM(indexSM6);  v_SM_6 = v_SM6d8(indexSM6);

v_SM7d8 = SM7d8(:,1);
indexSM7 = v_SM7d8 ~= 0;
x_SM_7 = xSM(indexSM7);  v_SM_7 = v_SM7d8(indexSM7);

%% ===== Single figure with 3 subplots =====
figure;
set(gcf, 'Position', [100 100 1200 400]);

% ---- Subplot 1: X = 5D/8 ----
subplot(1,3,1);
plot(x5d8, v5d8, 'ko', 'MarkerSize', 3, 'MarkerFaceColor', 'k'); hold on;
plot(v_M_5, x_M_5, 'b-', 'LineWidth', 2);
plot(v_MRF_5, x_MRF_5, 'r--', 'LineWidth', 2);
plot(v_SM_5, x_SM_5, 'g-.', 'LineWidth', 2);
xlabel('$V_x$', 'Interpreter', 'latex');
ylabel('$Z$', 'Interpreter', 'latex');
title('$X = \frac{5D}{8}$', 'Interpreter', 'latex');
legend({'Exp', 'Moving Wall', 'MRF', 'Sliding Mesh'}, ...
       'Location', 'north', 'Interpreter', 'latex');

% ---- Subplot 2: X = 6D/8 ----
subplot(1,3,2);
plot(x6d8, v6d8, 'ko', 'MarkerSize', 3, 'MarkerFaceColor', 'k'); hold on;
plot(v_M_6, x_M_6, 'b-', 'LineWidth', 2);
plot(v_MRF_6, x_MRF_6, 'r--', 'LineWidth', 2);
plot(v_SM_6, x_SM_6, 'g-.', 'LineWidth', 2);
xlabel('$V_x$', 'Interpreter', 'latex');
ylabel('$Z$', 'Interpreter', 'latex');
title('$X = \frac{6D}{8}$', 'Interpreter', 'latex');

% ---- Subplot 3: X = 7D/8 ----
subplot(1,3,3);
plot(x7d8, v7d8, 'ko', 'MarkerSize', 3, 'MarkerFaceColor', 'k'); hold on;
plot(v_M_7, x_M_7, 'b-', 'LineWidth', 2);
plot(v_MRF_7, x_MRF_7, 'r--', 'LineWidth', 2);
plot(v_SM_7, x_SM_7, 'g-.', 'LineWidth', 2);
xlabel('$V_x$', 'Interpreter', 'latex');
ylabel('$Z$', 'Interpreter', 'latex');
title('$X = \frac{7D}{8}$', 'Interpreter', 'latex');

% ---- Overall title ----
sgtitle('$V_x$ Comparison: Exp / Moving Wall / MRF / Sliding Mesh', 'Interpreter', 'latex');