clc
clear all

% Read experimental data
data5d8 = readmatrix("5d8.csv");
data6d8 = readmatrix("6d8.csv");
data7d8 = readmatrix("7d8.csv");

x5d8 = data5d8(:,1);
x6d8 = data6d8(:,1);
x7d8 = data7d8(:,1);

v5d8 = data5d8(:,2);
v6d8 = data6d8(:,2);
v7d8 = data7d8(:,2);

%% Coarse data
coarse5d8 = readmatrix("coarse5d8-data.csv");
coarse6d8 = readmatrix("coarse6d8-data.csv");
coarse7d8 = readmatrix("coarse7d8-data.csv");

xcor = coarse5d8(:,5);

v_coarse5d8 = coarse5d8(:,1);
indexC5 = v_coarse5d8 ~= 0;
x_C_5 = xcor(indexC5);
v_C_5 = v_coarse5d8(indexC5);

v_coarse6d8 = coarse6d8(:,1);
indexC6 = v_coarse6d8 ~= 0;
x_C_6 = xcor(indexC6);
v_C_6 = v_coarse6d8(indexC6);

v_coarse7d8 = coarse7d8(:,1);
indexC7 = v_coarse7d8 ~= 0;
x_C_7 = xcor(indexC7);
v_C_7 = v_coarse7d8(indexC7);

%% Medium Data
medium5d8 = readmatrix("medium5d8-data.csv");
medium6d8 = readmatrix("medium6d8-data.csv");
medium7d8 = readmatrix("medium7d8-data.csv");

xmed = medium5d8(:,5);

v_medium5d8 = medium5d8(:,1);
indexM5 = v_medium5d8 ~= 0;
x_M_5 = xmed(indexM5);
v_M_5 = v_medium5d8(indexM5);

v_medium6d8 = medium6d8(:,1);
indexM6 = v_medium6d8 ~= 0;
x_M_6 = xmed(indexM6);
v_M_6 = v_medium6d8(indexM6);

v_medium7d8 = medium7d8(:,1);
indexM7 = v_medium7d8 ~= 0;
x_M_7 = xmed(indexM7);
v_M_7 = v_medium7d8(indexM7);

%% Fine Data
fine5d8 = readmatrix("fine5d8-data.csv");
fine6d8 = readmatrix("fine6d8-data.csv");
fine7d8 = readmatrix("fine7d8-data.csv");

xfine = fine5d8(:,5);

v_fine5d8 = fine5d8(:,1);
index5 = v_fine5d8 ~= 0;
x_f_5 = xfine(index5);
v_f_5 = v_fine5d8(index5);

v_fine6d8 = fine6d8(:,1);
index6 = v_fine6d8 ~= 0;
x_f_6 = xfine(index6);
v_f_6 = v_fine6d8(index6);

v_fine7d8 = fine7d8(:,1);
index7 = v_fine7d8 ~= 0;
x_f_7 = xfine(index7);
v_f_7 = v_fine7d8(index7);

%% ===== Single figure with 3 subplots =====
figure;
set(gcf, 'Position', [100 100 1200 400]);

% ---- Subplot 1: X = 5D/8 ----
subplot(1,3,1);
plot(x5d8, v5d8, 'ko', 'MarkerSize', 3, 'MarkerFaceColor', 'k'); hold on;
plot(v_C_5, x_C_5, 'r-', 'LineWidth', 2);
plot(v_M_5, x_M_5, 'b-', 'LineWidth', 2);
plot(v_f_5, x_f_5, 'g-', 'LineWidth', 2);
xlabel('$V_x$', 'Interpreter', 'latex');
ylabel('$Z$', 'Interpreter', 'latex');
title('$X = \frac{5D}{8}$', 'Interpreter', 'latex');
legend({'Exp', 'Coarse', 'Medium', 'Fine'}, ...
       'Location', 'north', 'Interpreter', 'latex');

% ---- Subplot 2: X = 6D/8 ----
subplot(1,3,2);
plot(x6d8, v6d8, 'ko', 'MarkerSize', 3, 'MarkerFaceColor', 'k'); hold on;
plot(v_C_6, x_C_6, 'r-', 'LineWidth', 2);
plot(v_M_6, x_M_6, 'b-', 'LineWidth', 2);
plot(v_f_6, x_f_6, 'g-', 'LineWidth', 2);
xlabel('$V_x$', 'Interpreter', 'latex');
ylabel('$Z$', 'Interpreter', 'latex');
title('$X = \frac{6D}{8}$', 'Interpreter', 'latex');

% ---- Subplot 3: X = 7D/8 ----
subplot(1,3,3);
plot(x7d8, v7d8, 'ko', 'MarkerSize', 3, 'MarkerFaceColor', 'k'); hold on;
plot(v_C_7, x_C_7, 'r-', 'LineWidth', 2);
plot(v_M_7, x_M_7, 'b-', 'LineWidth', 2);
plot(v_f_7, x_f_7, 'g-', 'LineWidth', 2);
xlabel('$V_x$', 'Interpreter', 'latex');
ylabel('$Z$', 'Interpreter', 'latex');
title('$X = \frac{7D}{8}$', 'Interpreter', 'latex');

% ---- Overall title ----
sgtitle('Mesh Convergence Analysis: Coarse / Medium / Fine', 'Interpreter', 'latex');