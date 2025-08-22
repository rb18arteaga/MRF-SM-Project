%% Temporal convergence
vol_domain = 36.11; %m3
%             h0,   h1,        h2          h3
mesh_size = [inf, 85423959, 15558127, 3818019]

h = (vol_domain./mesh_size).^(1/3);

C_D = [0, 0.412, 0.402, 0.377]
C_L = [0, 0.07, 0.13, 0.063];


% mesh_refinement = [1/inf, 0.025,0.5,1];
% dt_original = time(2)-time(1);
% dt_values = dt_original.*mesh_refinement;
% std_acc_values = zeros(1,3); % To store standard deviations
% for i = 2:length(mesh_refinement) 
%     %Begins in 2 as time begins in 1/inf 
%     dtm = dt_values(i); % Current time step    
%     time_new = 0:dtm:max(time); % New time vector   
%     w_resampled = interp1(time, w_point, time_new, 'linear'); 
%     % Linear data interpolation    
%     acc_local = gradient(w_resampled, dtm); 
%     % Local acceleration with the new time step    
%     std_acc_values(i) = std(acc_local); % Standard deviation of acceleration
% end

L = C_D;
r21 = h(3)/h(2)
r32 = h(4)/h(3)

tol_wanted = 1e-10;
tol = 1;
q_old = 0.001;
cter = 1;

 
while tol > tol_wanted    
    p_est = (1 / log(r21)) * abs(log(abs((L(4) - L(3)) / (L(3) - L(2)))) + q_old); 
    
    q_new = log((r21^p_est - sign((L(4) - L(3)) / (L(3) - L(2)))) / ...        
        (r32^p_est - sign((L(4) - L(3)) / (L(3) - L(2)))));    
    tol = abs(q_new - q_old);     
    q_old = q_new;       
    cter = cter + 1; 
end
p = p_est 

% Richardson extrapolation para \(\Delta t \to 0\) 
L(1) = (r21^p * L(2) - L(3)) / (r21^p - 1)
% Relative errors
eps21 = abs((L(2) - L(1)) / L(1)) 
eps32 = abs((L(3) - L(2)) / L(2)) 
% Calculus of GCIs 
GCI21 = (1.25 * eps21 / (r21^p - 1))*100 
GCI32 = (1.25 * eps32 / (r32^p - 1))*100
% Asymptotic range
AR = (r21^(p)) * GCI21 / GCI32 


%visualisation NGS normalized grid spacing
figure(9);
plot(h,L,'-ko','MarkerSize',8,'MarkerFaceColor','k');
hold on;
plot(h(1)/h(2),L(1),'rd','MarkerSize',10,'MarkerFaceColor','r');
xlabel('$\Delta t$','Interpreter', 'latex'); 
ylabel('C_D', 'Interpreter', 'latex');
title('Mesh Convergence Analysis', 'Interpreter', 'latex'); 
legend('Computed Values', 'Extrapolated Value ($\Delta t \to 0$)', 'Location', 'SouthWest', 'Interpreter', 'latex'); 
grid on;
xticks(h);               
yticks(sort(L)); 

% %Results of downsample convergence
% disp('Results of downsample convergence analysis:');
% disp(table(downsample_rates', std_acc_downsampled', 'VariableNames', {'Downsampled Rates', 'Std_Acc_Downsampled'}));
% 
% %Results of temporal convergence analysis:
% disp('Results of temporal convergence analysis:');
% disp(table(h', L', 'VariableNames', {'Delta_t', 'Std_Acc_Local'}));
% disp(['Extrapolated value (Δt → 0): ' num2str(L(1))]); 
% disp(['Estimated order of convergence (p): ', num2str(p)]);
% disp(['GCI21: ', num2str(GCI21), ', GCI32: ', num2str(GCI32)]); disp(['Asymptotic Range (AR): ', num2str(AR)]);
% 
% disp('%%%%%')