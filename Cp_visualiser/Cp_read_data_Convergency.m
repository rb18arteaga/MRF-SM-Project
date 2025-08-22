clear all
clc

labelsData = ["Cp" "cell_element_type"	"cell_id"	"cell_type"	"cell_zone"	"density"	"dynamic_pressure"	"expr:C_D"	"expr:C_D_contact"	"expr:C_Lift"	"expr:C_Lift_contact"	"helicity"	"laminar_kinetic_energy"	"pressure"	"pressure_coefficient"	"q_criterion"	"raw_q_criterion"	"skin_friction_coef"	"strain_rate_mag"	"total_pressure"	"turb_diss_rate"	"turb_intensity"	"turb_kinetic_energy"	"turb_reynolds_number_rey"	"velocity:0"	"velocity:1"	"velocity:2"	"velocity_magnitude"	"viscosity_lam"	"viscosity_turb"	"vorticity_mag"	"wall_shear"	"x_velocity"	"y_plus"	"y_velocity"	"z_velocity"	"vtkValidPointMask"	"Points:0"	"Points:1"	"Points:2"]

tyre_Angle = linspace(0,360,361)';
coarseData = readmatrix("Coarse2-tyre-only-data.csv");
%mediumNewData = readmatrix("MediumNew-tyre-only-data.csv");
mediumOldData = readmatrix("MediumOld_tyre-only-data.csv");
fineData = readmatrix("FIne-tyre-only85-data.csv");

Cp_coarseData = coarseData(:,1);
%Cp_mediumNew_Data = mediumNewData(:,1);
Cp_mediumOld_Data = mediumOldData(:,1);
Cp_fineData = fineData(:,1);

exp_mears = readmatrix("mears-Exp-Data.csv");
Cp_exp_mears = exp_mears(:,2)
exp_tyre_angle = exp_mears(:,1)


%visualisation NGS normalized grid spacing
figure(1);
plot(exp_tyre_angle,Cp_exp_mears,'ko','MarkerSize',3,'MarkerFaceColor','k');
hold on;
plot(tyre_Angle,Cp_coarseData,'r' , LineWidth=1.5);
hold on;
%plot(tyre_Angle,Cp_mediumNew_Data,'-k',LineWidth=1.5);
hold on;
plot(tyre_Angle,Cp_mediumOld_Data,'g', LineWidth=1.5);
hold on;
plot(tyre_Angle,Cp_fineData, 'b' ,LineWidth=1.5);
xlabel('Theta','Interpreter', 'latex'); 
ylabel('C\_p', 'Interpreter', 'latex');
title('Mesh Convergence Analysis', 'Interpreter', 'latex'); 
legend('Cp-Exp-Mears','Cp-Coarse','Cp-Medium','Cp-Fine', 'Location', 'North', 'Interpreter', 'latex'); 
%legend('Cp-Exp-Mears','Cp-Coarse','Cp-MediumOld','Cp-MediumNEW','Cp-Fine', 'Location', 'North', 'Interpreter', 'latex'); 
grid on;
%xticks(h);               
%yticks(sort(L)); 

