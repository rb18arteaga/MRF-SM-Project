clc
clear all

labelsData = ["Cp" "cell_element_type"	"cell_id"	"cell_type"	"cell_zone"	"density"	"dynamic_pressure"	"expr:C_D"	"expr:C_D_contact"	"expr:C_Lift"	"expr:C_Lift_contact"	"helicity"	"laminar_kinetic_energy"	"pressure"	"pressure_coefficient"	"q_criterion"	"raw_q_criterion"	"skin_friction_coef"	"strain_rate_mag"	"total_pressure"	"turb_diss_rate"	"turb_intensity"	"turb_kinetic_energy"	"turb_reynolds_number_rey"	"velocity:0"	"velocity:1"	"velocity:2"	"velocity_magnitude"	"viscosity_lam"	"viscosity_turb"	"vorticity_mag"	"wall_shear"	"x_velocity"	"y_plus"	"y_velocity"	"z_velocity"	"vtkValidPointMask"	"Points:0"	"Points:1"	"Points:2"]

tyre_Angle = linspace(0,360,361)';
mediumOldData = readmatrix("MediumOld_tyre-only-data.csv");
MRF_Data = readmatrix("MRF-tyre-only-data.csv");
SM_Data = readmatrix("SM-tyre-only-data.csv");
%fineData = readmatrix("FIne-tyre-only85-data.csv");


Cp_mediumOld_Data = mediumOldData(:,1);
Cp_MRF_Data = MRF_Data(:,1);
% Crear un índice lógico para excluir ceros en y
indices = Cp_MRF_Data ~= 0;   % "~=" significa "no igual a cero"
% Filtrar los vectores
x_filtrado = tyre_Angle(indices);
y_filtrado = Cp_MRF_Data(indices);
Cp_SM_Data = SM_Data(:,1);
%Cp_fineData = fineData(:,1);

exp_mears = readmatrix("mears-Exp-Data.csv");
Cp_exp_mears = exp_mears(:,2)
exp_tyre_angle = exp_mears(:,1)


%visualisation NGS normalized grid spacing
figure(1);
plot(exp_tyre_angle,Cp_exp_mears,'ko','MarkerSize',3,'MarkerFaceColor','k');
hold on;
plot(tyre_Angle,Cp_mediumOld_Data,'r' , LineWidth=1.5);
hold on;
hold on;
plot(x_filtrado,y_filtrado,'g', LineWidth=2);
hold on;
plot(tyre_Angle,Cp_SM_Data, 'b' ,LineWidth=1.5);
xlabel('Theta','Interpreter', 'latex'); 
ylabel('Cp', 'Interpreter', 'latex');
title('Mesh Methods', 'Interpreter', 'latex'); 
legend('Cp-Exp-Mears','Cp-Medium','Cp-MRF','Cp-SM', 'Location', 'North', 'Interpreter', 'latex'); 
grid on;
%xticks(h);               
%yticks(sort(L)); 

