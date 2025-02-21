% SIDE CALCULATION

imported_resis_matrix - data_nodes_resistances;

ref_temp = readmatrix('input-ref-temperatures.xlsx');

T_combined = [T_all_steady,ref_temp(:,2)];

ref_temp(:,1) = [];

ref_temp(1) = [];
ref_temp(28) = [];



residual_calculated = abs(constVector-resisMatrix*T_steady);

residual_ref = abs(constVector - resisMatrix*ref_temp);

norm_calc = norm(residual_calculated)

norm_ref = norm(residual_ref)

norm_calc - norm_ref

T_SOR_HW = SOR_HW(-resisMatrix,constVector,70*ones(64,1),0.5);

power_nodes_ref_temps = -resisMatrix*ref_temp;

powers_combined = [power_nodes_ref_temps,constVector];

residual_calculated = norm(constVector-resisMatrix*T_steady)/norm(constVector)
resis_cond = cond(resisMatrix)
det_resisMatrix = det(resisMatrix)

T_combined(:,1) - T_combined(:,2)

% NOT WORKING => cgs(resisMatrix,constVector)




