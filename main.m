clc; clear; close all;

%% get DATA from xlsx files

% Node to Node Resistances
data_nodes_resistances = readmatrix('input-node-to-node-resistances.xlsx');

nodes_start = data_nodes_resistances(:,1);
nodes_end = data_nodes_resistances(:,2);
thermal_resis = data_nodes_resistances(:,3);

% Power losses of each Node
data_nodes_power_losses = readmatrix('input-nodes-power-losses.xlsx');

power_nodes = data_nodes_power_losses(:,2);

% Capacitance of each Node
data_nodes_capacitances = readmatrix('input-nodes-capacitances.xlsx');

capacitance_nodes = data_nodes_capacitances(:,2);

% Fixed Nodes & Temperatures

data_fixed_nodes = readmatrix('input-fixed-nodes-temperatures.xlsx');

fixed_nodes = data_fixed_nodes(:,1);

fixed_nodes_temp = data_fixed_nodes(:,2);



%% modify imported DATA



start_end = [nodes_start,nodes_end];
unique_nodes = unique(start_end(:).')';

unique_nodes = sort(unique_nodes);


% supress Cooling water Jacket power loss % capacitance

for j =1:length(fixed_nodes)

    fixedIndex = find(unique_nodes == fixed_nodes(j));

    power_nodes(fixedIndex) = [];

    capacitance_nodes(fixedIndex) = [];

end

capacitance_nodes = diag(capacitance_nodes);


%% Auto Assemble Resistance Matrix

unknown_nodes = unique_nodes;

unknown_nodes ( ismember(unique_nodes,fixed_nodes)) = [];


resisMatrix = zeros(length(unknown_nodes));



 % constVector = zeros(length(unknown_nodes),1);

 constVector = power_nodes ;


for i = 1 : length(nodes_start)
    % constIndex = 0;
    startIndex = find(unknown_nodes == nodes_start(i));
    endIndex = find(unknown_nodes == nodes_end(i));
    % assign zero to empty double
    if isempty(startIndex)
        startIndex =0;
    end
    if isempty(endIndex)
        endIndex =0;
    end

    % assign Resistance values to NONDIAGONAL unknown Elements
    if ~(ismember(nodes_start(i),fixed_nodes) || ismember(nodes_end(i),fixed_nodes))
        % constIndex = constIndex+1;

        resisMatrix(startIndex,endIndex) = thermal_resis(i)^-1;
        resisMatrix(endIndex,startIndex) = thermal_resis(i)^-1;

        % assign resis values to DIAGONAL elements
        resisMatrix(startIndex,startIndex) = resisMatrix(startIndex,startIndex) - thermal_resis(i)^-1;
        resisMatrix(endIndex,endIndex) = resisMatrix(endIndex,endIndex) - thermal_resis(i)^-1;
    else
        if ismember(startIndex,fixed_nodes)
            constIndex = endIndex ;
            fixedIndex = find(fixed_nodes==startIndex);
        else
            constIndex = startIndex ;
            fixedIndex = find(fixed_nodes==endIndex);
        end
        
        

        % assign Generated Power by fixed point and its Resistance
        constVector(constIndex) = fixed_nodes_temp(fixedIndex) *thermal_resis(i)^-1;
        % assign resis values to DIAGONAL elements
        resisMatrix(constIndex,constIndex) = resisMatrix(constIndex,constIndex) - thermal_resis(i)^-1;
    end

end


% CHECK LOGIC convert "zero" capacitances to "-Inf"

%(**************************** CHECK THE LOGIC **********************

for ii = 1:length(capacitance_nodes(:,1))
    if capacitance_nodes(ii,ii) == 0
        capacitance_nodes(ii,ii) = -Inf;
    end
end

%*********************************************************************)

%% STEADY-STATE SOLUTION

T_steady = -resisMatrix\constVector;

%% add fixed nodes to the final solution

T_all = T_steady;

align = [-1,-1];

for j =1:length(fixed_nodes)

   fixedIndex = find(unique_nodes == fixed_nodes(j));

   T_all = insert(fixed_nodes_temp(j),T_all,fixedIndex,align(j));

   

end

%% get the Nodes Lables and Make the Result TAble as xlsx

data_nodes_labels = readtable('input-nodes-labels.xlsx');

results_table = [data_nodes_labels(:,2),table(T_all)];

results_table. Properties. VariableNames = {'Nodes Labels','Steady Temperatures [°C]'};

writetable(results_table,'output-steady-temperatures.xlsx')

