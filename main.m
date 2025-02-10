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

capacitance_nodes_vector = data_nodes_capacitances(:,2);

% Fixed Nodes & Temperatures

data_fixed_nodes = readmatrix('input-fixed-nodes-temperatures.xlsx');

fixed_nodes = data_fixed_nodes(:,1);

fixed_nodes_temp = data_fixed_nodes(:,2);



%% modify imported DATA

% merge start & end nodes
start_end = [nodes_start,nodes_end];
unique_nodes = unique(start_end(:).')';

% sort nodes
 sort(unique_nodes);


% supress fixed nodes: cooling water jacket, power loss & capacitance

for j =1:length(fixed_nodes)

    fixedIndex = find(unique_nodes == fixed_nodes(j));

    power_nodes(fixedIndex) = [];

    capacitance_nodes_vector(fixedIndex) = [];

end




%% Auto Assemble Resistance Matrix

unknown_nodes = unique_nodes;

unknown_nodes ( ismember(unique_nodes,fixed_nodes)) = [];


resisMatrix = zeros(length(unknown_nodes));

 % constVector = zeros(length(unknown_nodes),1);

 constVector = power_nodes ;

for i = 1 : length(nodes_start)
  
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


%% STEADY-STATE SOLUTION

T_steady = -resisMatrix\constVector;

%% add fixed nodes to the final Steady solution

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

%% Prepare DATA for Transient Solution

nonzero_capacitannce_nodes = capacitance_nodes_vector;
ii = 0;
zero_index = double.empty;
nonzero_index = double.empty;
for i = 1:length(capacitance_nodes_vector(:,1))
    ii = ii+1;
    if capacitance_nodes_vector(i) == 0
       nonzero_capacitannce_nodes(ii) = [];
     %  nonzero_capacitannce_nodes(ii) = [];
       ii=ii-1;
       % find the index of ZERO capacitance elements
       zero_index = [zero_index;i];
    else
        % find the index of NONZERO capacitance elements
        nonzero_index = [nonzero_index;i];
    end
end

% REORDER matrixes of Resistance, Power Loss & Capacitance

modified_resistance_matrix = resisMatrix([nonzero_index;zero_index],:);
modified_const_vector = constVector([nonzero_index;zero_index],:);
modified_capacitances_vector = capacitance_nodes_vector([nonzero_index;zero_index],:);

% convert capacitance vector to matrix
modified_capacitance_matrix = diag(modified_capacitances_vector);

C1 = modified_capacitance_matrix(1:ii,1:ii);

A1 = modified_resistance_matrix(1:ii,1:ii);
A2 = modified_resistance_matrix(1:ii,ii+1:end);
A3 = modified_resistance_matrix(ii+1:end,1:ii);
A4 = modified_resistance_matrix(ii+1:end,ii+1:end);

B1 = modified_const_vector(1:ii);
B2 = modified_const_vector(ii+1:end);