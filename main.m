clc; clear; close all;

%% get DATA from xlsx files

% Node to Node Resistances
data_nodes_resistances = readmatrix('input-node-to-node-resistances.xlsx');

% THERMAL RESIS VALUES IS IMPORTTED THROUGH .RMF FILE SEE LINE [ 42 - 60 ]
% so 3rd lines commented

nodes_start = data_nodes_resistances(:,1);
nodes_end = data_nodes_resistances(:,2);
% % thermal_resis = data_nodes_resistances(:,3);

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


% merge start & end nodes
start_end = [nodes_start,nodes_end];
unique_nodes = unique(start_end(:).')';

% sort nodes
 sort(unique_nodes);


%% [NEW] get RESIS MATRIX from .RMF matrix data of MOTOR-CAD

data_rmf_resis_matrix = readmatrix('input-rmf-resis-matrix.xlsx');

% % % data_rmf_resis_matrix == data_rmf_resis_matrix';
% % % 
% % % if ~issymmetric (data_rmf_resis_matrix)
% % %     error('imported resistance matrix is NOT SYMMETRIC')
% % % end

k = 0;
imported_resis_matrix = zeros(1);
for i = 1:length(data_rmf_resis_matrix(1,:))
    for j = i:length(data_rmf_resis_matrix(:,1))
        if ~ismember(data_rmf_resis_matrix(i,j),[0,1e9,2e9,3e9,4e9])
            k = k+1;
            imported_resis_matrix(k,1) = unique_nodes(i);
            imported_resis_matrix(k,2) = unique_nodes(j);
            imported_resis_matrix(k,3) = data_rmf_resis_matrix(i,j);
        end
    end
end


nodes_start = imported_resis_matrix(:,1);
nodes_end = imported_resis_matrix(:,2);
thermal_resis = imported_resis_matrix(:,3);

%%% imported_resis_matrix - data_nodes_resistances


%% modify imported DATA

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

% T_all_steady includes the fixed nodes too
T_all_steady = T_steady;

align = [-1,-1];

for j =1:length(fixed_nodes)

   fixedIndex = find(unique_nodes == fixed_nodes(j));

   T_all_steady = insert(fixed_nodes_temp(j),T_all_steady,fixedIndex,align(j));

end

%% get the Nodes Lables and Make the Steady Result Table as .xlsx

data_nodes_labels = readtable('input-nodes-labels.xlsx');

results_table_steady = [data_nodes_labels(:,2),table(T_all_steady)];

results_table_steady. Properties. VariableNames = {'Nodes Labels','Steady Temperatures [°C]'};

writetable(results_table_steady,'output-steady-temperatures.xlsx')

%% Prepare DATA for Transient Solution

global A B

if ismember(0,capacitance_nodes_vector)

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

% REORDER matrixes of Power Loss & Capacitance

%%%% [THIS LINE IS WRONG] modified_resistance_matrix = resisMatrix([nonzero_index;zero_index],:);
modified_const_vector = constVector([nonzero_index;zero_index],:);
modified_capacitances_vector = capacitance_nodes_vector([nonzero_index;zero_index],:);

% convert capacitance vector to matrix
modified_capacitance_matrix = diag(modified_capacitances_vector);

C1 = modified_capacitance_matrix(1:ii,1:ii);

% modify (reorder) Constant vector

B1 = modified_const_vector(1:ii);
B2 = modified_const_vector(ii+1:end);

% MAKE MODIFIED resis matrix

modified_resistance_matrix = ...
    makeResisMatrix...
    ([nonzero_index;zero_index],data_nodes_resistances,fixed_nodes);

%%%%%%%%%%% this line a just added to check the idea of replacing
%%%%%%%%%%% the DIAGONAL OF RESISMATRIX WITH MODIFIED ONE  %%%%%%%%%%%

mainResisMatrixDiag = diag(resisMatrix);

modifiedResisMatrixDiag = mainResisMatrixDiag([nonzero_index;zero_index],:);

for i = 1:length(resisMatrix(:,1))
    modified_resistance_matrix(i,i) = modifiedResisMatrixDiag(i);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%% THE ABOVE IDEA SEEMS TO BE WORKING AS THE DATA IS PHYSICALLY RIGHT
%%%%%% %%%
A1 = modified_resistance_matrix(1:ii,1:ii);
A2 = modified_resistance_matrix(1:ii,ii+1:end);
A3 = modified_resistance_matrix(ii+1:end,1:ii);
A4 = modified_resistance_matrix(ii+1:end,ii+1:end);

A = C1\(A1-A2*(A4\A3));
 
B = C1\(B1-(A2*(A4\B2)));


else

    A = diag(capacitance_nodes_vector) \  resisMatrix;

    B = diag(capacitance_nodes_vector) \ constVector;
end

%% Transient Solution



% B = C1\(B1) ;

% Initialization and Solution

T0 = ones(ii,1)*70;

[t,T_transient_nonzero] = ode45(@lptn_ode,[0 10000],T0);

% find zero capacitance nodes Temperatures
T_transient_zero = -A4\A3*T_transient_nonzero';

T_transient_zero = T_transient_zero';

% combine NONZERO & ZERO Solutions of final timespan
T_transient_END = [T_transient_nonzero(end,:)';T_transient_zero(end,:)'];

%T_transient_END = T_transient_END([nonzero_index;zero_index],:);

modified_index = [nonzero_index;zero_index];

stored_indexes = double.empty;

for i =1 :length(modified_index)

    if i ~=modified_index(i) && ~ismember(i,stored_indexes)
        temp1 = T_transient_END(i);
        temp2 = T_transient_END(modified_index(i));
        stored_indexes = [stored_indexes;i,modified_index(i)];
        T_transient_END(i) = temp2;
        T_transient_END(modified_index(i)) = temp1;
        
    end

end

%T_transient_END = [T_transient_END(2:end);T_transient_END(1)];

%% add fixed nodes to the final Steady solution

% T_all_transient includes the fixed nodes too
T_all_transient = T_transient_END;

align = [-1,-1];

for j =1:length(fixed_nodes)

   fixedIndex = find(unique_nodes == fixed_nodes(j));

   T_all_transient = insert(fixed_nodes_temp(j),T_all_transient,fixedIndex,align(j));

end

 
%% get the Nodes Lables and Make the TRANSIENT Results Table as .xlsx

% THIS VARIABLE IS ALLREADY DEFINED 
% data_nodes_labels = readtable('input-nodes-labels.xlsx');

results_table_transient = [data_nodes_labels(:,2),table(T_all_transient)];

results_table_transient . Properties. VariableNames = {'Nodes Labels','End of TimeSpan Transient Temperatures [°C]'};

writetable(results_table_transient ,'output-transient-temperatures.xlsx')

%% calculate the difference between Transient and Steady Solution

%%% THIS PART WILL BE IMPLEMENTED LATER
% % % % max(T_all_steady-T_all_transient)
% % % 
% % % sol_combined = [T_all_steady,T_all_transient,abs((T_all_transient-T_all_steady)./T_all_steady)*100];
% % % 
% % % results_table_combined = array2table([data_nodes_labels(:,2),table(sol_combined)]);
% % % 
% % % results_table_combined . Properties. VariableNames = {'Nodes Labels','Steady Temperatures [°C]','Transient Temperatures [°C]','Error [%]'};
% % % 
% % % writetable(results_table_combined ,'output-combined-temperatures-error.xlsx')


