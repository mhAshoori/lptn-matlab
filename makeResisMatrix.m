function modified_resistance_matrix = ...
    makeResisMatrix...
    (modifiedIndexVector,startEndResisMatrix,fixed_nodes)

%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
%************* READ THIS ******************
% REORDER THE START,END,RESIS(64,3) MATRIX AND GIVE IT TO THE 
%*****************************************


%%%%%%%%%%%%%%%% THIS PART IS WRONG %%%%%%%%%%%%%%%%%%%%
%%% IT SHOULD REORDER THE 119 RESIS BUT MAKES 64 RESIS NUMBERS 

% modifiedStartEndResis = startEndResisMatrix(modifiedIndexVector,:);

%%%%%%%%% THE UPPER LINE IS SUPPRESSED AS IT'S WRONG %%%%%%%

nodes_start = startEndResisMatrix(:,1);
nodes_end = startEndResisMatrix(:,2);

thermal_resis = startEndResisMatrix(:,3);

%% modify imported DATA

% merge start & end nodes
start_end = [nodes_start,nodes_end];
unique_nodes = unique(start_end(:).')';

% No Sorting is Needed
% sort(unique_nodes);

% supress fixed nodes: cooling water jacket, power loss & capacitance

% for j =1:length(fixed_nodes)

  %  fixedIndex = find(unique_nodes == fixed_nodes(j));

   % power_nodes(fixedIndex) = [];

  %  capacitance_nodes_vector(fixedIndex) = [];

% end




%% Auto Assemble Resistance Matrix

unknown_nodes = unique_nodes;

unknown_nodes ( ismember(unique_nodes,fixed_nodes)) = [];

%%% INSTEAD WE RE ORDER 

unknown_nodes = unknown_nodes(modifiedIndexVector,:);



resisMatrix = zeros(length(unknown_nodes));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Assemble Loop

modified_resistance_matrix = zeros(length(unknown_nodes));

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

        modified_resistance_matrix(startIndex,endIndex) = thermal_resis(i)^-1;
        modified_resistance_matrix(endIndex,startIndex) = thermal_resis(i)^-1;

        % assign resis values to DIAGONAL elements
        modified_resistance_matrix(startIndex,startIndex) = resisMatrix(startIndex,startIndex) - thermal_resis(i)^-1;
        modified_resistance_matrix(endIndex,endIndex) = resisMatrix(endIndex,endIndex) - thermal_resis(i)^-1;
    else
        if ismember(startIndex,fixed_nodes)
            constIndex = endIndex ;
            fixedIndex = find(fixed_nodes==startIndex);
        else
            constIndex = startIndex ;
            fixedIndex = find(fixed_nodes==endIndex);
        end

        % assign Generated Power by fixed point and its Resistance
       %%%%%THIS LINE IS SUPRESSED, IT'S MADE OUT SIDE OF THE FUNCTION %%%%%%%%
       % constVector(constIndex) = fixed_nodes_temp(fixedIndex) *thermal_resis(i)^-1;
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
        % assign resis values to FIXED DIAGONAL elements
        modified_resistance_matrix(constIndex,constIndex) = resisMatrix(constIndex,constIndex) - thermal_resis(i)^-1;
    end

end


end

