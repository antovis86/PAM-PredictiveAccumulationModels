% Extract the number of rows of the table
numRows = height(models);

% Define the group size
groupSize = nsubj;

% Initalize the table 
resultsTable = table(); 

% Loop through the rows in groups of 50
for startRow = 1:groupSize:numRows
    endRow = min(startRow+groupSize-1, numRows); % Ensure endRow doesn't exceed numRows
    % Extract the rows for the current group
    groupData = models(startRow:endRow, :);

    % extract parameters values
    a = groupData.a(1);
    sigma = groupData.sigma(1);
    f = groupData.f(1);

    % Populate a temporary table
    temp_table = table();
    temp_table.a = a;
    temp_table.sigma = sigma;
    temp_table.f = f;
    temp_table.corr = {tapas_bayesian_parameter_average_av(groupData.model)};
    
    % Add the column 
    resultsTable = [resultsTable; temp_table]; % Append a row of NaNs

    % Do something with the result, e.g., display it or store it
    
    % You can also access the indices of the current group:
    fprintf('Processing rows %d to %d\n', startRow, endRow);
end
