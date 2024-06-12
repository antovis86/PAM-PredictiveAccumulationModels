%% RDM analysis

load(fullfile('D:\PAM\RDM_HGF\RDM_HGF\results\SIM_threshold.mat'))
load(fullfile('D:\PAM\RDM_HGF\RDM_HGF\results\table_threshold.mat'))

%(fullfile('/home/ccnl/neuroscience/codici/HGF_DDM/SIM_results_lnr/SIM_01_mu.mat'))

k=0;
results_b = table;
for b0x = 1:length(b0B_all)
    for b1x = 1:length(b1_all)
        for Vvx = 1:length(Vv_all)
            for Vix = 1:length(Vi_all)
                k=k+1;

                results_b.b0B_sim(k) = b0B_all(b0x);
                results_b.b1_sim(k) = b1_all(b1x);
                results_b.Vv_sim(k) = Vv_all(Vvx);
                results_b.Vi_sim(k) = Vi_all(Vix);


                tmp = (median(p_obs(b0x,b1x,Vvx,Vix,:,:),5));

                results_b.b0B(k) = abs(tmp(1)-b0B_all(b0x))/b0B_all(b0x);
                results_b.b1(k) = abs(tmp(2)-b1_all(b1x))/b1_all(b1x);
                results_b.Vv(k) = abs(tmp(3)-Vv_all(Vvx))/Vv_all(Vvx);
                results_b.Vi(k) = abs(tmp(4)-Vi_all(Vix))/Vi_all(Vix);

                tmp_2 =  squeeze(median(p_prc(b0x,b1x,Vvx,Vix,:,:),5));

                results_b.om2(k) = tmp_2(1);
                results_b.om3(k) = tmp_2(2);


            end
        end
    end
end


%%
% Define colors for each parameter

% Identify unique combinations of the 4 parameters

paramNames = {'b0B', 'b1', 'Vv', 'Vi'};

uniqueCombinations = unique(models(:,paramNames));

numSubplots = size(uniqueCombinations);
numSubplots = numSubplots(1);

% Calculate the number of rows and columns for the subplot grid
numCols = ceil(sqrt(numSubplots));
numRows = ceil(numSubplots / numCols);

colors = lines(height(uniqueCombinations)); % Using MATLAB's default colormap for 4 distinct colors


% Loop over each parameter
for paramIdx = 1:length(paramNames)
    figure;
    % Create a tiled layout with 12 rows (one for each unique combination)
    t = tiledlayout(numRows, numCols, 'TileSpacing', 'Compact', 'Padding', 'Compact');

    for i = 1:height(uniqueCombinations)
        % Filter the table for the current combination of parameters
        filteredTable = models(...
            models.b0B == uniqueCombinations.b0B(i) & ...
            models.b1 == uniqueCombinations.b1(i) & ...
            models.Vv == uniqueCombinations.Vv(i) & ...
            models.Vi == uniqueCombinations.Vi(i), :);

        valuesPerParam = [];

        % Iterate over each row of the filtered table
        for row = 1:height(filteredTable)
            % Access the cell array in the current row
            cellArray = filteredTable.model{row};

            % Extract the variable from the cell array
            variableStruct = cellArray.p_obs.p; % Adjust this index if necessary

            % Collect the values for the current parameter
            valuesPerParam = [valuesPerParam; variableStruct(:, paramIdx)];
        end

        % Plot the histogram
        nexttile(i);
        h = histogram(valuesPerParam, 'Normalization', 'pdf', 'FaceColor', colors(paramIdx, :), 'EdgeColor', 'none');
        hold on;

        % Compute the KDE
        [f, xi, bw] = ksdensity(valuesPerParam);


        f_normalized = f / max(f);


        plot(xi, f, 'LineWidth', 2, 'Color', 'k');


        groundTruthValue = uniqueCombinations{i, paramIdx};
        xline(groundTruthValue, 'r', 'LineWidth', 2);


        % Update title with true values of the parameters
        titleStr = sprintf('b0B=%.2f, b1=%.2f, Vv=%.2f, Vi=%.2f', ...
            uniqueCombinations.b0B(i), uniqueCombinations.b1(i), uniqueCombinations.Vv(i), uniqueCombinations.Vi(i));
        title(titleStr);

        hold off;
    end

    % Adjust overall title for the parameter
    sgtitle(sprintf('Histograms for Parameter %s with KDE', paramNames{paramIdx}));
end



%% DDM analysis

load(fullfile('D:\PAM\DDM_HGF\RESULTS\ddmHGF_v_sim_results.mat'))

% add two columns of zero to keep track of the unresgistered parameters
models = [models(:, 1:2), array2table(zeros(height(models), 2), 'VariableNames', {'bw', 'ba'}), models(:, 3:end)];


%%
% Define colors for each parameter

% Identify unique combinations of the 4 parameters

paramNames = {'a', 'v', 'bw','ba','bv', 'T'};

uniqueCombinations = unique(models(:,paramNames));

numSubplots = size(uniqueCombinations);
numSubplots = numSubplots(1);

% Calculate the number of rows and columns for the subplot grid
numCols = ceil(sqrt(numSubplots));
numRows = ceil(numSubplots / numCols);

colors = lines(height(uniqueCombinations)); % Using MATLAB's default colormap for 4 distinct colors


% Loop over each parameter
for paramIdx = 1:length(paramNames)
    figure;
    % Create a tiled layout with 12 rows (one for each unique combination)
    t = tiledlayout(numRows, numCols, 'TileSpacing', 'Compact', 'Padding', 'Compact');

    for i = 1:height(uniqueCombinations)
        % Filter the table for the current combination of parameters
        filteredTable = models(...
            models.a == uniqueCombinations.a(i) & ...
            models.v == uniqueCombinations.v(i) & ...
            models.bv == uniqueCombinations.bv(i) & ...
            models.T == uniqueCombinations.T(i), :);

        valuesPerParam = [];

        % Iterate over each row of the filtered table
        for row = 1:height(filteredTable)
            % Access the cell array in the current row
            cellArray = filteredTable.model{row};

            % Extract the variable from the cell array
            variableStruct = cellArray.p_obs.p; % Adjust this index if necessary

            % Collect the values for the current parameter
            valuesPerParam = [valuesPerParam; variableStruct(:, paramIdx)];
        end

        % Plot the histogram
        nexttile(i);
        h = histogram(valuesPerParam, 'Normalization', 'pdf', 'FaceColor', colors(paramIdx, :), 'EdgeColor', 'none');
        hold on;

        % Compute the KDE
        [f, xi, bw] = ksdensity(valuesPerParam);


        plot(xi, f, 'LineWidth', 2, 'Color', 'k');


        groundTruthValue = uniqueCombinations{i, paramIdx};
        xline(groundTruthValue, 'r', 'LineWidth', 2);


        % Update title with true values of the parameters
        titleStr = sprintf('a=%.2f, v=%.2f, bv=%.2f, T=%.2f', ...
            uniqueCombinations.a(i), uniqueCombinations.v(i), uniqueCombinations.bv(i), uniqueCombinations.T(i));
        title(titleStr);

        hold off;
    end
    % Adjust overall title for the parameter
    sgtitle(sprintf('Histograms for Parameter %s with KDE', paramNames{paramIdx}));

end

%% LNR analysis
load(fullfile('D:\PAM\LNR_HGF\LNR_HGF\lnrHGF_sim_results.mat'))

%%
% Define colors for each parameter

% Identify unique combinations of the 4 parameters

paramNames = {'bv', 'bi', 'b1', 'sigma'};

uniqueCombinations = unique(models(:,paramNames));

numSubplots = size(uniqueCombinations);
numSubplots = numSubplots(1);

% Calculate the number of rows and columns for the subplot grid
numCols = ceil(sqrt(numSubplots));
numRows = ceil(numSubplots / numCols);

colors = lines(height(uniqueCombinations)); % Using MATLAB's default colormap for 4 distinct colors


% Loop over each parameter
for paramIdx = 1:length(paramNames)
    figure;
    % Create a tiled layout with 12 rows (one for each unique combination)
    t = tiledlayout(numRows, numCols, 'TileSpacing', 'Compact', 'Padding', 'Compact');

    for i = 1:height(uniqueCombinations)
        % Filter the table for the current combination of parameters
        filteredTable = models(...
            models.bv == uniqueCombinations.bv(i) & ...
            models.bi == uniqueCombinations.bi(i) & ...
            models.b1 == uniqueCombinations.b1(i) & ...
            models.sigma == uniqueCombinations.sigma(i), :);

        valuesPerParam = [];

        % Iterate over each row of the filtered table
        for row = 1:height(filteredTable)
            % Access the cell array in the current row
            cellArray = filteredTable.model{row};

            % Extract the variable from the cell array
            variableStruct = cellArray.p_obs.p; % Adjust this index if necessary

            % Collect the values for the current parameter
            valuesPerParam = [valuesPerParam; variableStruct(:, paramIdx)];
        end

        % Plot the histogram
        nexttile(i);
        h = histogram(valuesPerParam, 'Normalization', 'pdf', 'FaceColor', colors(paramIdx, :), 'EdgeColor', 'none');
        hold on;

        % Compute the KDE
        [f, xi, bw] = ksdensity(valuesPerParam);


        f_normalized = f / max(f);


        plot(xi, f, 'LineWidth', 2, 'Color', 'k');


        groundTruthValue = uniqueCombinations{i, paramIdx};
        xline(groundTruthValue, 'r', 'LineWidth', 2);


        % Update title with true values of the parameters
        titleStr = sprintf('bv=%.2f, bi=%.2f, b1=%.2f, sigma=%.2f', ...
            uniqueCombinations.bv(i), uniqueCombinations.bi(i), uniqueCombinations.b1(i), uniqueCombinations.sigma(i));
        title(titleStr);

        hold off;
    end

    % Adjust overall title for the parameter
    sgtitle(sprintf('Histograms for Parameter %s with KDE', paramNames{paramIdx}));
end
