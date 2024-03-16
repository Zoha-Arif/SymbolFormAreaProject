%clear all; close all; clc;
%% ====================================================================== %%

% Set working directories.
rootdir = '/Volumes/LANDLAB/projects/sfa/';
blprojectid = 'proj-63e6b75dc538c16a8225593e-majortracts'; 
mp2rage = 'no';
wmmeasure = 'fa';

%generate column of tracts of interest ids
colorProfiles = '/Volumes/LANDLAB/projects/sfa/supportFiles/colorProfilesTOIS.csv';
colorProfiles = readtable(colorProfiles);
%% Correlate silhouette scores with math and reading scores

%insert local path of Tshort.csv and Tlong.csv file
d = readtable(fullfile(rootdir, 'supportfiles', ['TlongTOIS' '.csv']));
xyz = load(fullfile(rootdir, 'supportfiles', ['toisXYZ' '.mat']));
%convert struct into table; the table is named 'data' from the import
%script
xyz = xyz.data; 

nfig = 0; 
%% ====================================================================== %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Behavioral Data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

math = d.T2_WJIV_AppProb_raw - d.T1_WJIV_AppProb_raw;
read = d.T2_WJIV_LWID_raw - d.T1_WJIV_LWID_raw;

%%%%%%%%%%%%%% kmeans & silhouette score for various k-values %%%%%%%%%%%%%% 
% Define range of k values to try
kValues = 2;
% Iterate over different values of k

%generate column of tracts of interest ids
tractIDs = {'RNFA-RVWFA', 'LNFA-LVWFA'}; 

varNames = {'Subject', 'Tract', 'k', 'kmeans', 'silhouette', 'idx', 'data', 'source'}; 
varTypes = ["double", "char", "double", "cell", "double", "cell", "cell", "cell"]; 
numRows = length(varNames) * 4 * length(tractIDs); 
numCols = length(varNames); 
%create an empty table
optkTBL = table('Size',[numRows, numCols], 'VariableNames', varNames, 'VariableTypes', varTypes); 

rowNum = 1; 

clear silhouette; 

for i = 1:length(kValues)
    k = kValues(i);
    for t = 1:length(tractIDs)
    
        % Split the tract ID based on the '-' delimiter
        parts = strsplit(tractIDs{t}, '-');
    
        % Extract the tract names using regular expression
        tract1 = strjoin(parts(1));
        tract2 = strjoin(parts(2));
        
        tbl1 = table(xyz.subjectID, xyz.(tract1), 'VariableNames', {'Subject', 'Data'}); 
        tbl2 = table(xyz.subjectID, xyz.(tract2), 'VariableNames', {'Subject', 'Data'});
        
        tbl1.Source = repmat({tract1(2:end)}, size(tbl1, 1), 1);
        tbl2.Source = repmat({tract2(2:end)}, size(tbl2, 1), 1);
        
        tbl = vertcat(tbl1, tbl2);

        coordArray = []; 
        numRows = size(tbl, 1); 
        resultSubjectID = {};
        sourceArray = {}; 
        
        for f = 1:numRows
            % Extract the double array from the current cell
            currentDoubleArray = tbl.Data{f};
            
            % Extract subject ID from the table
            subjectID = tbl.Subject(f);  % Assuming 'SubjectID' is the name of your subject ID column
            
            sources = tbl.Source(f); 
            
            % Check if currentDoubleArray is not empty
            if ~isempty(currentDoubleArray) && size(currentDoubleArray, 2) == 3
                % Concatenate XYZ coordinates along the 1st dimension
                coordArray = [coordArray; currentDoubleArray];
                % Concatenate the subject ID
                resultSubjectID = [resultSubjectID; repmat(subjectID, size(currentDoubleArray, 1), 1)];
                %Concatenate the source array
                sourceArray = [sourceArray; repmat(sources, size(currentDoubleArray, 1), 1)];
            else
                warning('Array does not have 3 columns or is empty. Skipping...');
            end
        end
        
        coordTBL = table(resultSubjectID, coordArray(:,1), coordArray(:,2), coordArray(:,3), sourceArray, 'VariableNames', {'Subject', 'x_coords'... 
            'y_coords', 'z_coords', 'source'});
        coordTBL = coordTBL(~any(ismissing(coordTBL), 2), :);
        coordTBL = sortrows(coordTBL); 
        
        % Convert each cell element to a string
        charSubject = cellfun(@(x) num2str(x), coordTBL.Subject, 'UniformOutput', false);
        % Replace the column in the table with the char data
        coordTBL.Subject = charSubject;
        
        subIDs = unique(coordTBL.Subject); 
        for j = 1:size(subIDs)
            % Perform k-means clustering for the current participant
            
            %get x-coordinates
            xIDX = find(ismember(coordTBL.Subject, subIDs(j))); 
            x = coordTBL.x_coords(xIDX);
            
            %get y-coordinates
            y = coordTBL.y_coords(xIDX);
            
            %get z-coordinates
            z = coordTBL.z_coords(xIDX);
            
            %combine into a matrix
            data = [x, y, z];
            
            sourceCur = coordTBL.source(xIDX);
            
            % Perform k-means clustering
            [idx, centroids] = kmeans(data, k);
            
            % Calculate silhouette scores for each datapoint
            silhouette_values = silhouette(data, idx);
        
            % Check for NaN values in silhouette_values
            if any(isnan(silhouette_values))
                warning('Silhouette values contain NaNs. Try adjusting the number of clusters or preprocessing your data.');
            end

            % Calculate the average silhouette score
            average_silhouette_score = mean(silhouette_values);
            
            %Also get the average silhouette score for each cluster
            %in each clustering, keep track of which cluster each datapoint
            %belogns too; calculate what percentage of datapoints are in
            %each cluster
            
            C = num2cell([centroids(:,1), centroids(:,2), centroids(:,3)]);
            C = {cat(1, C{:})}; 
            
            %Save the results
            optkTBL.Subject(rowNum) = str2double(subIDs(j));
            optkTBL.Tract(rowNum) = {char(tractIDs(t))}; 
            optkTBL.k(rowNum) = k; 
            optkTBL.kmeans(rowNum) = C; 
            optkTBL.silhouette(rowNum) = average_silhouette_score; 
            optkTBL.idx(rowNum) = {idx}; 
            optkTBL.data(rowNum) = {data}; 
            optkTBL.source(rowNum) = {sourceCur};
            rowNum = rowNum + 1; 
            
            %clear idx;   
            %clear centroids; 
            
        end
    end
end

%%%%%%%%%%%%%% Find Optimal K-Value & Plot K-Means %%%%%%%%%%%%%% 
varNames = {'Subject', 'Tract', 'k', 'silhouette'}; 
varTypes = ["double", "char", "double", "double"]; 
numRows = length(varNames) * length(tractIDs); 
numCols = length(varNames); 

%create an empty table
bestk = table('Size',[numRows, numCols], 'VariableNames', varNames, 'VariableTypes', varTypes); 

rowNum = 1; 
subIDs = unique(optkTBL.Subject); 
for s = 1:length(subIDs)
    for g = 1:length(tractIDs)
        %get subject
        subIDX = find((optkTBL.Subject == subIDs(s)) & (strcmpi(optkTBL.Tract, {char(tractIDs(g))})));
        
        if ~(isempty(subIDX))
            tempTBL = table(); 
            tempTBL.k = optkTBL.k(subIDX);
            tempTBL.silhouette = optkTBL.silhouette(subIDX);
            tempTBL.kmeans = optkTBL.kmeans(subIDX);
            tempTBL.idx = optkTBL.idx(subIDX);
            tempTBL.data = optkTBL.data(subIDX);
            tempTBL.source = optkTBL.source(subIDX);
            
            %find the maximum silhouette score
            [maxSil, maxIDX] = max(tempTBL.silhouette); 
            rowWithMaxSilScore = tempTBL(maxIDX, :); 
        
            %save in table
            bestk.Subject(rowNum) = subIDs(s);
            bestk.Tract(rowNum) = {char(tractIDs(g))};
            bestk.k(rowNum) =  rowWithMaxSilScore.k; 
            bestk.silhouette(rowNum) =  rowWithMaxSilScore.silhouette; 
            rowNum = rowNum + 1; 
            
            %%%%%%%% PLOTS %%%%%%%%
            %print out figure for every subject so that it can be
            %visualized.
            %Plot the clusters
            %get x-coordinates
            
            data = rowWithMaxSilScore.data{1}; 
            sourceData = rowWithMaxSilScore.source{1};
            
            %centroids
            % Extract x, y, and z coordinates
            k = rowWithMaxSilScore.k; 
            centroids = rowWithMaxSilScore.kmeans{1}; 
            
            % Reshape the input array
            centroids = reshape(centroids, k, []);
            
            %group number
            idx = rowWithMaxSilScore.idx{1};
            colors = lines(k); 
            
            nfig = nfig + 1; 
            f = figure(nfig); 
            
            combinedData = table(data, string(sourceData), 'VariableNames', {'data', 'sourceData'}); 
            
            for v = 1:k
                if (v == 2)
                    hold on
                end
                tempData = combinedData(idx == v, :); 
                clusterData = combinedData.data(idx == v, :); 
                
                % Determine brightness for current cluster based on source
                brightness = ones(size(clusterData, 1), 1); % Initialize brightness
                for p = 1:size(clusterData, 1)
                    if strcmp(tempData.sourceData(p), "NFA")
                        brightness(p) = 1; % Light source
                    else
                        brightness(p) = 0.5; % Dark source
                    end
                end
                
                scatter3(clusterData(:, 1), clusterData(:,2), clusterData(:,3), [], brightness .* colors(v, :), 'filled'); 
            end
            
            clear v; 
            neonColors = min(colors + 0.5, 1); 
            for v = 1:k
                scatter3(centroids(v, 1), centroids(v, 2), centroids(v, 3), 200, neonColors(v, :), 'filled', 'Marker', 'o', 'MarkerEdgeColor', 'k',...
                    'LineWidth', 2);
            end
            
            title(['K-Means Clustering for ' num2str(subIDs(s)) ' for ' char(tractIDs(g))]);
            xlabel('X');
            ylabel('Y');
            zlabel('Z');
            hold off;
            
            %%%%%%%% END OF PLOTS %%%%%%%%
            
            %export figure as a png file
            mainpath = '/Volumes/LANDLAB/projects/sfa/supportFiles/final-plots-kmeans';
            filename = strjoin(['silHem_kmeans_' string(subIDs(s)) '_' tractIDs(g)]); % Define a filename for each plot
            fullFilePath = fullfile(mainpath, filename); % Create the full file path
            saveas(nfig, fullFilePath, 'png'); % Save the figure as a PNG file
            
        end
     end
end

%%%%%%%%%%%%%% Correlate Silhouette Scores with Math & Reading %%%%%%%%%%%%%% 
rsqTBLMath = table(tractIDs'); 

nrow = size(rsqTBLMath, 1); 

%Re-establishing rsqTBL
rsqTBLMath.beta = zeros(nrow, 1); 
rsqTBLMath.n = zeros(nrow, 1); 
rsqTBLMath.p_value = zeros(nrow, 1); 
rsqTBLMath.rmse = zeros(nrow, 1); 

%markerColor
markerColor = lines(length(tractIDs));

for a = 1:length(tractIDs)
    
    nfig = nfig + 1;
    f = figure(nfig);

    %startingx, startingy, width height
    f.Position = [1000 1000 800 700];
    
    hold on;
    
    %==================== Math ====================
    %plotting a linear regression model 
    MathTBL = table(d.subjectID, math, 'VariableNames', {'Subject', 'Data'});  
    MathTBL = MathTBL(~any(ismissing(MathTBL), 2), :);
    MathTBL = unique(MathTBL, 'rows'); 
    MathTBL = sortrows(MathTBL); 
    
    toiIDX = (strcmpi(bestk.Tract, {char(tractIDs(a))}));
    silTBL = table(bestk.Subject(toiIDX), bestk.silhouette(toiIDX), 'VariableNames', {'Subject', 'Data'});  
    
    %Merge Tables
    tbl = innerjoin(MathTBL, silTBL, 'Keys', 'Subject');
    
    %Rename Headers
    tbl.Properties.VariableNames{'Data_MathTBL'} = 'Math';
    tbl.Properties.VariableNames{'Data_silTBL'} = 'silhouette';
    
    tbl(any(ismissing(tbl), 2), :) = [];
    
    %z-score xVar and Math    
    silhouette = zscore(tbl.silhouette, [], 'omitnan');
    Math = zscore(tbl.Math, [], 'omitnan');
    
    % Subsetting tbl into tb so that outlier removals from this
    % point on do not affect all of the other tracts and behavioral
    % measures contained in tbl.
    tb = array2table([silhouette Math], 'VariableNames', {'silhouette', 'Math'});

    % Remove outliers, extreme z-scores, on either measure.
    %idx_remove = unique([find(abs(tb.xVar) >= 2.5); find(abs(tb.Math) >= 2.5)]);
    %if ~isempty(idx_remove); tb(idx_remove, :) = []; end
    %clear idx_remove
    
    %%defining the line to fit the model to
    Q = 'Math ~ silhouette';
    
    %generating the model
    mdl = fitlm(tb, Q);
    
    %get appropriate RGB color for tract by indexing into colorProfiles.csv
    %idx = find(strcmp(colorProfiles.NameOfTrack, char(tractIDs(t))) == 1);
    %markerColor = [colorProfiles.Red(idx)/255, colorProfiles.Green(idx)/255, colorProfiles.Blue(idx)/255];

    clf; 
    
    %plotting the model
    h = plot(mdl, 'Marker', 'o', 'MarkerEdgeColor', 'white', 'MarkerFaceColor', markerColor(a, :), 'MarkerSize', 12);
    
    %grab trendline and datapoints
    %dataHandle = findobj(h,'DisplayName','Data');
    %fitHandle = findobj(h,'DisplayName','Fit');
    dataHandle = h(1); 
    fitHandle = h(2); 

    %fill in confidence interval
    cbHandles = findobj(h,'DisplayName','Confidence bounds');
    cbHandles = findobj(h,'LineStyle',cbHandles.LineStyle, 'Color', cbHandles.Color);
   
    upperCBHandle = h(4);
    lowerCBHandle = h(3);
    
    xData = upperCBHandle.XData; 
    k = patch([xData xData(end:-1:1) xData(1)], [lowerCBHandle.YData upperCBHandle.YData(end:-1:1) lowerCBHandle.YData(1)], 'b');
    set(k, 'EdgeColor', 'none', 'FaceColor', [markerColor(a, 1)*0.55  markerColor(a, 2)*0.55 markerColor(a, 3)*0.55], 'FaceAlpha', '0.2')

    %style the trendline
    set(fitHandle, 'Color', [markerColor(a, 1) markerColor(a, 2) markerColor(a, 3)], 'LineWidth', 3)

    %w = plot(mdl, 'Marker', 'o', 'MarkerEdgeColor', 'white', 'MarkerFaceColor', markerColor, 'MarkerSize', 12);
    set(cbHandles, 'visible', 'off'); 
    pltLeg = legend('', '', '');
    set(pltLeg,'visible','off')
    %fitHandle2 = findobj(w,'DisplayName','Fit');
    %set(fitHandle2, 'Visible', 'off')

    %adding title and color to the model
    plotTitle = {char(tractIDs(t))};
    plotTitle = strjoin(['Simple Linear Model for', plotTitle]);
    title(plotTitle);
    
    n = size(Math, 1);
   
    rsqTBLMath.beta(a) = mdl.Coefficients{2,1}; 
    rsqTBLMath.n(a) = n; 
    rsqTBLMath.p_value(a) = mdl.Coefficients.pValue(2); 
    rsqTBLMath.rmse(a) = mdl.RMSE; 

    %set scale of y-axis
    %ylim([0.3 0.6])

    %add adjusted r squared to table.
    %rsqTableAdj.SimpleLin(t) = mdl.Rsquared.Adjusted;
    %rsqTableOrd.SimpleLin(t) = mdl.Rsquared.Ordinary;
    %aicTable.SimpleLin(t)= mdl.ModelCriterion.AIC;
    
    %===========================================================================
    % Set up plot and measure-specific details.
    capsize = 0;
    marker = 'o';
    linewidth = 1.5;
    linestyle = 'none';
    markersize = 100;
    
    fontname = 'Arial';
    fontsize = 50;
    fontangle = 'italic';
    
    yticklength = 0;
    xticklength = 0.02;

    % xaxis
    xax = get(gca, 'xaxis');
    xax.TickDirection = 'out';
    xax.TickLength = [xticklength xticklength];
    %minX = floor(min(xData) * 100) / 100;
    %maxX = floor(max(xData) * 100) / 100;
    %midX = floor(((min(xData) + max(xData))/2) * 100) / 100; 
    %set(gca, 'XLim', [minX maxX], 'XTick', [minX midX maxX]);
    set(gca, 'XLim', [-3 3], 'XTick', [-3 -2 -1 0 1 2 3]);
    xax.FontName = fontname;
    xax.FontSize = fontsize;

    % yaxis
    yax = get(gca,'yaxis');
    yax.TickDirection = 'out';
    yax.TickLength = [yticklength yticklength];
    set(gca, 'YLim', [-3 3], 'YTick', [-3 -2 -1 0 1 2 3]);
    yax.FontName = fontname;
    yax.FontSize = fontsize;
    yax.FontAngle = fontangle;
    
    %%%%%%To print on plots
    text(-2.8, 2.2, ['rmse = ' num2str(mdl.RMSE)])
    text(-2.8, 2.0, ['beta = ' num2str(mdl.Coefficients{2,1})])
    text(-2.8, 1.8, ['p = ' num2str(mdl.Coefficients.pValue(2))])
    text(-2.8, 1.6, ['n = ' num2str(n)]);

    %change figure background to white
    set(gcf, 'color', 'w')

    %===========================================================================
    
    hold off
    
    %adding title and color to the model
    plotTitle = {char(tractIDs(a))};
    plotTitle = strjoin(['Simple Linear Model for', plotTitle]);
    title(plotTitle);
    xlabel('Silhouette Scores');
    ylabel('Math Scores');
    
    %export figure as a png file
    mainpath = '/Volumes/LANDLAB/projects/sfa/supportFiles/final-plots-kmeans';
    filename = strjoin(['silHem_Math_', tractIDs(a)]); % Define a filename for each plot
    fullFilePath = fullfile(mainpath, filename); % Create the full file path
    saveas(nfig, fullFilePath, 'png'); % Save the figure as a PNG file
    
end

% ============================== Reading ============================== 
rsqTBLRead = table(tractIDs'); 

nrow = size(rsqTBLRead, 1); 

%Re-establishing rsqTBL
rsqTBLRead.beta = zeros(nrow, 1); 
rsqTBLRead.n = zeros(nrow, 1); 
rsqTBLRead.p_value = zeros(nrow, 1); 
rsqTBLRead.rmse = zeros(nrow, 1); 

for a = 1:length(tractIDs)
    
    nfig = nfig + 1;
    f = figure(nfig);

    %startingx, startingy, width height
    f.Position = [1000 1000 800 700];
    
    hold on;
    
    %==================== Read ====================
    %plotting a linear regression model 
    ReadTBL = table(d.subjectID, read, 'VariableNames', {'Subject', 'Data'});  
    ReadTBL = ReadTBL(~any(ismissing(ReadTBL), 2), :);
    ReadTBL = unique(ReadTBL, 'rows'); 
    ReadTBL = sortrows(ReadTBL); 
    
    toiIDX = (strcmpi(bestk.Tract, {char(tractIDs(a))}));
    silTBL = table(bestk.Subject(toiIDX), bestk.silhouette(toiIDX), 'VariableNames', {'Subject', 'Data'});  
    
    %Merge Tables
    tbl = innerjoin(ReadTBL, silTBL, 'Keys', 'Subject');
    
    %Rename Headers
    tbl.Properties.VariableNames{'Data_ReadTBL'} = 'Read';
    tbl.Properties.VariableNames{'Data_silTBL'} = 'silhouette';
    
    tbl(any(ismissing(tbl), 2), :) = [];
    
    %z-score xVar and Math    
    silhouette = zscore(tbl.silhouette, [], 'omitnan');
    Read = zscore(tbl.Read, [], 'omitnan');
    
    % Subsetting tbl into tb so that outlier removals from this
    % point on do not affect all of the other tracts and behavioral
    % measures contained in tbl.
    tb = array2table([silhouette Read], 'VariableNames', {'silhouette', 'Read'});

    % Remove outliers, extreme z-scores, on either measure.
    %idx_remove = unique([find(abs(tb.xVar) >= 2.5); find(abs(tb.Math) >= 2.5)]);
    %if ~isempty(idx_remove); tb(idx_remove, :) = []; end
    %clear idx_remove
    
    %%defining the line to fit the model to
    Q = 'Read ~ silhouette';
    
    %generating the model
    mdl = fitlm(tb, Q);
    
    %get appropriate RGB color for tract by indexing into colorProfiles.csv
    %idx = find(strcmp(colorProfiles.NameOfTrack, char(tractIDs(t))) == 1);
    %markerColor = [colorProfiles.Red(idx)/255, colorProfiles.Green(idx)/255, colorProfiles.Blue(idx)/255];

    clf; 
    
    %plotting the model
    h = plot(mdl, 'Marker', 'o', 'MarkerEdgeColor', 'white', 'MarkerFaceColor', markerColor(a, :), 'MarkerSize', 12);
    
    %grab trendline and datapoints
    %dataHandle = findobj(h,'DisplayName','Data');
    %fitHandle = findobj(h,'DisplayName','Fit');
    dataHandle = h(1); 
    fitHandle = h(2); 

    %fill in confidence interval
    cbHandles = findobj(h,'DisplayName','Confidence bounds');
    cbHandles = findobj(h,'LineStyle',cbHandles.LineStyle, 'Color', cbHandles.Color);
   
    upperCBHandle = h(4);
    lowerCBHandle = h(3);
    
    xData = upperCBHandle.XData; 
    k = patch([xData xData(end:-1:1) xData(1)], [lowerCBHandle.YData upperCBHandle.YData(end:-1:1) lowerCBHandle.YData(1)], 'b');
    set(k, 'EdgeColor', 'none', 'FaceColor', [markerColor(a, 1)*0.55  markerColor(a, 2)*0.55 markerColor(a, 3)*0.55], 'FaceAlpha', '0.2')

    %style the trendline
    set(fitHandle, 'Color', [markerColor(a, 1) markerColor(a, 2) markerColor(a, 3)], 'LineWidth', 3)

    %w = plot(mdl, 'Marker', 'o', 'MarkerEdgeColor', 'white', 'MarkerFaceColor', markerColor, 'MarkerSize', 12);
    set(cbHandles, 'visible', 'off'); 
    pltLeg = legend('', '', '');
    set(pltLeg,'visible','off')
    %fitHandle2 = findobj(w,'DisplayName','Fit');
    %set(fitHandle2, 'Visible', 'off')

    %adding title and color to the model
    plotTitle = {char(tractIDs(t))};
    plotTitle = strjoin(['Simple Linear Model for', plotTitle]);
    title(plotTitle);
    
    n = size(Read, 1);
   
    rsqTBLRead.beta(a) = mdl.Coefficients{2,1}; 
    rsqTBLRead.n(a) = n; 
    rsqTBLRead.p_value(a) = mdl.Coefficients.pValue(2); 
    rsqTBLRead.rmse(a) = mdl.RMSE; 

    %set scale of y-axis
    %ylim([0.3 0.6])

    %add adjusted r squared to table.
    %rsqTableAdj.SimpleLin(t) = mdl.Rsquared.Adjusted;
    %rsqTableOrd.SimpleLin(t) = mdl.Rsquared.Ordinary;
    %aicTable.SimpleLin(t)= mdl.ModelCriterion.AIC;
    
    %===========================================================================
    % Set up plot and measure-specific details.
    capsize = 0;
    marker = 'o';
    linewidth = 1.5;
    linestyle = 'none';
    markersize = 100;
    
    fontname = 'Arial';
    fontsize = 50;
    fontangle = 'italic';
    
    yticklength = 0;
    xticklength = 0.02;

    % xaxis
    xax = get(gca, 'xaxis');
    xax.TickDirection = 'out';
    xax.TickLength = [xticklength xticklength];
    %minX = floor(min(xData) * 100) / 100;
    %maxX = floor(max(xData) * 100) / 100;
    %midX = floor(((min(xData) + max(xData))/2) * 100) / 100; 
    %set(gca, 'XLim', [minX maxX], 'XTick', [minX midX maxX]);
    set(gca, 'XLim', [-3 3], 'XTick', [-3 -2 -1 0 1 2 3]);
    xax.FontName = fontname;
    xax.FontSize = fontsize;

    % yaxis
    yax = get(gca,'yaxis');
    yax.TickDirection = 'out';
    yax.TickLength = [yticklength yticklength];
    set(gca, 'YLim', [-3 3], 'YTick', [-3 -2 -1 0 1 2 3]);
    yax.FontName = fontname;
    yax.FontSize = fontsize;
    yax.FontAngle = fontangle;
    
    %%%%%%To print on plots
    text(-2.8, 2.2, ['rmse = ' num2str(mdl.RMSE)])
    text(-2.8, 2.0, ['beta = ' num2str(mdl.Coefficients{2,1})])
    text(-2.8, 1.8, ['p = ' num2str(mdl.Coefficients.pValue(2))])
    text(-2.8, 1.6, ['n = ' num2str(n)]);

    %change figure background to white
    set(gcf, 'color', 'w')

    %===========================================================================
    
    hold off
    
    %adding title and color to the model
    plotTitle = {char(tractIDs(a))};
    plotTitle = strjoin(['Simple Linear Model for', plotTitle]);
    title(plotTitle);
    xlabel('Silhouette Scores');
    ylabel('Read Scores');
    
    %export figure as a png file
    mainpath = '/Volumes/LANDLAB/projects/sfa/supportFiles/final-plots-kmeans';
    filename = strjoin(['silHem_Read_', tractIDs(a)]); % Define a filename for each plot
    fullFilePath = fullfile(mainpath, filename); % Create the full file path
    saveas(nfig, fullFilePath, 'png'); % Save the figure as a PNG file
    
end
