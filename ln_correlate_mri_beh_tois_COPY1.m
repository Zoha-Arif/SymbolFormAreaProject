%clear all; close all; clc;

%% ====================================================================== %%
% White matter measure
wmmeasure = 'fa'; % 'md', 'ad', 'md', 'rd', 'ndi', 'isovf', 'odi', 'map', 'T1', 'R1';

% Hemisphere
% This code will generate both hemispheres

% Set working directories.
rootdir = '/Volumes/LANDLAB/projects/sfa/';
blprojectid = 'proj-63e6b75dc538c16a8225593e-majortracts'; 
mp2rage = 'no';

%generate column of tracts of interest ids
nfig = 30; 
colorProfiles = '/Volumes/LANDLAB/projects/sfa/supportFiles/colorProfilesTOIS.csv';
colorProfiles = readtable(colorProfiles);
%% Correlate fa with math and reading scores

%insert local path of Tshort.csv and Tlong.csv file
d = readtable(fullfile(rootdir, 'supportfiles', ['TshortTOIS' '.csv']));
Tlong = readtable(fullfile(rootdir, 'supportfiles', ['TlongTOIS' '.csv']));

%% ====================================================================== %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Behavioral Data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

math = d.T2_WJIV_AppProb_raw - d.T1_WJIV_AppProb_raw;
read = d.T2_WJIV_LWID_raw - d.T1_WJIV_LWID_raw;
%% ====================================================================== %%

%============== Generate Math Plots ==============

%generate column of tracts of interest ids
tractIDs = {'lhvwfap', 'lhnfap', 'rhvwfap', 'rhnfap', 'lhvwfao', 'lhnfao', 'rhvwfao', 'rhnfao'}; 

rsqSimpleLinMath = table(tractIDs'); 

nrow = size(rsqSimpleLinMath, 1); 

%Re-establishing rsqSimpleLinRead
rsqSimpleLinMath.beta = zeros(nrow, 1); 
rsqSimpleLinMath.n = zeros(nrow, 1); 
rsqSimpleLinMath.p_value = zeros(nrow, 1); 
rsqSimpleLinMath.rmse = zeros(nrow, 1); 

for t = 1:length(tractIDs)
    
    found = false; 
    
    while all((strcmp(d.Properties.VariableNames, tractIDs{t})) == 0)
        t = t + 1; 
        if (t > length(tractIDs))
            found = true;  
            break; 
        end
    end
    
    if found
        break;
    end
    
    nfig = nfig + 1;
    f = figure(nfig); hold on;

    %startingx, startingy, width height
    f.Position = [1000 1000 800 700];
    
    hold on 
    
    %==================== Math ====================
    %plotting a nonlinear aggression model 
    MathTBL = table(d.subIDs, math, 'VariableNames', {'Subject', 'Data'});  
    MathTBL = MathTBL(~any(ismissing(MathTBL), 2), :);
    MathTBL = sortrows(MathTBL); 
    
    xVarTBL = table(d.subIDs, d.(char(tractIDs(t))), 'VariableNames', {'Subject', 'Data'});
    xVarTBL(xVarTBL.Data < 0.3, :)= [];
    xVarTBL = xVarTBL(~any(ismissing(xVarTBL), 2), :);
    xVarTBL = sortrows(xVarTBL); 
   
    sexTBL = table(d.subIDs, d.Sex, 'VariableNames', {'Subject', 'Data'});
    sexTBL = sexTBL(~any(ismissing(sexTBL), 2), :);
    sexTBL = sortrows(sexTBL); 
    
    %Merge Tables
    tblMerge = innerjoin(sexTBL, xVarTBL, 'Keys', 'Subject');
    tbl = innerjoin(tblMerge, MathTBL, 'Keys', 'Subject');
    
    %Rename Headers
    tbl.Properties.VariableNames{'Data_sexTBL'} = 'Sex';
    tbl.Properties.VariableNames{'Data_xVarTBL'} = 'xVar';
    tbl.Properties.VariableNames{'Data'} = 'Math';
    
    tbl(any(ismissing(tbl), 2), :) = [];
    
    Math = tbl.Math; 
    xVar = tbl.xVar; 
    Sex = tbl.Sex; 
    
    %Check correlation between Sex and Math
    P = 'Math ~ Sex'; 
    mdl1 = fitlm(tbl, P);
      
    %z-score xVar and Math    
    xVar = zscore(xVar, [], 'omitnan');
    Math = zscore(Math, [], 'omitnan');

    %%defining the line to fit the model to
    Q = 'Math ~ xVar';

    %generating the model
    mdl = fitlm(tbl, Q);
    
    %get appropriate RGB color for tract by indexing into colorProfiles.csv
    idx = find(strcmp(colorProfiles.NameOfTrack, char(tractIDs(t))) == 1);
    markerColor = [colorProfiles.Red(idx)/255, colorProfiles.Green(idx)/255, colorProfiles.Blue(idx)/255];

    clear f; 
    
    %======================================================================
    %Outliers
    
    outliers = [];
    % Examine model residuals: boxplot of raw residuals.
    figure(t + length(tractIDs)); k = figure('visible', 'off');
    m = mdl.Residuals.Raw;
    e = eps(max(m(:)));
    boxplot(m)
    % ylabel('Raw Residuals')
    % Suppress figure display.
    set(gcf,'Visible','off');              
    set(0,'DefaultFigureVisible','off');
%?
    % Get indices of the outliers.
    h1 = flipud(findobj(gcf,'tag','Outliers')); % flip order of handles
    for jj = 1 : length( h1 )
        x =  get( h1(jj), 'XData' );
        y =  get( h1(jj), 'YData' );
        for ii = 1 : length( x )
            if not( isnan( x(ii) ) )
                ix = find( abs( m(:,jj)-y(ii) ) < e );
                outliers = cat(1, outliers, ix);
                %                 text( x(ii), y(ii), sprintf( '\\leftarrowY%02d', ix ) )
            end
        end
    end
%?
    outliers = sort(outliers);
%?
    k = gcf; close(k);
    close;
%?
    %k = figure('visible', 'on');
    %set(gcf, 'Visible', 'off');              
    %set(0, 'DefaultFigureVisible', 'off');
%
%?

    %======================================================================

    clf;

    %Remove outliers
    tbl(outliers, :) = []; 

    %recalculate the model

    %generating the model
    mdl = fitlm(tbl, Q);

    %clear the figure
    clf(figure(nfig));

    %startingx, startingy, width height
    f.Position = [1000 1000 800 700];

    %plotting the model
    h = plot(mdl, 'Marker', 'o', 'MarkerEdgeColor', 'white', 'MarkerFaceColor', markerColor, 'MarkerSize', 12);
    %pltLeg = legend('', '', '');
    %set(pltLeg,'visible','off')

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
    set(k, 'EdgeColor', 'none', 'FaceColor', [markerColor(1)*0.55  markerColor(2)*0.55 markerColor(3)*0.55], 'FaceAlpha', '0.2')

    %style the trendline
    set(fitHandle, 'Color', [markerColor(1) markerColor(2) markerColor(3)], 'LineWidth', 3)

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
    xlabel('fa');
    %%%%%%Not a variable
    ylabel('Math Score');
    
    n = size(Math, 1);
    %text(-2.8, 2.2, ['rmse = ' num2str(mdl.RMSE)])
    %text(-2.8, 2.0, ['beta = ' num2str(mdl.Coefficients{2,1})])
    %text(-2.8, 1.8, ['p = ' num2str(mdl.Coefficients.pValue(2))])
    %text(-2.8, 1.6, ['n = ' num2str(n)]);
    
    rsqSimpleLinMath.beta(t) = mdl.Coefficients{2,1}; 
    rsqSimpleLinMath.n(t) = n; 
    rsqSimpleLinMath.p_value(t) = mdl.Coefficients.pValue(2); 
    rsqSimpleLinMath.rmse(t) = mdl.RMSE; 

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
    set(gca, 'XLim', [0.3 0.5], 'XTick', [0.3 0.4 0.5]);
    xax.FontName = fontname;
    xax.FontSize = fontsize;

    % yaxis
    yax = get(gca,'yaxis');
    yax.TickDirection = 'out';
    yax.TickLength = [yticklength yticklength];
    set(gca, 'YLim', [-5 10], 'YTick', [-5 0 5 10]);
    yax.FontName = fontname;
    yax.FontSize = fontsize;
    yax.FontAngle = fontangle;

    %change figure background to white
    set(gcf, 'color', 'w')

    %===========================================================================
    
    hold off
    
    %adding title and color to the model
    plotTitle = {char(tractIDs(t))};
    plotTitle = strjoin(['Simple Linear Model for', plotTitle]);
    title(plotTitle);
    xlabel(measure);
    ylabel('Math Scores');
    
    %export figure as a png file
    mainpath = '/Volumes/LANDLAB/projects/sfa/supportFiles/final-plots';
    filename = strjoin(['tois_Math_', tractIDs(t)]); % Define a filename for each plot
    fullFilePath = fullfile(mainpath, filename); % Create the full file path
    saveas(nfig, fullFilePath, 'png'); % Save the figure as a PNG file
    
end

%============== Generate Reading Plots ==============
rsqSimpleLinRead = table(tractIDs'); 

nrow = size(rsqSimpleLinRead, 1); 

%Re-establishing rsqSimpleLinRead
rsqSimpleLinRead = table(tractIDs'); 
rsqSimpleLinRead.beta = zeros(nrow, 1); 
rsqSimpleLinRead.n = zeros(nrow, 1); 
rsqSimpleLinRead.p_value = zeros(nrow, 1); 
rsqSimpleLinRead.rmse = zeros(nrow, 1); 

for t = 1:length(tractIDs)
    
    found = false; 
    
    while all((strcmp(d.Properties.VariableNames, tractIDs{t})) == 0)
        t = t + 1; 
        if (t > length(tractIDs))
            found = true;  
            break; 
        end
    end
    
    if found
        break;
    end
    
    
    nfig = nfig + 1;
    f = figure(nfig); hold on;

    %startingx, startingy, width height
    f.Position = [1000 1000 800 700];
    
    hold on 
    
    %==================== Reading ====================
    %plotting a nonlinear aggression model 
    ReadTBL = table(d.subIDs, read, 'VariableNames', {'Subject', 'Data'});  
    ReadTBL = ReadTBL(~any(ismissing(ReadTBL), 2), :);
    ReadTBL = sortrows(ReadTBL); 
    
    xVarTBL = table(d.subIDs, d.(char(tractIDs(t))), 'VariableNames', {'Subject', 'Data'});
    xVarTBL(xVarTBL.Data < 0.3, :)= [];
    xVarTBL = xVarTBL(~any(ismissing(xVarTBL), 2), :);
    xVarTBL = sortrows(xVarTBL); 
    
    sexTBL = table(d.subIDs, d.Sex, 'VariableNames', {'Subject', 'Data'});
    sexTBL = sexTBL(~any(ismissing(sexTBL), 2), :);
    sexTBL = sortrows(sexTBL); 
    
    %Merge Tables
    tblMerge = innerjoin(sexTBL, xVarTBL, 'Keys', 'Subject');
    tbl = innerjoin(tblMerge, ReadTBL, 'Keys', 'Subject');
    
    %Rename Headers
    tbl.Properties.VariableNames{'Data_sexTBL'} = 'Sex';
    tbl.Properties.VariableNames{'Data_xVarTBL'} = 'xVar';
    tbl.Properties.VariableNames{'Data'} = 'Read';
    
    tbl(any(ismissing(tbl), 2), :) = [];
    
    Read = tbl.Read; 
    xVar = tbl.xVar; 
    Sex = tbl.Sex; 
    
    %Check correlation between Sex and Math
    P = 'Read ~ Sex'; 
    mdl2 = fitlm(tbl, P);
      
    %z-score xVar and Math    
    xVar = zscore(xVar, [], 'omitnan');
    Read = zscore(Read, [], 'omitnan');

    %%defining the line to fit the model to
    Q = 'Read ~ xVar';

    %generating the model
    mdl = fitlm(tbl, Q);
    
    %get appropriate RGB color for tract by indexing into colorProfiles.csv
    idx = find(strcmp(colorProfiles.NameOfTrack, char(tractIDs(t))) == 1);
    markerColor = [colorProfiles.Red(idx)/255, colorProfiles.Green(idx)/255, colorProfiles.Blue(idx)/255];

    clear f; 

    %======================================================================
    %Outliers
    
    outliers = [];
    % Examine model residuals: boxplot of raw residuals.
    figure(t + length(tractIDs)); k = figure('visible', 'off');
    m = mdl.Residuals.Raw;
    e = eps(max(m(:)));
    boxplot(m)
    % ylabel('Raw Residuals')
    % Suppress figure display.
    set(gcf,'Visible','off');              
    set(0,'DefaultFigureVisible','off');
%?
    % Get indices of the outliers.
    h1 = flipud(findobj(gcf,'tag','Outliers')); % flip order of handles
    for jj = 1 : length( h1 )
        x =  get( h1(jj), 'XData' );
        y =  get( h1(jj), 'YData' );
        for ii = 1 : length( x )
            if not( isnan( x(ii) ) )
                ix = find( abs( m(:,jj)-y(ii) ) < e );
                outliers = cat(1, outliers, ix);
                %                 text( x(ii), y(ii), sprintf( '\\leftarrowY%02d', ix ) )
            end
        end
    end
%?
    
     outliers = sort(outliers);
%?
    k = gcf; close(k);
%?
    k = figure('visible', 'on');
    set(gcf, 'Visible', 'off');              
    set(0, 'DefaultFigureVisible', 'off');
    
    close;
%
%

    %======================================================================

    clf;

    %Remove outliers
    tbl(outliers, :) = []; 

    %recalculate the model

    %generating the model
    mdl = fitlm(tbl, Q);

    %clear the figure
    clf(figure(nfig));

    %startingx, startingy, width height
    f.Position = [1000 1000 800 700];

    %plotting the model
    h = plot(mdl, 'Marker', 'o', 'MarkerEdgeColor', 'white', 'MarkerFaceColor', markerColor, 'MarkerSize', 12);
    %pltLeg = legend('', '', '');
    %set(pltLeg,'visible','off')

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
    set(k, 'EdgeColor', 'none', 'FaceColor', [markerColor(1)*0.55  markerColor(2)*0.55 markerColor(3)*0.55], 'FaceAlpha', '0.2')

    %style the trendline
    set(fitHandle, 'Color', [markerColor(1) markerColor(2) markerColor(3)], 'LineWidth', 3)

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
    xlabel('fa');
    %%%%%%Not a variable
    ylabel('Read Score');
    
    n = size(Read, 1);
    %text(-2.8, 1.6, ['n = ' num2str(n)]);
    %text(-2.8, 2.2, ['rmse = ' num2str(mdl.RMSE)])
    %text(-2.8, 2.0, ['beta = ' num2str(mdl.Coefficients{2,1})])
    %text(-2.8, 1.8, ['p = ' num2str(mdl.Coefficients.pValue(2))])

    %set scale of y-axis
    %ylim([0.3 0.6])

    %add adjusted r squared to table.
    rsqSimpleLinRead.beta(t) = mdl.Coefficients{2,1}; 
    rsqSimpleLinRead.n(t) = n; 
    rsqSimpleLinRead.p_value(t) = mdl.Coefficients.pValue(2); 
    rsqSimpleLinRead.rmse(t) = mdl.RMSE; 

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
    set(gca, 'XLim', [0.3 0.5], 'XTick', [0.3 0.4 0.5]);
    xax.FontName = fontname;
    xax.FontSize = fontsize;

    % yaxis
    yax = get(gca,'yaxis');
    yax.TickDirection = 'out';
    yax.TickLength = [yticklength yticklength];
    set(gca, 'YLim', [-5 20], 'YTick', [-5 0 5 10 15 20]);
    yax.FontName = fontname;
    yax.FontSize = fontsize;
    yax.FontAngle = fontangle;

    %change figure background to white
    set(gcf, 'color', 'w')

    %===========================================================================
    
    hold off
    
    %adding title and color to the model
    plotTitle = {char(tractIDs(t))};
    plotTitle = strjoin(['Simple Linear Model for', plotTitle]);
    title(plotTitle);
    xlabel(measure);
    ylabel('Read Scores');
    
    %export figure as a png file
    mainpath = '/Volumes/LANDLAB/projects/sfa/supportFiles/final-plots';
    filename = strjoin(['tois_Read_', tractIDs(t)]); % Define a filename for each plot
    fullFilePath = fullfile(mainpath, filename); % Create the full file path
    saveas(nfig, fullFilePath, 'png'); % Save the figure as a PNG file
    
end

%============== Export rsqTable as a csv ==============
%local path to save table: 
%mainpath = '/Users/land/Desktop/projectTrackProfiles/supportFiles';

%table_path_format_rsqTAdj = fullfile(mainpath, 'rsqTableAdj.csv');
%table_path_format_rsqTOrd = fullfile(mainpath, 'rsqTableOrd.csv');
%table_path_format_aicTable = fullfile(mainpath, 'aicTable.csv');

%funally, save tables
%writetable(rsqTableAdj, table_path_format_rsqTAdj);
%writetable(rsqTableOrd, table_path_format_rsqTOrd);
%writetable(aicTable, table_path_format_aicTable);