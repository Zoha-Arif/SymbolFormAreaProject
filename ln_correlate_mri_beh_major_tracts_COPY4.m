%clear all; close all; clc;

%% ====================================================================== %%
% White matter measure
wmmeasure = 'fa';

% Hemisphere
hemisphere = 'left'; %left, right, both

% Set working directories.
%w_measures = {'fa', 'md'}; % 'md', 'ad', 'md', 'rd', 'ndi', 'isovf', 'odi', 'map', 'T1', 'R1';
rootdir = '/Volumes/LANDLAB/projects/sfa/';
blprojectid = 'proj-63e6b75dc538c16a8225593e-majortracts'; 
mp2rage = 'no';

%generate column of tracts of interest ids
tractIDs1 = ["vof" "parc" "tpc" "mdlfang" "mdlfspl" "ilf" "ifof" "slf1and2" "slf3" "fat"];
tractIDs2 = ["VOF" "pArc" "TPC" "MDLFang" "MDLFspl" "ILF" "IFOF" "SLF1And2" "SLF3" "Aslant"];
nfig = 30; 
colorProfiles = '/Volumes/LANDLAB/projects/sfa/supportFiles/colorProfiles.csv';
colorProfiles = readtable(colorProfiles);

%% Correlate fa with math and reading scores

%insert local path of Tshort.csv and Tlong.csv file
d = readtable(fullfile(rootdir, 'supportfiles', ['Tshort' '.csv']));
Tlong = readtable(fullfile(rootdir, 'supportfiles', ['Tlong' '.csv']));

%% ====================================================================== %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Behavioral Data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

math = d.T2_WJIV_AppProb_raw - d.T1_WJIV_AppProb_raw;
read = d.T2_WJIV_LWID_raw - d.T1_WJIV_LWID_raw;
%% ====================================================================== %%

%============== Generate Math Plots ==============

%generate column of tracts of interest ids
tractIDs = {'VOF', 'pArc', 'TPC', 'MDLFang', 'MDLFspl', ...
    'ILF', 'IFOF', 'SLF1And2', 'SLF3', 'Aslant'}; 

%Re-establishing rsqSimpleLinMath
rsqSimpleLinMath = table(tractIDs'); 
nrow = size(rsqSimpleLinMath, 1); 
rsqSimpleLinMath.beta = zeros(nrow, 1); 
rsqSimpleLinMath.n = zeros(nrow, 1); 
rsqSimpleLinMath.p_value = zeros(nrow, 1); 
rsqSimpleLinMath.rmse = zeros(nrow, 1); 

%Re-establishing rsqSimpleLinMathSex
rsqSimpleLinMSex = table(tractIDs'); 
rsqSimpleLinMSex.beta = zeros(nrow, 1); 
rsqSimpleLinMSex.n = zeros(nrow, 1); 
rsqSimpleLinMSex.p_value = zeros(nrow, 1); 
rsqSimpleLinMSex.rmse = zeros(nrow, 1); 

for t = 1:length(tractIDs)
    
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
    
    xVarTBL = table(d.subIDs, d.(char(strcat(hemisphere, tractIDs(t)))), 'VariableNames', {'Subject', 'Data'});
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
      
    %z-score xVar and Math    
    xVar = zscore(xVar, [], 'omitnan');
    Math = zscore(Math, [], 'omitnan');

    %Check correlation between Sex and Math
    P = 'Math ~ Sex'; 
    mdlMSex = fitlm(tbl, P);
    
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

    %clear figure
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
    %%%%%%%To print on plots
    %text(-2.8, 2.2, ['rmse = ' num2str(mdl.RMSE)])
    %text(-2.8, 2.0, ['beta = ' num2str(mdl.Coefficients{2,1})])
    %text(-2.8, 1.8, ['p = ' num2str(mdl.Coefficients.pValue(2))])
    %text(-2.8, 1.6, ['n = ' num2str(n)]);
    
    rsqSimpleLinMath.beta(t) = mdl.Coefficients{2,1}; 
    rsqSimpleLinMath.n(t) = n; 
    rsqSimpleLinMath.p_value(t) = mdl.Coefficients.pValue(2); 
    rsqSimpleLinMath.rmse(t) = mdl.RMSE; 
    
    rsqSimpleLinMSex.beta(t) = mdlMSex.Coefficients{2,1}; 
    rsqSimpleLinMSex.n(t) = n; 
    rsqSimpleLinMSex.p_value(t) = mdlMSex.Coefficients.pValue(2); 
    rsqSimpleLinMSex.rmse(t) = mdlMSex.RMSE; 

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
    filename = strjoin(['major_', hemisphere, '_Math_', tractIDs(t)]); % Define a filename for each plot
    fullFilePath = fullfile(mainpath, filename); % Create the full file path
    saveas(nfig, fullFilePath, 'png'); % Save the figure as a PNG file
    
end

%============== Generate Reading Plots ==============

%Re-establishing rsqSimpleLinRead
rsqSimpleLinRead = table(tractIDs'); 
nrow = size(rsqSimpleLinRead, 1); 
rsqSimpleLinRead.beta = zeros(nrow, 1); 
rsqSimpleLinRead.n = zeros(nrow, 1); 
rsqSimpleLinRead.p_value = zeros(nrow, 1); 
rsqSimpleLinRead.rmse = zeros(nrow, 1); 

%Re-establishing rsqSimpleLinMathSex
rsqSimpleLinRSex = table(tractIDs'); 
rsqSimpleLinRSex.beta = zeros(nrow, 1); 
rsqSimpleLinRSex.n = zeros(nrow, 1); 
rsqSimpleLinRSex.p_value = zeros(nrow, 1); 
rsqSimpleLinRSex.rmse = zeros(nrow, 1); 

for t = 1:length(tractIDs)
    
    nfig = nfig + 1;
    f = figure(nfig);

    %startingx, startingy, width height
    f.Position = [1000 1000 800 700];
    
    hold on 
    
    %==================== Reading ====================
    %plotting a nonlinear aggression model 
    ReadTBL = table(d.subIDs, read, 'VariableNames', {'Subject', 'Data'});  
    ReadTBL = ReadTBL(~any(ismissing(ReadTBL), 2), :);
    ReadTBL = sortrows(ReadTBL); 
    
    xVarTBL = table(d.subIDs, d.(char(strcat(hemisphere, tractIDs(t)))), 'VariableNames', {'Subject', 'Data'});
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
      
    %z-score xVar and Math    
    xVar = zscore(xVar, [], 'omitnan');
    Read = zscore(Read, [], 'omitnan');
    
    %Check correlation between Sex and Math
    P = 'Read ~ Sex'; 
    mdlRSex = fitlm(tbl, P);

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
    
    rsqSimpleLinRSex.beta(t) = mdlRSex.Coefficients{2,1}; 
    rsqSimpleLinRSex.n(t) = n; 
    rsqSimpleLinRSex.p_value(t) = mdlRSex.Coefficients.pValue(2); 
    rsqSimpleLinRSex.rmse(t) = mdlRSex.RMSE; 

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
    filename = strjoin(['major_', hemisphere, '_Read_', tractIDs(t)]); % Define a filename for each plot
    fullFilePath = fullfile(mainpath, filename); % Create the full file path
    saveas(nfig, fullFilePath, 'png'); % Save the figure as a PNG file
    
end

%%%%%%%%%%%%%%MULTIPLE VARIABLE MODEL -- READING %%%%%%%%%%%%%%%%%%%%%
%%defining the line to fit the model to
 ReadTBL = table(d.subIDs, read, 'VariableNames', {'Subject', 'Data'});  
 ReadTBL = ReadTBL(~any(ismissing(ReadTBL), 2), :);
 ReadTBL = sortrows(ReadTBL); 
 
 %pArc
 pArcTBL = table(d.subIDs, d.(char(strcat(hemisphere, tractIDs(2)))), 'VariableNames', {'Subject', 'Data'});
 pArcTBL = pArcTBL(~any(ismissing(pArcTBL), 2), :);
 pArcTBL = sortrows(pArcTBL); 
 mask = ismember(ReadTBL.Subject, pArcTBL.Subject);
 ReadTBL = ReadTBL(mask, :); 
 mask2 = ismember(pArcTBL.Subject, ReadTBL.Subject);
 pArcTBL = pArcTBL(mask2, :); 
 
 %tpc
 tpcTBL = table(d.subIDs, d.(char(strcat(hemisphere, tractIDs(3)))), 'VariableNames', {'Subject', 'Data'});
 tpcTBL = tpcTBL(~any(ismissing(tpcTBL), 2), :);
 tpcTBL = sortrows(tpcTBL); 
 mask = ismember(ReadTBL.Subject, tpcTBL.Subject);
 ReadTBL = ReadTBL(mask, :); 
 mask2 = ismember(tpcTBL.Subject, ReadTBL.Subject);
 tpcTBL = tpcTBL(mask2, :); 
 
 %MDLfang
 MDLfangTBL = table(d.subIDs, d.(char(strcat(hemisphere, tractIDs(4)))), 'VariableNames', {'Subject', 'Data'});
 MDLfangTBL = MDLfangTBL(~any(ismissing(MDLfangTBL), 2), :);
 MDLfangTBL = sortrows(MDLfangTBL); 
 mask = ismember(ReadTBL.Subject, MDLfangTBL.Subject);
 ReadTBL = ReadTBL(mask, :); 
 mask2 = ismember(MDLfangTBL.Subject, ReadTBL.Subject);
 MDLfangTBL = MDLfangTBL(mask2, :); 
 
 %MDLFspl
 MDLFsplTBL = table(d.subIDs, d.(char(strcat(hemisphere, tractIDs(5)))), 'VariableNames', {'Subject', 'Data'});
 MDLFsplTBL = MDLFsplTBL(~any(ismissing(MDLFsplTBL), 2), :);
 MDLFsplTBL = sortrows(MDLFsplTBL); 
 mask = ismember(ReadTBL.Subject, MDLFsplTBL.Subject);
 ReadTBL = ReadTBL(mask, :); 
 mask2 = ismember(MDLFsplTBL.Subject, ReadTBL.Subject);
 MDLFsplTBL = MDLFsplTBL(mask2, :); 
 
 tbl = table(ReadTBL.Data, pArcTBL.Data, tpcTBL.Data, MDLfangTBL.Data, MDLFsplTBL.Data, 'VariableNames', {'Read', 'pArc', 'TPC', 'MDLfang', 'MDLFspl'}); 
 tbl(any(ismissing(tbl), 2), :) = [];

 Read = tbl.Read; 
 pArc = tbl.pArc; 
 TPC = tbl.TPC; 
 MDLfang = tbl.MDLfang; 
 MDLFspl = tbl.MDLFspl; 
    
 pArc = zscore(pArc, [], 'omitnan');
 TPC = zscore(TPC, [], 'omitnan');
 MDLfang = zscore(MDLfang, [], 'omitnan');
 MDLFspl = zscore(MDLFspl, [], 'omitnan');
 Read = zscore(Read, [], 'omitnan');


Q = 'Read ~ pArc + TPC + MDLfang + MDLFspl';

%generating the model
mdlRead = fitlm(tbl, Q);

%%%%%%%%%%%%%%MULTIPLE VARIABLE MODEL -- MATH %%%%%%%%%%%%%%%%%%%%%
%%defining the line to fit the model to
 MathTBL = table(d.subIDs, math, 'VariableNames', {'Subject', 'Data'});  
 MathTBL = MathTBL(~any(ismissing(MathTBL), 2), :);
 MathTBL = sortrows(MathTBL); 
 
 %pArc
 pArcTBL = table(d.subIDs, d.(char(strcat(hemisphere, tractIDs(2)))), 'VariableNames', {'Subject', 'Data'});
 pArcTBL = pArcTBL(~any(ismissing(pArcTBL), 2), :);
 pArcTBL = sortrows(pArcTBL); 
 mask = ismember(MathTBL.Subject, pArcTBL.Subject);
 MathTBL = MathTBL(mask, :); 
 mask2 = ismember(pArcTBL.Subject, MathTBL.Subject);
 pArcTBL = pArcTBL(mask2, :); 
 
 %tpc
 tpcTBL = table(d.subIDs, d.(char(strcat(hemisphere, tractIDs(3)))), 'VariableNames', {'Subject', 'Data'});
 tpcTBL = tpcTBL(~any(ismissing(tpcTBL), 2), :);
 tpcTBL = sortrows(tpcTBL); 
 mask = ismember(MathTBL.Subject, tpcTBL.Subject);
 MathTBL = MathTBL(mask, :); 
 mask2 = ismember(tpcTBL.Subject, MathTBL.Subject);
 tpcTBL = tpcTBL(mask2, :); 
 
 %MDLfang
 MDLfangTBL = table(d.subIDs, d.(char(strcat(hemisphere, tractIDs(4)))), 'VariableNames', {'Subject', 'Data'});
 MDLfangTBL = MDLfangTBL(~any(ismissing(MDLfangTBL), 2), :);
 MDLfangTBL = sortrows(MDLfangTBL); 
 mask = ismember(MathTBL.Subject, MDLfangTBL.Subject);
 MathTBL = MathTBL(mask, :); 
 mask2 = ismember(MDLfangTBL.Subject, MathTBL.Subject);
 MDLfangTBL = MDLfangTBL(mask2, :); 
 
 %MDLFspl
 MDLFsplTBL = table(d.subIDs, d.(char(strcat(hemisphere, tractIDs(5)))), 'VariableNames', {'Subject', 'Data'});
 MDLFsplTBL = MDLFsplTBL(~any(ismissing(MDLFsplTBL), 2), :);
 MDLFsplTBL = sortrows(MDLFsplTBL); 
 mask = ismember(MathTBL.Subject, MDLFsplTBL.Subject);
 MathTBL = MathTBL(mask, :); 
 mask2 = ismember(MDLFsplTBL.Subject, MathTBL.Subject);
 MDLFsplTBL = MDLFsplTBL(mask2, :); 
 
 tbl = table(MathTBL.Data, pArcTBL.Data, tpcTBL.Data, MDLfangTBL.Data, MDLFsplTBL.Data, 'VariableNames', {'Math', 'pArc', 'TPC', 'MDLfang', 'MDLFspl'}); 
 tbl(any(ismissing(tbl), 2), :) = [];

 Math = tbl.Math; 
 pArc = tbl.pArc; 
 TPC = tbl.TPC; 
 MDLfang = tbl.MDLfang; 
 MDLFspl = tbl.MDLFspl; 
    
 pArc = zscore(pArc, [], 'omitnan');
 TPC = zscore(TPC, [], 'omitnan');
 MDLfang = zscore(MDLfang, [], 'omitnan');
 MDLFspl = zscore(MDLFspl, [], 'omitnan');
 Math = zscore(Math, [], 'omitnan');


Q = 'Math ~ pArc + TPC + MDLfang + MDLFspl';

%generating the model
mdlMath = fitlm(tbl, Q);

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