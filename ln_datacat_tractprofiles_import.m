%redirecting folder to correct path. clear.
clear all; clc;

%(!) % define measure (!)
measure = 'fa';

%fetching all subject's subfolders and storing them in subfolders
%insert local path to track profiles here by changing value of mainpath.
mainpath = '/Volumes/LANDLAB/projects/sfa/supportfiles/proj-63e6b75dc538c16a8225593e-majortracts'; 
cd '/Volumes/LANDLAB/projects/sfa/supportfiles/proj-63e6b75dc538c16a8225593e-majortracts'; 
topLevelFolder = dir(mainpath);

subfolders = topLevelFolder([topLevelFolder(:).isdir]); %find names of all subject's subfolders

% (!) Uncomment line if running script on data for the first time(!)
% (!) Removes the unnecessary '.' and '..' files (!)
% subfolders = subfolders(arrayfun(@(x) x.name(1), subfolders) ~= '.');

% Keep only names that are subject folders.
subfolders = subfolders(arrayfun(@(x) x.name(1), subfolders) == '2');

%tables is a cell array to store each sub's table in. We can use explicit loops or cellfun to
%apply the same code to each table. 
tables = {}; 
Headers = {'subjectID','structureID', 'nodeID', 'fa', 'x_coords', 'y_coords', 'z_coords'};
Tlong = cell2table(cell(0,7),'VariableNames', Headers);

%insert local path to csv here by changing value of mainpath.
% Load beh data
data_beh = '/Volumes/LANDLAB/projects/sfa/supportFiles/ln_data_beh_n68_20231031 copy.csv'; 
data_beh = readtable(data_beh);

% MISSING brainlife: 2029, 2126
% Identify outliers for removal.
beh_outlier = []; % 
missing_mridata = []; % beh data, but no mri data
missing_tracts = []; % tractography failures
missing_run2 = []; % technical issues during scan resulted in no run-2
bad_testretest_fa = [2009 2069 2077]; % test-retest analysis on FA revealed diffusion data to be unreliable
bad_testretest_snr = [2167]; % test-retest analysis on SNR revealed diffusion data to be unreliable, note that SNR for 2001 and 2009 was also outrageously high, 2029 was also very high
brain_anomaly = []; % bain anomaly 
remove = cat(2, beh_outlier, missing_mridata, missing_tracts, brain_anomaly, missing_run2, ...
    bad_testretest_fa, bad_testretest_snr);

idxRemove = ismember(data_beh.subID, remove); 
data_beh(idxRemove, :) = [];

beh = sortrows(data_beh); clear data_beh;

%============== Generate Tlong ==============

for i = 1:(length(subfolders)) 
    %get subpath to each subject's folder
    subpath = fullfile(mainpath, subfolders(i).name); %pathToSubfolder
    
    %cd into subject's folder
    cd(subpath); 
    
    %find path of folder containing .csv file
    subpathCSV = strcat(subpath, '/run-1'); 
    
    
    %cd into file containing .csv file
    if exist(subpathCSV, 'dir')
        cd(subpathCSV);
    
        %grab .csv file and convert file path to string
        csvfile = dir("tractmeasures.csv");
        if length(csvfile) == 1
            csvFullPaths = fullfile({csvfile.folder}, {csvfile.name});
            csvFinalPath = string(csvFullPaths(1:1));
        else
            csvfile = dir("output_FiberStats.csv");
            csvFullPaths = fullfile({csvfile.folder}, {csvfile.name});
            csvFinalPath = string(csvFullPaths(1:1));
        end
    
        %read table into variable T
        T = readtable(csvFinalPath);
        T = T((T.nodeID >= 20 & T.nodeID <= 180), :);

        Tnew = table; 
        Tnew.subjectID = T.subjectID; 
        Tnew.structureID = T.structureID;
        Tnew.nodeID = T.nodeID;
        Tnew.fa = T.fa; 
        Tnew.x_coords = T.x_coords;
        Tnew.y_coords = T.y_coords;
        Tnew.z_coords = T.z_coords;

        Tnew(any(ismissing(Tnew),2), :) = [];
    
        if(~isempty(Tnew.fa))
            Tlong = vertcat(Tnew, Tlong);
        end

        %add table to cell array tables
        %tables{i} = Tnew;

        %cd out of subject's directory for next loop iteration
        cd ..
    end 
    
    cd ..

end

%sort Tlong
Tlong = sortrows(Tlong, 'subjectID'); 
%create new columns in Tlong for Age, Sex, Handedness, and Bin
Tlong.Age = repmat("NA", height(Tlong), 1);
Tlong.Sex = repmat("NA", height(Tlong), 1);
Tlong.Handedness = repmat("NA", height(Tlong), 1);
Tlong.Bins = repmat("NA", height(Tlong), 1);

nrow = size(Tlong, 1);
Tlong.Age = zeros(nrow, 1);    %0
Tlong.Sex = zeros(nrow, 1);    %0
Tlong.Handedness = zeros(nrow, 1);  %0

%fill in Age, Sex, and Handedness data in Tlong by copying from csv. 
idx = ismember(beh.subID, unique(Tlong.subjectID));
T4 = beh(idx, :); 

subIDs = unique(Tlong.subjectID);

clear i; 
for i = 1:length(subIDs)

    idx = find(Tlong.subjectID == subIDs(i));
    idx2 = find(T4.subID == subIDs(i));
    
    if(~isempty(idx2))
        Tlong.Age(idx) = repmat(T4.T1_Age(idx2), [1 length(idx)]);
        Tlong.Sex(idx) = repmat(T4.Sex_0M_1F(idx2), [1 length(idx)]);
        %Tlong.Handedness(idx) = repmat(T4.T1_Handedness_0R1L2A(idx2), [1 length(idx)]);
        
        %Tests at T1
        Tlong.T1_WJIV_AppProb_raw(idx) = repmat(T4.T1_WJIV_AppProb_raw(idx2), [1 length(idx)]);
        Tlong.T1_WJIV_LWID_raw(idx) = repmat(T4.T1_WJIV_LWID_raw(idx2), [1 length(idx)]);
        
        %Tests at T2
        Tlong.T2_WJIV_AppProb_raw(idx) = repmat(T4.T2_WJIV_AppProb_raw(idx2), [1 length(idx)]);
        Tlong.T2_WJIV_LWID_raw(idx) = repmat(T4.T2_WJIV_LWID_raw(idx2), [1 length(idx)]);
        
        %Tests at T3
        %Tlong.T3_WJIV_AppProb_raw(idx) = repmat(T4.T3_WJIV_AppProb_raw(idx2), [1 length(idx)]);
        %Tlong.T3_WJIV_LWID_raw(idx) = repmat(T4.T3_WJIV_LWID_raw(idx2), [1 length(idx)]);
        
    end
        
    clear idx; 
    clear idx2;
end

%============== Cleaning out Tlong ==============

%replace all NaN in sex with 0
columnToReplace = 'Sex';
columnData = Tlong.Sex;
nanIndices = isnan(columnData); 
columnData(nanIndices) = 0; 
Tlong.Sex = columnData; 

%remove all rows that are missing behavioral data (indicated by 0)
columnsToCheck = {'T1_WJIV_AppProb_raw', 'T1_WJIV_LWID_raw', 'T2_WJIV_AppProb_raw'...
    'T2_WJIV_LWID_raw'}; 
zeroRows = any(Tlong{:, columnsToCheck} == 0, 2); 
rowsToRemove = find(zeroRows); 
Tlong(rowsToRemove, :) = []; 

%remove all rows that are missing behavioral data (indicated by NaN)
nanRows = any(ismissing(Tlong{:, columnsToCheck}, 2)); 
Tlong(nanRows, :) = []; 

%remove all rows that are missing behavioral data (indicated by 9999)
nineRows = any(Tlong{:, columnsToCheck} == 9999, 2); 
rowsToRemove = find(nineRows); 
Tlong(rowsToRemove, :) = []; 

%============== Generate Tshort ==============
tractIDs = unique(Tlong.structureID);
subIDs = unique(Tlong.subjectID); 

Tshort = table(subIDs);  

for s = 1:length(subIDs)
    Sidx = find(Tlong.subjectID == subIDs(s));
    
    for t = 1:length(tractIDs)
        Tidx = find(strcmp(Tlong.structureID, char(tractIDs(t))) == 1);
        STidx = intersect(Tidx, Sidx);
        faTract = Tlong(STidx, measure);
        faTract = table2array(faTract);
        Tshort.(char(tractIDs(t)))(s) = repmat(mean(faTract, 1), [1 length(s)]); 
    end
    
    Tshort.Sex(s) = Tlong.Sex(Sidx(1));
     
    Tshort.T1_WJIV_AppProb_raw(s) = Tlong.T1_WJIV_AppProb_raw(Sidx(1));
    Tshort.T1_WJIV_LWID_raw(s) = Tlong.T1_WJIV_LWID_raw(Sidx(1));
    
    Tshort.T2_WJIV_AppProb_raw(s) = Tlong.T2_WJIV_AppProb_raw(Sidx(1));
    Tshort.T2_WJIV_LWID_raw(s) = Tlong.T2_WJIV_LWID_raw(Sidx(1));
    
    %Tshort.T3_WJIV_AppProb_raw(s) = Tlong.T3_WJIV_AppProb_raw(Sidx(1));
    %Tshort.T3_WJIV_LWID_raw(s) = Tlong.T3_WJIV_LWID_raw(Sidx(1));
    
end

%Remove all data from Tshort that has an fa value lower than 0.30


%============== Export Tshort and Tlong as a csv ==============
%local path to save table: 
mainpath = '/Volumes/LANDLAB/projects/sfa/supportfiles'; 

table_path_format_tshort = fullfile(mainpath, 'Tshort.csv');
table_path_format_tlong = fullfile(mainpath, 'Tlong.csv');

%finally, save tables
writetable(Tshort, table_path_format_tshort);
writetable(Tlong, table_path_format_tlong);
