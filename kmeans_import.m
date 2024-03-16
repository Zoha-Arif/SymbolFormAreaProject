%clear all; close all; clc;
%% ====================================================================== %%
%read in data

%fetching all subject's subfolders and storing them in subfolders
%insert local path to track profiles here by changing value of mainpath.
mainpath = '/Volumes/LANDLAB/projects/sfa/supportfiles/2019_rois'; 
cd '/Volumes/LANDLAB/projects/sfa/supportfiles/2019_rois'; 
topLevelFolder = dir(mainpath);

subfolders = topLevelFolder([topLevelFolder(:).isdir]); %find names of all subject's subfolders

% (!) Uncomment line if running script on data for the first time(!)
% (!) Removes the unnecessary '.' and '..' files (!)
% subfolders = subfolders(arrayfun(@(x) x.name(1), subfolders) ~= '.');

% Keep only names that are subject folders.
subfolders = subfolders(arrayfun(@(x) x.name(1), subfolders) == '2');

%Create table to store all of the data in it.
tables = {}; 
Headers = {'subjectID','RNFA', 'LNFA', 'RVWFA', 'LVWFA'};
data = cell2table(cell(0,5),'VariableNames', Headers);
data.subjectID = num2cell(data.subjectID); 
data.RNFA = num2cell(data.RNFA);  
data.LNFA = num2cell(data.LNFA); 
data.RVWFA = num2cell(data.RVWFA); 
data.LVWFA = num2cell(data.LVWFA); 

for i = 1:(length(subfolders)) 
    %get subpath to each subject's folder
    subpath = fullfile(mainpath, subfolders(i).name); %pathToSubfolder
    
    %cd into subject's folder
    cd(subpath); 
    
    %find path of folder containing .csv file
    subpathNII = strcat(subpath, '/1_22_2024.round_2.ventral_temporal_parietal_occipital_rois.no_inflation.run-1/rois'); 
    
    %cd into file containing .csv file
    if (exist(subpathNII, 'dir') && ~isempty(subpathNII))
        
        cd(subpathNII);
    
        %%%%% Right NFA -> parietal %%%%% 
         %grab nii.gz file and read it in
        fileNames = {'track_freesurfer-12125_to_rh.nfa_exact_1mm_RAS_FiberEndpoint.nii.gz'...
            'track_freesurfer-12126_to_rh.nfa_exact_1mm_RAS_FiberEndpoint.nii.gz'...
            'track_freesurfer-12127_to_rh.nfa_exact_1mm_RAS_FiberEndpoint.nii.gz' ...
            'track_freesurfer-12157_to_rh.nfa_exact_1mm_RAS_FiberEndpoint.nii.gz'}; 
        
        niiData = cell(1, numel(fileNames)); 
        
        combinedData = []; 
        
        for j = 1:numel(fileNames)
            
            try
                niiData{j} = niftiread(fileNames{j}); 
                
                %make sure the sizes of the two matrices are the same
                sizeNII = size(niiData{j});
                sizeCombined = size(combinedData); 
                
                if (~isequal(sizeNII, sizeCombined) && (j ~= 1))
                    
                    if (sizeNII < sizeCombined)
                        paddingSize = sizeCombined - sizeNII;
                        niiData{j} = padarray(niiData{j}, paddingSize, 0, 'post');
                    else
                        paddingSize = sizeNII - sizeCombined;
                        combinedData = padarray(combinedData, paddingSize, 0, 'post');
                    end
                    
                end
                
                %combine the data
                if (j == 1)
                    combinedData = niiData{j}; 
                else
                    combinedData = combinedData + niiData{j}; 
                end
                
                
            catch ME
                disp(['Failed to read ', fileNames{j}, ': ', ME.message]); 
            end
        end
        
        %find x, y, and z coordinates
        [row, col, pag] = ind2sub(size(combinedData), find(combinedData));
        
        %store data in table
        data.RNFA{i} = [row, col, pag];
        data.subjectID{i} = str2double(subfolders(i).name); 
       
        clear niiData; clear combinedData; 
        clear row; clear col; clear pag; clear j; 
        
        %%%%% Left NFA -> parietal %%%%% 
         %grab nii.gz file and read it in
        fileNames = {'track_freesurfer-11125_to_lh.nfa_exact_1mm_RAS_FiberEndpoint.nii.gz'...
            'track_freesurfer-11126_to_lh.nfa_exact_1mm_RAS_FiberEndpoint.nii.gz'...
            'track_freesurfer-11127_to_lh.nfa_exact_1mm_RAS_FiberEndpoint.nii.gz' ...
            'track_freesurfer-11157_to_lh.nfa_exact_1mm_RAS_FiberEndpoint.nii.gz'}; 
        
        niiData = cell(1, numel(fileNames)); 
        
        combinedData = []; 
        
        for j = 1:numel(fileNames)
            
            try
                niiData{j} = niftiread(fileNames{j}); 
                
                %make sure the sizes of the two matrices are the same
                sizeNII = size(niiData{j});
                sizeCombined = size(combinedData); 
                
                if (~isequal(sizeNII, sizeCombined) && (j ~= 1))
                    
                    if (sizeNII < sizeCombined)
                        paddingSize = sizeCombined - sizeNII;
                        niiData{j} = padarray(niiData{j}, paddingSize, 0, 'post');
                    else
                        paddingSize = sizeNII - sizeCombined;
                        combinedData = padarray(combinedData, paddingSize, 0, 'post');
                    end
                    
                end
                
                %combine the data
                if (j == 1)
                    combinedData = niiData{j}; 
                else
                    combinedData = combinedData + niiData{j}; 
                end
                
                
            catch ME
                disp(['Failed to read ', fileNames{j}, ': ', ME.message]); 
            end
        end
        
        %find x, y, and z coordinates
        [row, col, pag] = ind2sub(size(combinedData), find(combinedData));
        
        %store data in table
        data.LNFA{i} = [row, col, pag];  
       
        clear niiData; clear combinedData; 
        clear row; clear col; clear pag; clear j; 
        
        %%%%% Right VWFA -> parietal %%%%% 
         %grab nii.gz file and read it in
        fileNames = {'track_freesurfer-12125_to_rh.vwfa_exact_1mm_RAS_FiberEndpoint.nii.gz'...
            'track_freesurfer-12126_to_rh.vwfa_exact_1mm_RAS_FiberEndpoint.nii.gz'...
            'track_freesurfer-12127_to_rh.vwfa_exact_1mm_RAS_FiberEndpoint.nii.gz' ...
            'track_freesurfer-12157_to_rh.vwfa_exact_1mm_RAS_FiberEndpoint.nii.gz'}; 
        
        niiData = cell(1, numel(fileNames)); 
        
        combinedData = []; 
        
        for j = 1:numel(fileNames)
            
            try
                niiData{j} = niftiread(fileNames{j}); 
                
                %make sure the sizes of the two matrices are the same
                sizeNII = size(niiData{j});
                sizeCombined = size(combinedData); 
                
                if (~isequal(sizeNII, sizeCombined) && (j ~= 1))
                    
                    if (sizeNII < sizeCombined)
                        paddingSize = sizeCombined - sizeNII;
                        niiData{j} = padarray(niiData{j}, paddingSize, 0, 'post');
                    else
                        paddingSize = sizeNII - sizeCombined;
                        combinedData = padarray(combinedData, paddingSize, 0, 'post');
                    end
                    
                end
                
                %combine the data
                if (j == 1)
                    combinedData = niiData{j}; 
                else
                    combinedData = combinedData + niiData{j}; 
                end
                
                
            catch ME
                disp(['Failed to read ', fileNames{j}, ': ', ME.message]); 
            end
        end
        
        %find x, y, and z coordinates
        [row, col, pag] = ind2sub(size(combinedData), find(combinedData));
        
        %store data in table
        data.RVWFA{i} = [row, col, pag];  
       
        clear niiData; clear combinedData; 
        clear row; clear col; clear pag; clear j; 
        
        %%%%% Left VWFA -> parietal %%%%% 
         %grab nii.gz file and read it in
        fileNames = {'track_freesurfer-11125_to_lh.vwfa_exact_1mm_RAS_FiberEndpoint.nii.gz'...
            'track_freesurfer-11126_to_lh.vwfa_exact_1mm_RAS_FiberEndpoint.nii.gz'...
            'track_freesurfer-11127_to_lh.vwfa_exact_1mm_RAS_FiberEndpoint.nii.gz' ...
            'track_freesurfer-11157_to_lh.vwfa_exact_1mm_RAS_FiberEndpoint.nii.gz'}; 
        
        niiData = cell(1, numel(fileNames)); 
        
        combinedData = []; 
        
        for j = 1:numel(fileNames)
            
            try
                niiData{j} = niftiread(fileNames{j}); 
                
                %make sure the sizes of the two matrices are the same
                sizeNII = size(niiData{j});
                sizeCombined = size(combinedData); 
                
                if (~isequal(sizeNII, sizeCombined) && (j ~= 1))
                    
                    if (sizeNII < sizeCombined)
                        paddingSize = sizeCombined - sizeNII;
                        niiData{j} = padarray(niiData{j}, paddingSize, 0, 'post');
                    else
                        paddingSize = sizeNII - sizeCombined;
                        combinedData = padarray(combinedData, paddingSize, 0, 'post');
                    end
                    
                end
                
                %combine the data
                if (j == 1)
                    combinedData = niiData{j}; 
                else
                    combinedData = combinedData + niiData{j}; 
                end
                
                
            catch ME
                disp(['Failed to read ', fileNames{j}, ': ', ME.message]); 
            end
        end
        
        %find x, y, and z coordinates
        [row, col, pag] = ind2sub(size(combinedData), find(combinedData));
        
        %store data in table
        data.LVWFA{i} = [row, col, pag];  
       
        clear niiData; clear combinedData; 
        clear row; clear col; clear pag; clear j; 
  
        %cd out of subject's directory for next loop iteration
        cd ..
        cd ..
    end 
    
    cd ..
    
end

%niiInfo = niftiinfo(niiFilePath);

% Export the table to a CSV file
mainpath = '/Volumes/LANDLAB/projects/sfa/supportFiles/';
filename = 'xyz_coords_tois.csv';
fullfilePath = fullfile(mainpath, filename); % Create the full file path
writetable(data, fullfilePath);

