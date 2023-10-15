function runFlexbusAnalysis(data_repo)


% get list of all folders we need to iterate through
items = dir(data_repo);
for i = 1:numel(items)
  % Check if the item is a directory and not '.' or '..'
  if items(i).isdir && ~ismember(items(i).name, {'.', '..'})
    % load V_DF and cus_experience files from this folder

    thisFolder = [data_repo, items(i).name, filesep];
    folderItems = dir(thisFolder);
    allFiles = string({folderItems.name});

    V_DF_file = fullfile(thisFolder, allFiles(find(contains(allFiles,'V_DF'))));
    data = readtable(V_DF_file);  % You can use readmatrix for newer MATLAB versions



    cus_file = fullfile(thisFolder, allFiles(find(contains(allFiles,'cus_'))));
    cus_data = readtable(cus_file);  % You can use readmatrix for newer MATLAB versions

  end
end

% customers are in the .txt file NOT in numerical order
% customers appear IN SAME ORDER in V_DF.csv
% cus number is in col 1 of V_DF
% x,y, coordinates are in cols 3,4 of V_DF if want to check

% in cus_experience we have customers in same order using numbers from V_DF
% cus	pickup_bus_node	drop_off_transit_node	ride_time	cus_walking_time
% 734	      704	          728	              4.616379326	6.108753383

