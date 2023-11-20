% data should be in ONE master folder that contains N subfolders
% each with a different 'run' of some sort
% This script will iterate through EACH results subfolder, extract data for
% that 'run' and generate summary statistics.

% sort out folder and file names
% repo_flexbus = [get_repo_folder, 'Flexbus3_v0.8.4\'];
% data_to_load = [repo_flexbus, 'results\grid_results_v2\'];
% data_to_load = 'Q:\REPOS\UrbanMorph\data\';
data_to_load = [get_repo_folder, 'UrbanMorph\data_bicentric\'];

% % ****** REMOVE NOT VIABLE RUNS ********
% badDataDir = 'Q:\REPOS\UrbanMorph\failed_data\';
% items = dir(data_to_load);
% % Check if the item is a directory and not '.' or '..'
% items = items([items.isdir]); items = items(~ismember({items.name}, {'.', '..'}));
% nInstance = numel(items); % this is the instance folder
% for i = 1:nInstance % this is the instance folder
%   thisFolder = [data_to_load, items(i).name, filesep];
%   folderItems = dir(thisFolder); allFiles = string({folderItems.name});
%   cus_file = allFiles(contains(allFiles,'cus_'));
%   route_files =  allFiles(startsWith(allFiles,'route_detail_'));
%   if isempty(cus_file) || isempty(route_files) % we have some results in this instance
%     % move this folder to failed_instances
%     movefile(thisFolder,badDataDir);
%   end
% end

items = dir(data_to_load);
% Check if the item is a directory and not '.' or '..'
items = items([items.isdir]); items = items(~ismember({items.name}, {'.', '..'}));
nFolders = numel(items); % this is the instance folder

nZ = zeros(nFolders,1);
R_instance = table(nZ,nZ,nZ,nZ,nZ,nZ,nZ,nZ,nZ,nZ,nZ,...
  'VariableNames', {'id','nPax','Pax_maxRadius', 'BS_separation',...
  'maxWalkingDist', 'FleetSize', 'VehKms', 'VehEmptyKms', 'CusDirectKms', 'CusTravelledKms','sSeparation'});
R_bus = table([],[],[],[],[],[],[],[],[],[], 'VariableNames', {'id','nPax','Pax_maxRadius', 'BS_separation','maxWalkingDist',...
  'VehKms','VehEmptyKms', 'CusDirectKms', 'CusTravelledKms','MaxOcc'});
R_Pax = table([],[],[],[],[],[],[],[],[],[], 'VariableNames', {'id','nPax','Pax_maxRadius', 'BS_separation','maxWalkingDist',...
  'paxID','cusID','CusDirectKms','rideTime','walkTime'});

for i = 1:nFolders % this is the instance folder
  thisFolder = [data_to_load, items(i).name, filesep];
  folderItems = dir(thisFolder);
  allFiles = string({folderItems.name});

  cus_file = allFiles(contains(allFiles,'cus_')); route_files =  allFiles(startsWith(allFiles,'route_detail_'));
  if ~isempty(cus_file) && ~isempty(route_files) % we have some results in this instance

    % load all the usual scenario tables and params
    params_file = fullfile(thisFolder, allFiles(contains(allFiles,'parameters'))); p = readtable(params_file);
    filename = fullfile(thisFolder, allFiles(contains(allFiles,'passenger'))); T_Pax = readtable(filename);

    % allP = [allP;p]; % save these params
    T_C = readtable(fullfile(thisFolder, cus_file)); T_Pax = [T_Pax, T_C]; %#ok<AGROW>
    paxOK = T_Pax.ride_time>0; % these are the served customers

    this_Pax = table(i*ones(p.nPax,1), p.nPax*ones(p.nPax,1), p.Pax_maxRadius*ones(p.nPax,1), p.BS_separation*ones(p.nPax,1),p.maxWalkingDist*ones(p.nPax,1),...
      T_Pax.passenger_ID, T_Pax.cus, T_Pax.DistanceFromStation, T_Pax.ride_time, T_Pax.cus_walking_time,....
      'VariableNames', {'id','nPax','Pax_maxRadius', 'BS_separation','maxWalkingDist','paxID','cusID','CusDirectKms','rideTime','walkTime'});
    R_Pax = [R_Pax;this_Pax]; %#ok<AGROW>

    this_fleetSize = numel(route_files);
    % new instance so we zero all the totals
    total_veh_kms = 0; total_empty_kms = 0; total_cus_kms = 0; total_direct_kms = 0;
    for ff = 1:this_fleetSize % one for each bus. Could do multiple loops
      R = readtable(fullfile(thisFolder,route_files{ff}));  % You can use readmatrix for newer MATLAB versions
      this_leg_kms = sqrt(diff(R.x).^2 + diff(R.y).^2);
      this_leg_occ = R.num_passenger(1:end-1);
      this_cus_kms = this_leg_kms'*this_leg_occ; % empty bus automatically not counted
      this_empty_kms = this_leg_kms'*~this_leg_occ; % empty bus automatically not counted

      C = string(R.cus_set); C = C(~startsWith(C,"Int64")); % extract customer numbers on this bus
      this_bus_cus = cell2mat(cellfun(@str2double, regexp(C(:)', '\d+', 'match'), 'UniformOutput', false));

      this_direct_kms = sum(T_Pax.DistanceFromStation(ismember(T_Pax.cus, this_bus_cus)));
      total_veh_kms = total_veh_kms + sum(this_leg_kms);
      total_empty_kms = total_empty_kms + this_empty_kms;
      total_cus_kms = total_cus_kms + this_cus_kms;
      total_direct_kms = total_direct_kms + this_direct_kms;
      R_bus(end+1,:) = {i,p.nPax,p.Pax_maxRadius,p.BS_separation,p.maxWalkingDist,sum(this_leg_kms),this_empty_kms,this_direct_kms,this_cus_kms,max(this_leg_occ)}; %#ok<*SAGROW>

    end % loop bus routes
    R_instance{i,:} = [i,p.nPax,p.Pax_maxRadius,p.BS_separation,p.maxWalkingDist,this_fleetSize,total_veh_kms,total_empty_kms,total_direct_kms,total_cus_kms,p.stationSeparation];
  end % instance (non-empty)
end

save BICENTRIC_PROCESSED