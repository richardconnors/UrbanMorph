% data should be in ONE master folder that contains N subfolders
% each with a different 'run' of some sort
% This script will iterate through EACH results subfolder, extract data for
% that 'run' and generate summary statistics.

% sort out folder and file names
% repo_flexbus = [get_repo_folder, 'Flexbus3_v0.8.4\'];
% data_to_load = [repo_flexbus, 'results\grid_results_v2\'];
data_to_load = 'Q:\REPOS\UrbanMorph\data\';


% get list of all folders we need to iterate through & sift them into  viable and non-viable folders
badDataDir = 'Q:\REPOS\UrbanMorph\failed_data\';
items = dir(data_to_load);
% Check if the item is a directory and not '.' or '..'
items = items([items.isdir]); items = items(~ismember({items.name}, {'.', '..'}));
nInstance = numel(items); % this is the instance folder
for i = 1:nInstance % this is the instance folder
  thisFolder = [data_to_load, items(i).name, filesep];
  folderItems = dir(thisFolder); allFiles = string({folderItems.name});
  cus_file = allFiles(contains(allFiles,'cus_'));
  route_files =  allFiles(startsWith(allFiles,'route_detail_'));
  if isempty(cus_file) || isempty(route_files) % we have some results in this instance
    % move this folder to failed_instances
    movefile(thisFolder,badDataDir);
  end
end

items = dir(data_to_load);
% Check if the item is a directory and not '.' or '..'
items = items([items.isdir]); items = items(~ismember({items.name}, {'.', '..'}));
nInstance = numel(items); % this is the instance folder
allP = table;
for i = 1:nInstance % this is the instance folder
  thisFolder = [data_to_load, items(i).name, filesep];
  folderItems = dir(thisFolder); allFiles = string({folderItems.name});
  cus_file = allFiles(contains(allFiles,'cus_'));
  route_files =  allFiles(startsWith(allFiles,'route_detail_'));
  if ~isempty(cus_file) && ~isempty(route_files) % we have some results in this instance

    % load all the usual scenario tables and params
    params_file = fullfile(thisFolder, allFiles(contains(allFiles,'parameters'))); p = readtable(params_file);
    filename = fullfile(thisFolder, allFiles(contains(allFiles,'passenger'))); T_Passenger = readtable(filename);
    % filename = fullfile(thisFolder, allFiles(contains(allFiles,'busFleet'))); T_busFleet = readtable(filename);
    % filename = fullfile(thisFolder, allFiles(contains(allFiles,'busStopXY'))); T_busStop = readtable(filename);
    % filename = fullfile(thisFolder, allFiles(contains(allFiles,'chargerXY')));  T_Charger = readtable(filename);
    % filename = fullfile(thisFolder, allFiles(contains(allFiles,'stationXY'))); T_Station = readtable(filename);
    % filename = fullfile(thisFolder, allFiles(contains(allFiles,'transitTimetable'))); allStationDeps = readtable(filename);
    % filename = fullfile(thisFolder, allFiles(contains(allFiles,'depotXY'))); T_depot = readtable(filename);

    allP = [allP;p];
    cus_data = readtable(fullfile(thisFolder, cus_file));  % You can use readmatrix for newer MATLAB versions
    allCus = cus_data.cus;
    % these customers should be in the same order as in T_Passenger so we
    % can match the ride time etc
    T_Passenger.RideTime = cus_data.ride_time; T_Passenger.WalkTime = cus_data.cus_walking_time;
    this_rideTime = cus_data.ride_time; this_rideTime(~this_rideTime) = NaN;
    paxOK = cus_data.ride_time>0; % these are the served customers
    
    busOccupancyData = []; % for histogram of bus occupancy per trip
    servedCus = []; % which customers are served/missed
    allTours = {}; % want to plot a line for each tour
    tourDist = []; tourLegs_kms = []; tourLegs_Occ = [];
    fleetSize(i) = numel(route_files);
    for ff = 1:fleetSize(i) % one for each bus?
      R = readtable(fullfile(thisFolder,route_files{ff}));  % You can use readmatrix for newer MATLAB versions
      busOccupancyData = [busOccupancyData; R.num_passenger];
      allTours{end+1} = [R.x, R.y]; %#ok<*SAGROW>

      tourLegs_kms = [tourLegs_kms; sqrt(diff(R.x).^2 + diff(R.y).^2)];
      tourLegs_Occ = [tourLegs_Occ; R.num_passenger(1:end-1)];
      tourDist = [tourDist; sum(tourLegs_kms)];
      

      C = string(R.cus_set); C = C(~startsWith(C,"Int64"));
      pattern = '\d+';
      % Loop through each element of the string array
      for ii = 1:numel(C)
        % Find all numerical values in the current string
        matches = regexp(C(ii), pattern, 'match');
        % Convert the matched values to a numeric array and append them to the vector
        servedCus = [servedCus, matches]; %#ok<*AGROW>
      end
      servedCus = unique(str2double(servedCus));
    end % loop bus routes
    totalTourDist(i) = sum(tourDist);
    nServed(i) = sum(paxOK);
    totalDirectDistance(i) = sum(T_Passenger.DistanceFromStation(paxOK));
  end % instance (non-empty)
end

% Total veh kms
figure;
h1 = scatter3(allP.nPax, allP.Pax_maxRadius, 2*totalDirectDistance);
hold on
xlabel('nCustomers'); ylabel('City Radius'); zlabel('KPI')
h2 = scatter3(allP.nPax, allP.Pax_maxRadius, totalTourDist);
legend({'Total Direct Distance','Total Tour Distance'})
cmap=colormap_generator(2);
h1.SizeData = 40;
h1.MarkerFaceColor = cmap(1,:);
h1.MarkerFaceAlpha = 0.7;
h2.SizeData = 40;
h2.MarkerFaceColor = cmap(2,:);
h2.MarkerFaceAlpha = 0.7;

% relative vehkms compared with taxi
figure;
h1 = scatter3(allP.nPax, allP.Pax_maxRadius, totalDirectDistance./totalTourDist);
xlabel('nCustomers'); ylabel('City Radius'); zlabel('KPI')
title('Taxi Distance / DRT distance')
cmap=colormap_generator(2);
h1.SizeData = 40;
h1.MarkerFaceColor = cmap(1,:);
h1.MarkerFaceAlpha = 0.7;

% fleet size
figure;
h1 = scatter3(allP.nPax, allP.Pax_maxRadius, fleetSize);
xlabel('nCustomers'); ylabel('City Radius'); zlabel('KPI')
title('Fleet Size')
cmap=colormap_generator(2);
h1.SizeData = 40;
h1.MarkerFaceColor = cmap(1,:);
h1.MarkerFaceAlpha = 0.7;

fh = boxPlot3D(fleetSize, allP.nPax, allP.Pax_maxRadius);
xlabel('nCustomers'); ylabel('City Radius'); zlabel('Fleet Size')

% how busy or well used are buses?
% maybe plot nuber of passenger kms compared with num bus kms?


% chargeable kms vs bus kms


