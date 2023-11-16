% data should be in ONE master folder that contains N subfolders
% each with a different 'run' of some sort
% This script will iterate through EACH results subfolder, extract data for
% that 'run' and generate summary statistics.

% sort out folder and file names
repo_flexbus = [get_repo_folder, 'Flexbus3_v0.8.4\'];
data_to_load = [repo_flexbus, 'results\urbanMorph_results\everything_singleStation\'];

% get list of all folders we need to iterate through
items = dir(data_to_load);
% Check if the item is a directory and not '.' or '..'
items = items([items.isdir]); items = items(~ismember({items.name}, {'.', '..'}));
nInstance = numel(items); % this is the instance folder
% save "coordinates" 
% nPax  BS_separation maxWalkingDist Pax_minRadius Pax_maxRadius nStation
allP = table;
for i = 1:nInstance % this is the instance folder
  thisFolder = [data_to_load, items(i).name, filesep];
  folderItems = dir(thisFolder); allFiles = string({folderItems.name});
  cus_file = allFiles(contains(allFiles,'cus_'));
  route_files =  allFiles(startsWith(allFiles,'route_detail_'));
  if ~isempty(cus_file) && ~isempty(route_files) % we have some results in this instance

    % load all the usual scenario tables and params
    params_file = fullfile(thisFolder, allFiles(contains(allFiles,'parameters')));
    p = readtable(params_file);
    filename = fullfile(thisFolder, allFiles(contains(allFiles,'passenger')));
    T_Passenger = readtable(filename);
    filename = fullfile(thisFolder, allFiles(contains(allFiles,'busFleet')));
    T_busFleet = readtable(filename);
    filename = fullfile(thisFolder, allFiles(contains(allFiles,'busStopXY')));
    T_busStop = readtable(filename);
    filename = fullfile(thisFolder, allFiles(contains(allFiles,'chargerXY')));
    T_Charger = readtable(filename);
    filename = fullfile(thisFolder, allFiles(contains(allFiles,'stationXY')));
    T_Station = readtable(filename);
    filename = fullfile(thisFolder, allFiles(contains(allFiles,'transitTimetable')));
    allStationDeps = readtable(filename);
    filename = fullfile(thisFolder, allFiles(contains(allFiles,'depotXY')));
    T_depot = readtable(filename);

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
    for ff = 1:length(route_files) % one for each bus?
      R = readtable(fullfile(thisFolder,route_files{ff}));  % You can use readmatrix for newer MATLAB versions
      busOccupancyData = [busOccupancyData; R.num_passenger];
      allTours{end+1} = [R.x, R.y]; %#ok<*SAGROW>

      tourLegs_kms = [tourLegs_kms; sqrt(diff(R.x).^2 + diff(R.y).^2)];
      tourLegs_Occ = [tourLegs_Occ; R.num_passenger(1:end-1)];
      tourDist = [tourDist; sum(leg_kms)];
      

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

figure;
h1 = scatter3(allP.nPax, allP.Pax_maxRadius, totalDirectDistance./totalTourDist);
xlabel('nCustomers'); ylabel('City Radius'); zlabel('KPI')
title('Taxi Distance / DRT distance')
cmap=colormap_generator(2);
h1.SizeData = 40;
h1.MarkerFaceColor = cmap(1,:);
h1.MarkerFaceAlpha = 0.7;



return

figure('pos',[ 162,524,1511,420])
boxplot(T_Passenger{ind,7:end}'); ylabel('Ride Distance')
yyaxis right
h = scatter(1:100,T_Passenger{ind,6}); ylabel('Distance to station')
title('Distribution of ride time with distance from station')

figure;
histogram(busOccupancyData(busOccupancyData>0))
title('Bus Occupancy'); ylabel('Number of Trips')

figure; % Plot all tours as patches
h_bs = scatter(T_busStop.busStop_X,T_busStop.busStop_Y,'+');
h_bs.MarkerEdgeColor = 0.8*ones(3,1);
hold on
h_st = scatter(0,0,'gs'); h_st.MarkerFaceColor = 'g';
h_pax = scatter(T_Passenger.passenger_X,T_Passenger.passenger_Y,'r.');
fh = cell(numel(allTours),1);
for tt = 1:numel(allTours)
  thisTour = allTours{tt};
  fh{tt} = fill(thisTour(:,1),thisTour(:,2),0.9*[1,1,1]);
  fh{tt}.FaceAlpha = 0.2;
  fh{tt}.FaceColor = 0.6*[1,1,1];
  fh{tt}.EdgeColor = "none";
end

return

figure;
h_bs = scatter(T_busStop.busStop_X,T_busStop.busStop_Y,'+');
h_bs.MarkerEdgeColor = 0.8*ones(3,1);
hold on
h_st = scatter(0,0,'gs'); h_st.MarkerFaceColor = 'g';
h_pax = scatter(T_Passenger.passenger_X,T_Passenger.passenger_Y,'r.');
for i=1:10
  h_pax = scatter(T_Passenger.passenger_X(pSchedule(:,i)),T_Passenger.passenger_Y(pSchedule(:,i)),'bo');
  title(sprintf('Day %02d',i));
  pause(0.5)
  delete(h_pax)
end
