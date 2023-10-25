% data should be in ONE master folder that contains N subfolders
% each with a different 'run' of some sort
% This script will iterate through EACH results subfolder, extract data for
% that 'run' and generate summary statistics.

% sort out folder and file names
repo_flexbus = [get_repo_folder, 'Flexbus3_v0.7\'];
data_to_load = [repo_flexbus, 'data\results\BS_separation\'];


items = dir(data_to_load);
items = items([items.isdir]);
indexToKeep = ~ismember({items.name}, {'.', '..'});
items = items(indexToKeep);
for i = 1:numel(items)
  thisFolder = [data_to_load, items(i).name, filesep];
  matFiles = dir([thisFolder,'*.mat']);
  load([thisFolder ,matFiles(1).name]) % should be the workspace
  disp('=================================')
  disp(['Loading ' matFiles(1).name])

  rFiles = dir([thisFolder,'route_detail*.csv']);
  cFile = dir([thisFolder,'cus_experience*.csv']);
  vdfFile = dir([thisFolder,'V_DF.csv']);

  T = T_Passenger;
  busOccupancyData = []; % for histogram of bus occupancy per trip
  allTours = {}; % want to plot a line for each tour

  if ~isempty(vdfFile) && ~isempty(cFile) && numel(rFiles)>0
    % display brief param descriptions
    fprintf('nPassengers = %d \n', nPassengers);
    fprintf('Pax_maxRadius = %.1f \n', Pax_maxRadius)
    fprintf('Pax_minRadius = %.1f \n', Pax_minRadius)
    fprintf('demandPeakness = %.1f \n', demandPeakness)

    vdf = readtable([thisFolder, vdfFile.name]);
    this_pax = vdf(strcmp(vdf.v_type,'customer'),:); % has columns of x and y data
    % we assume that this_pax are T_passengers in order.
    % lets double check
    tolerance = 0.001;
    XareEqual = all(abs(this_pax.x - T_Passenger.passenger_X) < tolerance);
    YareEqual = all(abs(this_pax.y - T_Passenger.passenger_Y) < tolerance);
    if ~XareEqual || ~YareEqual
      disp('Customers in V_DF do not match T_Passenger')
    end
    pax_vdfID = this_pax.id; % correspond to original customer numbers 1:N

    cus_data = readtable([thisFolder, cFile.name]);  % You can use readmatrix for newer MATLAB versions
    rideTime = cus_data.ride_time;
    walkTime = cus_data.cus_walking_time;
    rejectedCusInd = ~rideTime;
    rideTime(rejectedCusInd) = NaN; % customer was rejected
    walkTime(rejectedCusInd) = NaN;
    totalDistance = 0; totalPaxTime = 0;
    totalDriveTime = 0; % will include charging, waiting etc
    for ff = 1:numel(rFiles) % one for each bus?
      this_route = readtable(fullfile(thisFolder,rFiles(ff).name));  % You can use readmatrix for newer MATLAB versions

      % total distance travelled
      route_xy = [this_route.x,this_route.y];
      diffs = diff(route_xy, 1, 1);
      distances = sqrt(sum(diffs.^2, 2)); % between consecutive points
      % Sum all the distances to get the total distance traveled
      totalDistance = totalDistance + sum(distances);
      totalDriveTime = totalDriveTime + this_route.t_arr(end);

      % ride time for each customer
      % so inherently includes multiplying travel time by occupancy
      totalPaxTime = totalPaxTime + sum(this_route.ride_time);

      allTours{end+1} = [this_route.x, this_route.y]; %#ok<*SAGROW>

    end
  else
    disp(['No solutions in ', items(i).name])
  end
end %loop through all instance folders


return

% what do we want to analyse?
[~,ind] = sort(T.DistanceFromStation);
figure('pos',[ 162,524,1511,420])
boxplot(T{ind,7:end}'); ylabel('Ride Distance')
yyaxis right
h = scatter(1:100,T{ind,6}); ylabel('Distance to station')
title('Distribution of ride time with distance from station')

% occupancy of bus on each legs
figure;
histogram(busOccupancyData(busOccupancyData>0))
title('Bus Occupancy'); ylabel('Number of Trips')

% Plot all tours as patches on same figure
figure; 
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

% scatter plot of bus stops, passenger locations, etc
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
