
% generateScenario allows the inputs listed below to be specified.
% You can set values here, and run this script to call the function
% Outputs are written to comma delimited text files (easy to change this)
% Output text file names currently hardcoded as below, saved to PWD
% % transitTimetable.txt
% % stationXY.txt
% % chargerXY.txt
% % paxData.txt
% % busStopXY.txt

% % All distances in KM
% % Timetable of departures for each station is hardcoded
% % For each passenger departure time is randomly assigned
% % Demand distribution is controlled and can go from uniform -> peaked
% nTransitStations = 2;
% Station_separation = 5; % distance between multiple stations
% nPassengers = 100; % number passengers to generate
% % Passengers will be generated uniformly distributed in an annulus around each station
% % minRadius and maxRadius can have single value = used for all stations
% % OR can be a row vector size [1 x nTransitStations]
% Pax_minRadius = 0.5; % min radius away from station for passenger locations
% Pax_maxRadius = 2.0; % max radius around station for passenger locations
% paxSeparation = 0.05; % min distance between passengers
% BS_separation = 1; % spacing for grid of potential bus stop locations
% maxWalkingDist = BS_separation*1.05/sqrt(2); %
% nCharger  = 4; % how many chargers PER TOWN/TRANSIT STATION.
% charger_radius = 1.5; % chargers located on circle around each town centre
% demandPeakness = 0; % 0 = uniform. 1 = peaked in middle
% PLOTFLAG = 1; % plots the network data


nTransitStations = 1;
Station_separation = 1; % distance between multiple stations
nPassengers = 100; % number passengers to generate
% Passengers will be generated uniformly distributed in an annulus around each station
% minRadius and maxRadius can have single value = used for all stations
% OR can be a row vector size [1 x nTransitStations]
Pax_minRadius = 1.5; % min radius away from station for passenger locations
Pax_maxRadius = 5.0; % max radius around station for passenger locations
paxSeparation = 0.05; % min distance between passengers
BS_separation = 1; % spacing for grid of potential bus stop locations
maxWalkingDist = 1.2; % BS_separation*sqrt(2)/2; %
nCharger  = 4; % how many chargers PER TOWN/TRANSIT STATION.
charger_radius = 2; % chargers located on circle around each town centre
demandPeakness = 0; % 0 = uniform. 1 = peaked in middle
PLOTFLAG = 0; % plots the network data


% here generate the total scenario with all passengers
[T_busFleet, T_Passenger,T_busStop,T_Charger,T_Station,allStationDeps,T_depot] = generateScenario(nTransitStations, Station_separation, nPassengers, Pax_minRadius, Pax_maxRadius,...
  paxSeparation, maxWalkingDist, BS_separation, nCharger, charger_radius, demandPeakness, PLOTFLAG);

% now subsample the T_passenger and write to yumeng style files.
% will need to generate different filenames?
repo_folder = 'Q:\REPOS\Flexbus3_v0.7\';
[~, hostname] = system('hostname');
if contains(hostname, 'PP0695')
  repo_folder = 'C:\Users\richard.connors\Documents\REPOS\Flexbus3_v0.7\';
end

nDays = 50;
nSample = 10;
pSchedule = zeros(nPassengers,nDays);
for i = 1:nDays
  this_Pax = randperm(nPassengers,nSample);
  pSchedule(this_Pax,i) = 1;
  T_thisPax = T_Passenger(this_Pax,:);
  % do I need to re-number passengers to be consecutive?
  
  saveFolder = [repo_folder, 'data\test\', sprintf('P%d_%dD%d_%d',nPassengers,nSample,nDays,i)];
  if ~isfolder(saveFolder), mkdir(saveFolder); end
  saveToYumengFormat(saveFolder,T_busFleet, T_thisPax ,T_busStop,T_Charger,T_Station,allStationDeps,T_depot,maxWalkingDist)
end




