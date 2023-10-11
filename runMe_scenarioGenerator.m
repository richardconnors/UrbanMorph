
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
maxWalkingDist = 1.2;%BS_separation*sqrt(2)/2; %
nCharger  = 4; % how many chargers PER TOWN/TRANSIT STATION.
charger_radius = 2; % chargers located on circle around each town centre
demandPeakness = 0; % 0 = uniform. 1 = peaked in middle
PLOTFLAG = 1; % plots the network data




% the plot shows 
% station & chargers
% passengers (with text showing their desired departure number)
% the grid of (potential) bus stop locations

generateScenario(nTransitStations, Station_separation, nPassengers, Pax_minRadius, Pax_maxRadius,...
  paxSeparation, maxWalkingDist, BS_separation, nCharger, charger_radius, demandPeakness, PLOTFLAG)

