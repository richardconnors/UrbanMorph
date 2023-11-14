function [T_busFleet, T_Passenger,T_busStop,T_Charger,T_Station,allStationDeps,T_depot] = generateScenario(params)

nStation = params.nStation;
Station_separation = params.stationSeparation;
nPax = params.nPax;
Pax_minRadius = params.Pax_minRadius;
Pax_maxRadius = params.Pax_maxRadius;
paxSeparation = params.paxSeparation;
maxWalkingDist = params.maxWalkingDist;
BS_separation = params.BS_separation;
nCharger = params.nCharger;
charger_radius = params.chargerRadius;
demandPeakness = params.demandPeakness;

% [T_busFleet, T_Passenger,T_busStop,T_Charger,T_Station,allStationDeps,T_depot] =  generateScenario(params)
% Generate scenarios for MEVRST algorithms to tackle
% area of interest will be a rectangle its size will be derived from the location of data points below
% set the centre at [0,0]

% ====== BUS FLEET
busType = 1; maxPax = 10; maxKWH = 35.8; % kWh - does this correspond to 100% SOC or maxSOC?
minSOC = 10; % percent
maxSOC = 80; % percent
SOC = 20; consumption = 0.24; %kW.km
nBus = ceil(0.5*nPax*0.7/maxPax);
busFleet1 = table(busType, maxPax, maxKWH, minSOC, maxSOC, SOC, consumption);
busFleet1 = repmat(busFleet1,nBus,1);
% distribute SOC amongst buses
busFleet1.SOC = linspace(20,80,nBus)'; % even spacing of bus SOC from 20 -> 80

busType = 2; maxPax = 20; maxKWH = 53.7; % kWh - does this correspond to 100% SOC or maxSOC?
minSoC = 10; % percent
maxSoC = 80; % percent
SoC = 20; consumption = 0.29; %kW.km
nBus = ceil(0.5*nPax*0.7/maxPax);
busFleet2 = table(busType,maxPax,maxKWH,minSOC,maxSOC,SOC,consumption);
busFleet2 = repmat(busFleet2,nBus,1);
busFleet2.SOC = linspace(20,80,nBus)'; % even spacing of bus SOC from 20 -> 80
T_busFleet = [busFleet1;busFleet2];

% ====== START network geometry
rng default % make randomness repeatable
x0 = 0; y0 = 0; % centre (could geolocate this later)

% ====== GENERATE LOCATIONS FOR TRANSIT, CHARGERS, PAX
% ====== how many stations = population centres? Must be >=1
station_XY = [x0,y0];
if nStation>1
  % random location around circle would use theta = rand(1, nTowns-1)*2*pi;
  % equidistant round circle
  theta = mod((2*pi/(nStation-1))*(1:nStation-1),2*pi);
  [this_x,this_y] = pol2cart(theta, Station_separation);
  station_XY = vertcat(station_XY,[this_x(:),this_y(:)]);
end
station_X = round(station_XY(:,1),5); % dont need more than 5 decimal places!
station_Y = round(station_XY(:,2),5);
station_ID = (1:nStation)';
T_Station = table(station_ID, station_X,station_Y);

% ========================================
% ======== DEPARTURES FROM EACH STATION ??
bufferT = 15;
headway1 = 20;
depOffset = 10;
headwayX = 15;
% for each transit station we have a list of departure times
% generate here for now - could load from elsewhere
Tstart = datetime(2022,10,3,6,0,0); Tstart.Format = 'yyyy-MM-dd HH:mm';
Tend = datetime(2022,10,3,10,00,0); Tend.Format = 'yyyy-MM-dd HH:mm';

% First station timetable
headway = headway1;
thisClockTimes = (Tstart:minutes(headway):Tend)';  % clock time departures
nDeps = size(thisClockTimes,1);
TstartRel = headway;
thisRelTimes = (TstartRel:headway:nDeps*headway)'; % relative time departures in minutes
allStationDeps = table(ones(nDeps,1),thisClockTimes,thisRelTimes-bufferT,thisRelTimes,'VariableNames',{'StationID','DepClockTime','ei','li'});

for t = 2:nStation
  headway = headwayX;
  Tstart = Tstart+minutes(depOffset); % increments by depOffset each time round the loop
  TstartRel = TstartRel + depOffset; % increments by depOffset each time round the loop
  thisClockTimes = (Tstart+minutes(depOffset):minutes(headway):Tend)';
  nDeps = size(thisClockTimes,1);
  thisRelTimes = TstartRel+(0:nDeps-1)'*headway;
  thisStationDeps = table(t*ones(nDeps,1),thisClockTimes,thisRelTimes-bufferT,thisRelTimes,'VariableNames',{'StationID','DepClockTime','ei','li'});
  allStationDeps = [allStationDeps;thisStationDeps];
end

allStationDeps = sortrows(allStationDeps,'li');
allStationDeps.layer = (1:height(allStationDeps))';

% ====== chargers
charger_XY = station_XY; % if only 1 charger per station
if nCharger>1
  theta = mod((2*pi/(nCharger-1))*(1:nCharger-1),2*pi);
  [this_x,this_y] = pol2cart(theta, charger_radius);
  for t = 1:nStation
    % equidistant round circle
    charger_XY = vertcat(charger_XY,[station_XY(t,1)+this_x(:),station_XY(t,2)+this_y(:)]);
  end
end
nCharger = size(charger_XY,1);
charger_ID = (1:nCharger)';
charger_X = round(charger_XY(:,1),5); % dont need more than 5 decimal places!
charger_Y = round(charger_XY(:,2),5);
% chargerX = [0;2];chargerY = [0;0];
charger_rate = 0.83*ones(nCharger,1); % kW/minute
T_Charger = table(charger_ID,charger_X,charger_Y,charger_rate);

% ========================================
% ====== generate passenger locations
% nPax is split amongst the stations
% For each station we generate nPax/nStations around that station (min/max
% radii) AND these pax go to THAT station.
pax_XY = []; pax_stationID = []; DistanceFromStation = [];
maxIter = 1000; % for algorithm ensuring pax min separation
this_Pax_minRadius = Pax_minRadius(1);
this_Pax_maxRadius = Pax_maxRadius(1);
for t = 1:nStation
  % generatePassengers(x0, y0, minRadius, maxRadius, nPax, minPaxDist, maxIter)
  if length(Pax_minRadius)==nStation
    this_Pax_minRadius = Pax_minRadius(t);
  end
  if length(Pax_maxRadius)==nStation
    this_Pax_maxRadius = Pax_maxRadius(t);
  end
  this_XY = generatePassengers(station_XY(t,1), station_XY(t,2), this_Pax_minRadius, this_Pax_maxRadius, floor(nPax./nStation), paxSeparation, maxIter);
  this_N = size(this_XY,1);
  this_ID = t*ones(this_N,1);
  pax_XY = [pax_XY;this_XY];
  pax_stationID = [pax_stationID;this_ID(:)]; %#ok<*AGROW>

  this_DfromO = sqrt((this_XY(:,1) - repmat(station_XY(t,1),this_N,1)).^2 +...
    (this_XY(:,2) - repmat(station_XY(t,2),this_N,1)).^2);
  DistanceFromStation = [DistanceFromStation; this_DfromO];
end
nPax = size(pax_XY,1);
% if we want to send passengers to random station locations
% pax_stationID = randi(nStation,[nPax,1]);

% ====== end GENERATE LOCATIONS FOR PAX, TRANSIT, CHARGERS




% ========================================
% ========================================
% ====== GENERATE DEMAND ASSIGNMENTS
% generate passenger destinations = transit station and train dep time
% assume pax go to local transit station
% assume same demand profile for each station

% demandPeakness: 0 = uniform, 1 = peaked demand  2 = all on single dep
pax_depT = nan(nPax,1);
for t = 1:nStation
  % distribute demand (pax with this station ID)
  % among departures from THIS station only
  % according to demand peakedness
  thisDeps = allStationDeps{allStationDeps.StationID==t,'li'};
  nDeps = length(thisDeps);
  pax_i = pax_stationID==t; % pax using this station
  if demandPeakness >=0 && demandPeakness<= 1
    % to get desired demand profile
    % what prob for each departure time?
    pUnif = ones(nDeps,1)*1/nDeps;
    pNorm = normpdf(1:nDeps, mean(1:nDeps),1); pNorm = pNorm./sum(pNorm);
    pDeps = (1-demandPeakness)*pUnif + demandPeakness*pNorm(:); % prob of a pax taking each departure
  else
    pDeps = zeros(nDeps,1); pDeps(end) = 1;
  end
  thisPaxDeps = gendist(pDeps,1,sum(pax_i));
  pax_depT(pax_i) = thisPaxDeps(:);
end
% give the fields nice names for the table
passenger_ID = (1:nPax)';
passenger_X = round(pax_XY(:,1),5);
passenger_Y = round(pax_XY(:,2),5);
passenger_StationID = pax_stationID;
passenger_DepartureTime = pax_depT;

T_Passenger = table(passenger_ID, passenger_X, passenger_Y, passenger_StationID, passenger_DepartureTime, DistanceFromStation);

% ====== end GENERATE DEMAND ASSIGNMENTS
% ========================================

% ========================================
% ====== GENERATE BUS STOPS
% Grid of all legally viable bus stop locations
% BS grid separation needs to match max walking distance
if 0.5*BS_separation*sqrt(2)>maxWalkingDist+eps % to deal with rounding error
  disp('Warning: max BS separation > walking distance');
end

% get bounding box encompassing pax locations
minX = min(pax_XY(:,1)) - BS_separation;
maxX = max(pax_XY(:,1)) + BS_separation;
minY = min(pax_XY(:,2)) - BS_separation;
maxY = max(pax_XY(:,2)) + BS_separation;
% Generate a regular square grid of points within the bounding box
[BS_X,BS_Y] = meshgrid(minX:BS_separation:maxX, minY:BS_separation:maxY);
% Convert the grid points to a list of [x, y] coordinates
BS_XY = [BS_X(:),BS_Y(:)];

% for each customer find nearest meeting point
% check within walking distance
[nearestBS,dist] = dsearchn(BS_XY,pax_XY);
paxMissed = dist>maxWalkingDist;

% remove potential bus stops more than walking distance away from any pax
% for each bus stop search amongst ALL pax locations to find nearest pax
% index those BS with nearest pax within maxWalkingDist
[~,dist] = dsearchn(pax_XY,BS_XY);
BSinReach = dist<=maxWalkingDist;
% check that every BS that was the closest to a person is still included
if ~all(BSinReach(nearestBS))
  disp('Closest BS has been removed'); keyboard
end

nBS = sum(BSinReach);
busStop_ID = (1:nBS)';
busStop_X = round(BS_XY(BSinReach,1),5);
busStop_Y = round(BS_XY(BSinReach,2),5);
T_busStop = table(busStop_ID,busStop_X,busStop_Y);
% ====== end GENERATE BUS STOPS

depot_ID = 1;
depot_X = round(mean(busStop_X),5);
depot_Y = round(mean(busStop_Y),5);
% for Arlon
% depot_X = 0; depot_Y = -1;
T_depot = table(depot_ID,depot_X,depot_Y);

end
