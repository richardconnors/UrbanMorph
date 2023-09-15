function generateScenario(nStation, Station_separation, nPax, Pax_minRadius, Pax_maxRadius, paxSeparation, maxWalkingDist, BS_separation, nCharger, charger_radius, demandPeakness, PLOTFLAG)
% Generate scenarios for MEVRST algorithms to tackle
% area of interest will be a rectangle its size will be derived from the location of data points below
% set the centre at [0,0]

arguments % declare arguments and supply default values
  nStation (1,1) double = 3;
  Station_separation (1,1) double = 3; % distance between multiple centres
  nPax (1,1) double = 50;
  Pax_minRadius (1,:) double = 0.6; % min radius around town centre for passenger locations
  Pax_maxRadius (1,:) double = 2; % max radius around town centre for passenger locations
  paxSeparation (1,1) double = 50; % min distance between passengers
  maxWalkingDist (1,1) double = 1.5; %
  BS_separation (1,1) double = .25; % spacing for grid of potential bus stop locations
  nCharger (1,1) double = 4; % how many chargers PER TOWN.
  charger_radius (1,1) double = 1; % radius of circle around each town centre
  demandPeakness  (1,1) double = 0; % 0 = uniform. 1 = peaked in middle
  PLOTFLAG (1,1) logical = 1; % plot data (1) or not (0)?
end

% ===== SAVE LOCATION
saveFolder = sprintf('P%dS%dC%dDP%.1f',nPax,nStation,nCharger,demandPeakness);
if ~isfolder(saveFolder), mkdir(saveFolder); end


% ====== PARAMETERS
params.maxWalkingDist = round(maxWalkingDist,5);
params.walkingSpeed = 0.085; % km/minute
params.busSpeed = 0.83; % km/minute
params.busStopServiceTime = 0.5; % minutes
params.maxTransitWaitingTime = 15; % minutes
params.nPax = nPax;
params.paxSeparation = paxSeparation;
params.nStation = nStation;
params.nCharger = nCharger;
params.chargerRadius = charger_radius;
params.demandPeakness = demandPeakness;
params = orderfields(params);
writetable(struct2table(params),[saveFolder '\parameters.txt'],'Delimiter',',');
writestruct(params,[saveFolder '\parameters.xml']);


% ====== BUS FLEET
busType = 1;
maxPax = 10;
maxKWH = 35.8; % kWh - does this correspond to 100% SOC or maxSOC?
minSOC = 10; % percent
maxSOC = 80; % percent
SOC = 20;
consumption = 0.24; %kW.km
nBus = ceil(0.5*nPax*0.7/maxPax);
busFleet1 = table(busType, maxPax, maxKWH, minSOC, maxSOC, SOC, consumption);
busFleet1 = repmat(busFleet1,nBus,1);
% distribute SOC amongst buses
busFleet1.SOC = linspace(20,80,nBus)'; % even spacing of bus SOC from 20 -> 80

busType = 2;
maxPax = 20;
maxKWH = 53.7; % kWh - does this correspond to 100% SOC or maxSOC?
minSoC = 10; % percent
maxSoC = 80; % percent
SoC = 20;
consumption = 0.29; %kW.km
nBus = ceil(0.5*nPax*0.7/maxPax);
busFleet2 = table(busType,maxPax,maxKWH,minSOC,maxSOC,SOC,consumption);
busFleet2 = repmat(busFleet2,nBus,1);
busFleet2.SOC = linspace(20,80,nBus)'; % even spacing of bus SOC from 20 -> 80
busFleet = [busFleet1;busFleet2];
writetable(busFleet,[saveFolder '\busFleet.txt'],'Delimiter',',');

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

% ======== DEPARTURES FROM EACH STATION ??
% for each transit station we have a list of departure times
% generate here for now - could load from elsewhere
Tstart = datetime(2022,10,3,6,0,0); Tstart.Format = 'yyyy-MM-dd HH:mm';
Tend = datetime(2022,10,3,10,00,0); Tend.Format = 'yyyy-MM-dd HH:mm';
headway = minutes(20);
T_depTimes = table((Tstart:headway:Tend)','VariableNames',{'Station01'}); % dep times at Station 1
nDeps = size(T_depTimes,1);
for t = 2:nStation
  thisDeps = (Tstart+(t-1)*minutes(10):headway:Tend)';
  thisDeps = [thisDeps;NaT(nDeps-length(thisDeps),1)]; % ensure all cols same length
  thisCol = table(thisDeps,'VariableNames',{strcat('Station',num2str(t,'%02d'))});
  T_depTimes = [T_depTimes,thisCol];
end

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
charger_rate = 0.83*ones(nCharger,1); % kW/minute
T_Charger = table(charger_ID,charger_X,charger_Y,charger_rate);

% ====== generate passenger locations
% Create a set of points.
pax_XY = []; pax_stationID = [];
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
  this_ID = t*ones(size(this_XY,1),1);
  pax_XY = [pax_XY;this_XY];
  pax_stationID = [pax_stationID;this_ID(:)]; %#ok<*AGROW>
end
nPax = size(pax_XY,1);
% ====== end GENERATE LOCATIONS FOR PAX, TRANSIT, CHARGERS

% ====== GENERATE DEMAND ASSIGNMENTS
% generate passenger destinations = transit station and train dep time
% assume pax go to local transit station
% assume same demand profile for each station

% demandPeakness: 0 = uniform, 1 = very peaked demand
pax_depT = [];
for t = 1:nStation
  thisDeps = T_depTimes{:,t}; thisDeps=thisDeps(~isnat(thisDeps));
  nDeps = length(thisDeps);

  % what prob for each departure to get desired demand profile?
  pUnif = ones(nDeps,1)*1/nDeps;
  pNorm = normpdf(1:nDeps, mean(1:nDeps),1); pNorm = pNorm./sum(pNorm);
  pDeps = (1-demandPeakness)*pUnif + demandPeakness*pNorm(:);
  thisPaxDeps = gendist(pDeps,1,sum(pax_stationID==t));
  pax_depT = [pax_depT; thisPaxDeps(:)];

end
passenger_ID = (1:nPax)'; 
passenger_X = round(pax_XY(:,1),5);
passenger_Y = round(pax_XY(:,2),5);
passenger_StationID = pax_stationID;
passenger_DepartureTime = pax_depT;

T_Passenger = table(passenger_ID, passenger_X, passenger_Y, passenger_StationID, passenger_DepartureTime);


% ====== end GENERATE DEMAND ASSIGNMENTS


% ====== GENERATE BUS STOPS
% Grid of all legally viable bus stop locations
% BS grid separation needs to match max walking distance
if 0.5*BS_separation*sqrt(2)>maxWalkingDist,   disp('Warning: BS separation >> walking distance'); end

% generate BS over grid extending to all pax locations
% extent = max(abs(pax_XY(:))); % distance in m across all pax locations
% BS_x = x0-extent-BS_separation:BS_separation:x0+extent+BS_separation;
% BS_y = y0-extent-BS_separation:BS_separation:y0+extent+BS_separation;
biggestX = max(abs(pax_XY(:,1)));
BS_x = unique([0:-BS_separation:-biggestX-BS_separation,0:BS_separation:biggestX+BS_separation]);
BS_x = BS_x(BS_x>=min(pax_XY(:,1))-BS_separation);
BS_x = BS_x(BS_x<=max(pax_XY(:,1))+BS_separation);

biggestY = max(abs(pax_XY(:,2)));
BS_y = unique([0:-BS_separation:-biggestY-BS_separation,0:BS_separation:biggestY+BS_separation]);
BS_y = BS_y(BS_y>=min(pax_XY(:,2))-BS_separation);
BS_y = BS_y(BS_y<=max(pax_XY(:,2))+BS_separation);


[BS_X,BS_Y] = meshgrid(BS_x,BS_y);
BS_XY = [BS_X(:),BS_Y(:)];

% remove potential bus stops more than walking distance away from pax
% for each bus stop search amongst ALL pax locations to find nearest pax
% index those BS with nearest pax within maxWalkingDist
[~,dist] = dsearchn(pax_XY,BS_XY);
BSinCH = dist<=maxWalkingDist;
% find bus stops inside convex hull of pax locations (with 100m tolerance)
% BSinCH = inhull(BS_XY,pax_XY,[],100);
nBS = sum(BSinCH);

busStop_ID = (1:nBS)';
busStop_X = round(BS_XY(BSinCH,1),5);
busStop_Y = round(BS_XY(BSinCH,2),5);
T_busStop = table(busStop_ID,busStop_X,busStop_Y);
% ====== end GENERATE BUS STOPS

depot_ID = 1;
depot_X = round(mean(busStop_X),5);
depot_Y = round(mean(busStop_Y),5);
T_depot = table(depot_ID,depot_X,depot_Y);



% % ====== SAVE TO TEXT FILES
% % note by default saves 15 decimal places which seems excessive hence use round
% writetable(T_depTimes,[saveFolder '\transitTimetable.txt'],'Delimiter',',');
% writetable(T_Station,[saveFolder '\stationXY.txt'],'Delimiter',',')
% writetable(T_Charger,[saveFolder '\chargerXY.txt'],'Delimiter',',')
% % pax data columns: x-cood, y-cood, stationID, depTimeIndex
% writetable(T_Passenger,[saveFolder '\passengerData.txt'],'Delimiter',',')
% writetable(T_busStop,[saveFolder '\busStopXY.txt'],'Delimiter',',')
% writetable(T_depot,[saveFolder '\depotXY.txt'],'Delimiter',',')

% ====== end SAVE TO TEXT FILES



% ======= PLOTTING DATA IF DESIRED =====================
if PLOTFLAG
%   figure; hh = histogram(pax_depT);
%   hh.BinEdges = 0.5:nDeps+1.5; 
%   title(sprintf('Departure Time Distribution [%.1f]', demandPeakness))

  figure; clf;
  % bus stop grid not being considered
  scatter(BS_XY(~BSinCH,1),BS_XY(~BSinCH,2),25,0.9*[1,1,1],'+');
  hold on
  % bus stops being considered
  scatter(BS_XY(BSinCH,1),BS_XY(BSinCH,2),25,0.1*[1,1,1],'+');
  % paseenger locations with number of departure labelled
  scatter(pax_XY(:,1),pax_XY(:,2),200,'r.');
  % plot departure number if not too many pax
  if nPax < 100
    dx = 0; dy = .02; % displacement so the text does not overlay the data points
    text(pax_XY(:,1)+dx,pax_XY(:,2)+dy, cellstr(num2str(pax_depT(:))));
  end
  % transit stations
  scatter(station_XY(:,1),station_XY(:,2),150,'filled','square','blue');
  % charger locations
  scatter(charger_XY(:,1),charger_XY(:,2),80,'filled','^','green');
  axis equal
  legend({sprintf('(BusStop) [%d]',size(BS_XY,1)-nBS),...
    sprintf('BusStop [%d]',nBS), sprintf('Pax [%d]',nPax),'Station',...
    sprintf('Charger [%d]',nCharger)});

  if nStation>1
    title(sprintf('Station Gap = %.2f. Bus Stop Gap= %.2f. \n Walk Dist = %.2f. Pax radius =  [%.2f,%.2f]',...
      Station_separation,BS_separation,maxWalkingDist,Pax_minRadius(1), Pax_maxRadius(1)))
  else
    title(sprintf('Bus Stop Gap= %.2f. Walk Dist = %.2f. Pax radius = [%.2f,%.2f]',...
      BS_separation,maxWalkingDist,Pax_minRadius(1), Pax_maxRadius(1)))
  end

  % draw circles around each passenger location
  % h_circles = viscircles(pax_XY, maxWalkingDist*ones(size(pax_XY,1),1), 'color', [1 0.9 0.9]);

end
