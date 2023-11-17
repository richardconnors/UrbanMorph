function name = scenarioName(params)


% params.maxWalkingDist = round(maxWalkingDist,5);
% params.walkingSpeed = 0.085; % km/minute
% params.busSpeed = 0.83; % km/minute
% params.busStopServiceTime = 0.5; % minutes
% params.maxTransitWaitingTime = 15; % minutes
% params.nPax = nPassengers;
% params.Pax_minRadius = Pax_minRadius;
% params.Pax_maxRadius = Pax_maxRadius;
% params.paxSeparation = paxSeparation;
% params.BS_separation = BS_separation;
% params.nStation = nStations;
% params.nCharger = nCharger;
% params.chargerRadius = charger_radius;
% params.demandPeakness = demandPeakness;

name = sprintf('P%03d_S%d_R%04.1f_W%03.1f_MP%04.1f_SS%04.1f_DP%03.1f',params.nPax,params.nStation,params.Pax_maxRadius,...
  params.maxWalkingDist,params.BS_separation,params.stationSeparation, params.demandPeakness);