function fh = plotScenario(T_Passenger,T_busStop,T_Charger,T_Station,T_depot)

% needs to be based on outputs of
% [T_busFleet, T_Passenger,T_busStop,T_Charger,T_Station,allStationDeps,T_depot] = ...
%        generateScenario(saveData_UrbanMorph,nStations, Station_separation, nPassengers, Pax_minRadius, Pax_maxRadius,...
%        paxSeparation, maxWalkingDist, BS_separation, nCharger, charger_radius, demandPeakness, PLOTFLAG);

% ========================================
% ======= PLOTTING DATA IF DESIRED =======
% figure; hh = histogram(pax_depT);
% hh.BinEdges = 0.5:nDeps+1.5;
% title(sprintf('Departure Time Distribution [%.1f]', demandPeakness))


fh = figure;  % Create a visible figure
hold on
% bus stops being considered
scatter(T_busStop.busStop_X,T_busStop.busStop_Y,25,0.1*[1,1,1],'+');
nBS = height(T_busStop);
legendCell = {sprintf('Meeting Points [%d]',nBS)};

% transit stations display each one a different colour
nStation = height(T_Station);
if nStation>1, cmap = colormap_generator(nStation);
  for i = 1:nStation
    sh = scatter(T_Station.station_X(i),T_Station.station_Y(i),150,'filled','square','blue');
    sh.MarkerEdgeColor = cmap(i,:);
    sh.MarkerFaceColor = cmap(i,:);
    pax_ind = T_Passenger.passenger_StationID == i;
    ph = scatter(T_Passenger.passenger_X(pax_ind),T_Passenger.passenger_Y(pax_ind),200,'.');
    ph.MarkerEdgeColor = cmap(i,:);
    ph.MarkerFaceColor = cmap(i,:);
    legendCell = [legendCell, {'Station',sprintf('Customers [%d]',sum(pax_ind))}]; %#ok<*AGROW>
  end
else
  sh = scatter(T_Station.station_X,T_Station.station_Y,150,'filled','square','blue');
  pax_ind = T_Passenger.passenger_StationID == 1;
  ph = scatter(T_Passenger.passenger_X(pax_ind),T_Passenger.passenger_Y(pax_ind),200,'.');
  ph.MarkerEdgeColor = 'r';
  ph.MarkerFaceColor = 'r';
  
  legendCell = [legendCell, {'Station',sprintf('Customers [%d]',sum(pax_ind))}];
end

% depot location
dh = scatter(T_depot.depot_X,T_depot.depot_Y,80,'filled','o','mag');
dh.MarkerEdgeColor = 'k';
dh.MarkerFaceColor = 'm';
legendCell = [legendCell, {'Depot'}];

% charger locations
ch = scatter(T_Charger.charger_X,T_Charger.charger_Y,80,'filled','^','green');
nCharger = height(T_Charger);
legendCell = [legendCell, {sprintf('Chargers [%d]',nCharger)}];

axis equal; legend(legendCell)
xlabel('X coordinate (km)'); ylabel('Y coordinate (km)')

% if any(paxMissed)
%   disp('These customers cant reach any meeting point:');
%   disp(pax_XY(paxMissed,:))
%   scatter(T_Passenger.passenger_X(paxMissed),T_Passenger.passenger_Y(paxMissed),200,'k*');
% end


% TITLE that notes various settings
% if nStation>1
%   title(sprintf('Station Gap = %.2f. Meeting Point Separation= %.2f. Max  Walking Dist = %.2f. Customer radius = [%.2f,%.2f]',...
%     Station_separation,BS_separation,maxWalkingDist,Pax_minRadius(1), Pax_maxRadius(1)))
% else
%   title(sprintf('Meeting Point Separation= %.2f. Max  Walking Dist = %.2f. Customer radius = [%.2f,%.2f]',...
%     BS_separation,maxWalkingDist,Pax_minRadius(1), Pax_maxRadius(1)))
% end

% draw circles around each passenger location of max walking distance
% h_circles = viscircles(pax_XY, maxWalkingDist*ones(size(pax_XY,1),1), 'color', [1 0.9 0.9]);
