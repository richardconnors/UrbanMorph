load GRIDTESTS_LOADED

% Total veh kms
figure;
h1 = scatter3(R_instance.nPax, R_instance.Pax_maxRadius, R_instance.VehEmptyKms);
hold on
xlabel('nCustomers'); ylabel('City Radius'); zlabel('Empty VehKms')
% h2 = scatter3(R_instance.nPax, R_instance.Pax_maxRadius, R_instance.VehKms./R_instance.CusDirectKms);
% legend({'Total Direct Distance','Total Tour Distance'})
cmap=colormap_generator(2);
h1.SizeData = 40;
h1.MarkerFaceColor = cmap(1,:);
h1.MarkerFaceAlpha = 0.7;
h2.SizeData = 40;
h2.MarkerFaceColor = cmap(2,:);
h2.MarkerFaceAlpha = 0.7;

% fleet size
figure;
h1 = scatter3(R_instance.nPax, R_instance.Pax_maxRadius, R_instance.FleetSize);
xlabel('nCustomers'); ylabel('City Radius'); zlabel('fleet size')
cmap=colormap_generator(2);
h1.SizeData = 40;
h1.MarkerFaceColor = cmap(1,:);
h1.MarkerFaceAlpha = 0.7;

fh1 = boxPlot3D(R_bus.MaxOcc, R_bus.nPax, R_bus.Pax_maxRadius);
xlabel('nCustomers'); ylabel('City Radius'); zlabel('Occupancy')

fh2 = boxPlot3D(R_bus.CusTravelledKms./R_bus.VehKms, R_bus.nPax, R_bus.Pax_maxRadius);
xlabel('nCustomers'); ylabel('City Radius'); zlabel('Bus Utilization Ratio')

% chargeable kms vs bus kms
fh3 = boxPlot3D(R_bus.CusDirectKms./R_bus.VehKms, R_bus.nPax, R_bus.Pax_maxRadius);
xlabel('nCustomers'); ylabel('City Radius'); zlabel('Chargeable Ratio')


% ======================================================
% ======================================================
% individual passenger experience variability
cityS = unique(R_instance.Pax_maxRadius);
cityP = unique(R_instance.nPax);
for i = 1:numel(cityP)
  figure('Units', 'normalized', 'OuterPosition', [0, 0, 1, 1]);;
  for j= 1:numel(cityS)
    thisR = R_Pax(and(R_Pax.nPax == cityP(i),R_Pax.Pax_maxRadius == cityS(j)),:);

    subplot(2,3,j)
    bh = boxchart(thisR.CusDirectKms,thisR.rideTime);
    xlabel('Direct Distance to Transit (km)')
    ylabel('Ride Time (minutes)')
    title(sprintf('Box Plot of Customer Ride Time Variability [#Cus = %d, KM = %03.1f]',cityP(i),cityS(j)))
    bh.BoxWidth = max(thisR.CusDirectKms)*0.2/30;
    hold on
    bh2 = boxchart(thisR.CusDirectKms,thisR.walkTime);
    bh2.BoxWidth = max(thisR.CusDirectKms)*0.2/30;
    legend({'Ride Time','Walk Time'},'location','northwest')
    axis tight;
    
%     thisR2 = thisR(:,["paxID","rideTime","walkTime","CusDirectKms"]);
%     stdDevData = grpstats(thisR,"paxID",["std","mean"]);
%     figure; scatter(stdDevData.mean_CusDirectKms,stdDevData.std_rideTime)
  end
end



% ======================================================
% ======================================================

%   % create matrix with a column for each individual ordered by direct
%   % dist
%   % assume nPax is true
%   [thisPax, idx] = unique(thisR.paxID); % remember more Pax than cityP
%   % what order is by distance
%   thisDD = thisR.CusDirectKms(idx);
%   [~,ddidx] = sort(thisDD,'ascend');
%   paxOrder = thisPax(ddidx);
%   dd = thisR.CusDirectKms(idx(ddidx)); % shows that the ordering is ok
%   nRuns = max(thisR.id) - min(thisR.id) + 1;
%   rideMatrix = NaN(nRuns,numel(thisPax));
%
%   for kk = 1:numel(thisPax) % for each person
%     thisRides = thisR.rideTime(thisR.paxID == paxOrder(kk));
%     rideMatrix(1:numel(thisRides),kk) = thisRides;
%   end
