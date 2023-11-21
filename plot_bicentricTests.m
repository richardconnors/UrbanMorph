
load BICENTRIC_PROCESSED

SS = unique(R_instance.sSeparation);
i5 = R_instance.Pax_maxRadius == 5;
i15 = R_instance.Pax_maxRadius == 15;

figure; 
boxchart(R_instance.sSeparation(i5),R_instance.FleetSize(i5)); hold on
boxchart(R_instance.sSeparation(i15),R_instance.FleetSize(i15)); hold on
xlabel('Station Separation (km)'); ylabel('Fleet Size'); title('OPERATOR')
legend({'Radius = 5km','Radius = 15km'})
ax = gca; ax.YLim(1) = 0;

fh2 = boxPlot3D(R_bus.MaxOcc, R_bus.nPax, R_bus.Pax_maxRadius);
xlabel('nCustomers'); ylabel('City Radius'); zlabel('Maximum Occupancy')
title('OPERATOR')

fh3 = boxPlot3D(R_bus.CusTravelledKms./R_bus.VehKms, R_bus.nPax, R_bus.Pax_maxRadius);
xlabel('nCustomers'); ylabel('City Radius'); zlabel('Bus Utilization Ratio')
title('OPERATOR')

% chargeable kms vs bus kms
fh4 = boxPlot3D(R_bus.CusDirectKms./R_bus.VehKms, R_bus.nPax, R_bus.Pax_maxRadius);
xlabel('nCustomers'); ylabel('City Radius'); zlabel('Chargeable Ratio')
title('OPERATOR')


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
    % plot the 0.83 km/min line and the 1.5 multiplier
    uh = line([0,cityS(j)],[0,1.5*cityS(j)./p.busSpeed]);
    uh.LineWidth = 2; uh.Color = 'k'; uh.LineStyle = ":";
    lh = line([0,cityS(j)],[0,cityS(j)./p.busSpeed]);
    lh.LineWidth = 2; lh.Color = 'k'; lh.LineStyle = "--";
    ax = gca; ax.XLim = [0,cityS(j)]; ax.YLim = [0,max(prctile(thisR.walkTime,95),1.2*1.5*cityS(j)./p.busSpeed)];
    
    legend({'Ride Time','Walk Time','Max Rerouting','Direct Service'},'location','northwest')
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
