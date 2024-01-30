load GRIDTESTS_PROCESSED

cityS = unique(R_instance.Pax_maxRadius); nS = numel(cityS);
cityP = unique(R_instance.nPax); nP = numel(cityP);

figure('Units', 'normalized', 'OuterPosition', [0, 0, 1, 1]);
for i = nP:-1:1
  thisR = R_instance(R_instance.nPax == cityP(i),:);
  bh = boxchart(thisR.Pax_maxRadius,thisR.FleetSize);
  hold on
  bh.BoxWidth = 2; bh.MarkerStyle = 'none';

  yp = polyfit(sqrt(thisR.Pax_maxRadius),thisR.FleetSize,1)
  xFit = linspace(min(thisR.Pax_maxRadius), max(thisR.Pax_maxRadius), 100);
  yFit = polyval(yp, sqrt(xFit));
  fh = plot(xFit, yFit, 'r-', 'HandleVisibility', 'off');

end
xlabel('Max Radius of Operational Area (km)')
ylabel('Fleet Size')
title(sprintf('Box Plot of Fleet Size Variability'))
legend(arrayfun(@(n) sprintf('Customers %d', n), flipud(cityP), 'UniformOutput', false),...
  'Location','northwest');


figure;
for i = 1:nS
  thisR = R_instance(R_instance.Pax_maxRadius== cityS(i),:);
  bh = boxchart(thisR.nPax,thisR.FleetSize);
  hold on
  bh.BoxWidth = 20;
  yp = polyfit(thisR.nPax,thisR.FleetSize,1)
  xFit = linspace(min(thisR.nPax), max(thisR.nPax), 100);
  yFit = polyval(yp, xFit);
  fh = plot(xFit, yFit, 'r-', 'HandleVisibility', 'off');

end
xlabel('Number of Customers')
ylabel('Fleet Size')
title(sprintf('Box Plot of Fleet Size Variability'))
legend(arrayfun(@(n) sprintf('Max Radius %.1f', n), cityS, 'UniformOutput', false),...
  'Location','northwest');





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

fh1 = boxPlot3D(R_instance.FleetSize, R_instance.nPax, R_instance.Pax_maxRadius);
xlabel('nCustomers'); ylabel('City Radius'); zlabel('Fleet Size')
title('OPERATOR')

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
    bh = boxchart(thisR.CusDirectKms,thisR.rideTime,'Notch','on');
    xlabel('Direct Distance to Transit (km)')
    ylabel('Ride Time (minutes)')
    title(sprintf('Box Plot of Customer Ride Time Variability [#Cus = %d, KM = %03.1f]',cityP(i),cityS(j)))
    bh.BoxWidth = max(thisR.CusDirectKms)*0.1/30;
    hold on
    bh2 = boxchart(thisR.CusDirectKms,thisR.walkTime,'notch','on');
    bh2.BoxWidth = max(thisR.CusDirectKms)*0.1/30;


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



cityS = unique(R_instance.Pax_maxRadius);
cityP = unique(R_instance.nPax);
for i = 1:numel(cityP)
  figure('Units', 'normalized', 'OuterPosition', [0, 0, 1, 1]);;
  for j= 1:numel(cityS)
    thisR = R_Pax(and(R_Pax.nPax == cityP(i),R_Pax.Pax_maxRadius == cityS(j)),:);

    thisR.NormRT = p.busSpeed*thisR.rideTime./thisR.CusDirectKms;
    subplot(2,3,j)
    S = grpstats(thisR,"CusDirectKms",["mean","std"],"DataVars",["NormRT"]);
    h = scatter(S.CusDirectKms,S.std_NormRT);
    %     hold on;
    %     sh(2) = scatter(S.CusDirectKms,S.mean_rideTime);
    %     sh(3) = scatter(S.CusDirectKms,S.max_rideTime);
    %     ind = S.std_rideTime>0
    %     sh(4) = scatter(S.CusDirectKms,S.mean_rideTime-S.std_rideTime);
    %     sh(5) = scatter(S.CusDirectKms,S.mean_rideTime+S.std_rideTime);
    title(sprintf('Std Dev of Customer Ride Time Variability [#Cus = %d, KM = %03.1f]',cityP(i),cityS(j)))
    ylabel('Std Dev of individual normalised ride time')
    h.Marker = 'o';
    h.SizeData = 5;
    ax = gca; ax.YLim = [0,0.6];
    
  end
end
% ======================================================
% ======================================================
cityS = unique(R_instance.Pax_maxRadius);
cityP = unique(R_instance.nPax);
for i = 1:numel(cityP)
  figure('Units', 'normalized', 'OuterPosition', [0, 0, 1, 1]);;
  for j= 1:numel(cityS)
    thisR = R_Pax(and(R_Pax.nPax == cityP(i),R_Pax.Pax_maxRadius == cityS(j)),:);
    subplot(2,3,j)
    S = grpstats(thisR,"CusDirectKms",["min","mean","max","std"],"DataVars",["rideTime"]);
    sh(1) = scatter(S.CusDirectKms,S.std_rideTime);
    %     hold on;
    %     sh(2) = scatter(S.CusDirectKms,S.mean_rideTime);
    %     sh(3) = scatter(S.CusDirectKms,S.max_rideTime);
    %     ind = S.std_rideTime>0
    %     sh(4) = scatter(S.CusDirectKms,S.mean_rideTime-S.std_rideTime);
    %     sh(5) = scatter(S.CusDirectKms,S.mean_rideTime+S.std_rideTime);
    title(sprintf('Std Dev of Customer Ride Time Variability [#Cus = %d, KM = %03.1f]',cityP(i),cityS(j)))

    h_min.Marker = 'v';
    h_mean.Marker = 'o';
    h_max.Marker = '^';
    h_mean.SizeData = 5;
    h_min.SizeData = 5;
    h_max.SizeData = 5;

  end
end

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
