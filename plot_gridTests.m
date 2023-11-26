% load GRIDTESTS_PROCESSED
load FIXED80_PROCESSED

cityS = unique(R_instance.Pax_maxRadius);
cityP = unique(R_instance.nPax);
nS = numel(cityS);
nP = numel(cityP);


% ============== BUS LEVEL BOXPLOTS ==============
% how does fleet size vary with population density?
figure; legendtext = [];
for i = 1:nP
  thisR = R_bus(R_bus.nPax == cityP(i),:);
  thisArea = (pi*thisR.Pax_maxRadius.^2) - (pi*p.Pax_minRadius.^2);
  bh = boxchart(thisR.nPax./thisArea,thisR.CusDirectKms./thisR.VehKms);
  bh.BoxWidth = 0.1;

  legendtext = [legendtext;sprintf('nPax = %03d',cityP(i))]; %#ok<*AGROW>
  hold on;
end
xlabel('Population Density [p/km2]'); ylabel('Chargeable Ratio');
legend(legendtext)

% ============== INSTANCE LEVEL BOXPLOTS ==============
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

% chargeable kms vs bus kms
fh5 = boxPlot3D(R_instance.FleetSize./R_instance.nPax, R_instance.nPax, R_instance.Pax_maxRadius);
xlabel('nCustomers'); ylabel('City Radius'); zlabel('Fleet Size Per Person')
title('OPERATOR')

% how does fleet size vary with population density?
figure; legendtext = [];
for i = 1:nP
  thisR = R_instance(R_instance.nPax == cityP(i),:);
  thisArea = (pi*thisR.Pax_maxRadius.^2) - (pi*p.Pax_minRadius.^2);
  bh = boxchart(thisR.nPax./thisArea,thisR.FleetSize./thisR.nPax);
  bh.BoxWidth = 0.1;
  legendtext = [legendtext;sprintf('nPax = %03d',cityP(i))];
  hold on;
end
xlabel('Population Density [p/km2]'); ylabel('Buses/Pax'); title('Fleet Size Per Person')
legend(legendtext)


% ======================================================
% ============== PERSON LEVEL BOXPLOTS =================
% ======================================================
% individual passenger experience variability - boxplot for each pop/radius
% every passenger contributes on vertical boxplot element
for i = 1:numel(cityP)
  figure('Units', 'normalized', 'OuterPosition', [0, 0, 1, 1]);;
  for j= 1:numel(cityS)
    thisR = R_Pax(and(R_Pax.nPax == cityP(i),R_Pax.Pax_maxRadius == cityS(j)),:);

    subplot(2,3,j)
    bh = boxchart(thisR.CusDirectKms,thisR.rideTime);
    xlabel('Direct Distance to Transit (km)')
    ylabel('Ride Time (minutes)')
    title(sprintf('Box Plot of Customer Ride Time Variability [#Cus = %d, KM = %03.1f]',cityP(i),cityS(j)))
    bh.BoxWidth = max(thisR.CusDirectKms)*0.1/30;
    hold on
    bh2 = boxchart(thisR.CusDirectKms,thisR.walkTime);
    bh2.BoxWidth = max(thisR.CusDirectKms)*0.1/30;

    % plot the 0.83 km/min line and the 1.5 multiplier
    uh = line([0,cityS(j)],[0,1.5*cityS(j)./p.busSpeed]);
    uh.LineWidth = 2; uh.Color = 'k'; uh.LineStyle = ":";
    lh = line([0,cityS(j)],[0,cityS(j)./p.busSpeed]);
    lh.LineWidth = 2; lh.Color = 'k'; lh.LineStyle = "--";
    ax = gca; ax.XLim = [0,cityS(j)];
    ax.YLim = [0,max(max(thisR.walkTime),1.2*1.5*cityS(j)./p.busSpeed)];

    legend({'Ride Time','Walk Time','Max Rerouting','Direct Service'},'location','northwest')
    %     thisR2 = thisR(:,["paxID","rideTime","walkTime","CusDirectKms"]);
    %     stdDevData = grpstats(thisR,"paxID",["std","mean"]);
    %     figure; scatter(stdDevData.mean_CusDirectKms,stdDevData.std_rideTime)
  end
end

% ======================================================
% ======================================================
% How much of max possible variability does each person experience?
cmap = colormap_generator(10); 
for i = 1:numel(cityP)
  figure('Units', 'normalized', 'OuterPosition', [0,0,0.95,0.95]);
  for j= 1:numel(cityS)
    subplot(2,3,j)
    thisR = R_Pax(and(R_Pax.nPax == cityP(i),R_Pax.Pax_maxRadius == cityS(j)),:);
    % for each person calculate their min and max ride time.
    % if detour factor is 1.5 then normalise their experience as values
    % [1.0, 1.5] there may be some values just outside this range.
    minRT = thisR.CusDirectKms./p.busSpeed;
    maxRT = 1.5*thisR.CusDirectKms./p.busSpeed;
    nzdRT = p.busSpeed*thisR.rideTime./thisR.CusDirectKms;
    thisR = addvars(thisR,nzdRT, 'NewVariableNames', 'NzdRideTime');

    % want to color the box chart by the fixed and variable customers
    % G = groupsummary(thisR, "paxID",["min","max","std"],"rideTime");
    G = groupsummary(thisR, "paxID");
    fi = G.paxID(G.GroupCount==max(G.GroupCount));
    
    yh = yline([1,1.5]); [yh(:).LineWidth] = deal(2); [yh(:).HandleVisibility] = deal({"off"});
    yh(1).LineStyle = "--"; yh(2).LineStyle = ":";
    hold on
    title(sprintf('Customer Ride Time Variability [#Cus = %d, KM = %03.1f]',cityP(i),cityS(j)))

    boxW = 0.1*cityS(j)/30;
    bf = boxchart(thisR.CusDirectKms(ismember(thisR.paxID,fi)),thisR.NzdRideTime(ismember(thisR.paxID,fi)));
    bf.BoxWidth = boxW; bf.MarkerStyle = 'none'; bf.BoxFaceColor = cmap(1,:);
    bu = boxchart(thisR.CusDirectKms(~ismember(thisR.paxID,fi)),thisR.NzdRideTime(~ismember(thisR.paxID,fi)));
    bu.BoxWidth = boxW; bu.MarkerStyle = 'none'; bu.BoxFaceColor = cmap(3,:);
    ax = gca; ax.XLim = [1.5,cityS(j)]; ax.YLim = [0.6,1.9];
    legend({"Regular","Random"},'location','northwest')
  end
end






cityS = unique(R_instance.Pax_maxRadius);
cityP = unique(R_instance.nPax);
for i = 1:numel(cityP)
  figure('Units', 'normalized', 'OuterPosition', [0, 0, 1, 1]);;
  for j= 1:numel(cityS)
    thisR = R_Pax(and(R_Pax.nPax == cityP(i),R_Pax.Pax_maxRadius == cityS(j)),:);
    subplot(2,3,j)
    S = grpstats(thisR,"CusDirectKms",["min","mean","max","std"],"DataVars",["rideTime"]);
    sh(1) = scatter(S.CusDirectKms,S.min_rideTime); hold on;
    sh(2) = scatter(S.CusDirectKms,S.mean_rideTime);
    sh(3) = scatter(S.CusDirectKms,S.max_rideTime);
    %     ind = S.std_rideTime>0
    sh(4) = scatter(S.CusDirectKms,S.mean_rideTime-S.std_rideTime);
    sh(5) = scatter(S.CusDirectKms,S.mean_rideTime+S.std_rideTime);
    h_min.Marker = 'v';
    h_mean.Marker = 'o';
    h_max.Marker = '^';
    h_mean.SizeData = 5;
    h_min.SizeData = 5;
    h_max.SizeData = 5;

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
