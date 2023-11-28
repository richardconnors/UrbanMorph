load GRIDTESTS_PROCESSED
cityS = unique(R_instance.Pax_maxRadius);
cityP = unique(R_instance.nPax);
nS = numel(cityS);
nP = numel(cityP);

i = 4; j = 5;

% ======================================================
% ============== PERSON LEVEL BOXPLOTS =================
% ======================================================
% individual passenger experience variability - boxplot for each pop/radius
% every passenger contributes on vertical boxplot element
figure('Units', 'normalized', 'OuterPosition', [0, 0, 1, 1]);;
thisR = R_Pax(and(R_Pax.nPax == cityP(i),R_Pax.Pax_maxRadius == cityS(j)),:);

thisPax = 5;
figure;
boxchart(thisR.CusDirectKms(thisR.paxID==thisPax),thisR.rideTime(thisR.paxID==thisPax))
ax = gca;   ax.YLim(1) = 0;
xlabel('Direct Distance to Transit (km)')
ylabel('Ride Time (minutes)')
title(sprintf('Box Plot of Customer Ride Time Variability [CusID = %d]',thisPax))


bh = boxchart(thisR.CusDirectKms,thisR.rideTime);
xlabel('Direct Distance to Transit (km)')
ylabel('Ride Time (minutes)')
title(sprintf('Box Plot of Customer Ride Time Variability [#Cus = %d, KM = %03.1f]',cityP(i),cityS(j)))
bh.BoxWidth = max(thisR.CusDirectKms)*0.1/30;
hold on
bh2 = boxchart(thisR.CusDirectKms,thisR.walkTime);
bh2.BoxWidth = max(thisR.CusDirectKms)*0.1/30;

% % plot the 0.83 km/min line and the 1.5 multiplier
% uh = line([0,cityS(j)],[0,1.5*cityS(j)./p.busSpeed]);
% uh.LineWidth = 2; uh.Color = 'k'; uh.LineStyle = ":";
% lh = line([0,cityS(j)],[0,cityS(j)./p.busSpeed]);
% lh.LineWidth = 2; lh.Color = 'k'; lh.LineStyle = "--";
ax = gca; ax.XLim = [0,cityS(j)];
ax.YLim = [0,max(max(thisR.walkTime),1.2*1.5*cityS(j)./p.busSpeed)];

legend({'Ride Time','Walk Time','Max Rerouting','Direct Service'},'location','northwest')
%     thisR2 = thisR(:,["paxID","rideTime","walkTime","CusDirectKms"]);
%     stdDevData = grpstats(thisR,"paxID",["std","mean"]);
%     figure; scatter(stdDevData.mean_CusDirectKms,stdDevData.std_rideTime)


% ======================================================
% ============== PERSON LEVEL BOXPLOTS =================
% ======================================================
% individual passenger experience variability - boxplot for each pop/radius
% every passenger contributes on vertical boxplot element
for i = 1:numel(cityP)
  figure('Units', 'normalized', 'OuterPosition', [0, 0, .95, .95]);
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
