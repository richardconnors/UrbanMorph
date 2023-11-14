
c50u = T_Passenger.passenger_DepartureTime;

c100u = T_Passenger.passenger_DepartureTime;

c200u = T_Passenger.passenger_DepartureTime;

figure; 

h200 = histogram(c200); hold on;
h100 = histogram(c100);
h50 = histogram(c50);

h200u = histogram(c200u); h200u.FaceColor = 0.9*[1,1,1];
h100u = histogram(c100u); h100u.FaceColor = 0.6*[1,1,1];
h50u = histogram(c50u); h50u.FaceColor = 0.3*[1,1,1];

nDeps = size(allStationDeps,1);
h50.BinEdges = 0.5:nDeps+1.5;
h100.BinEdges = 0.5:nDeps+1.5;
h200.BinEdges = 0.5:nDeps+1.5;
title('Departure Time Distribution (Peak)')

xticks(1:13); % Set the x-axis ticks to represent hours
xticklabels({'06:00', '06:20', '06:40',...
  '07:00', '07:20', '07:40',...
  '08:00', '08:20', '08:40',...
  '09:00', '09:20', '09:40',...
  '10:00'}); % Set the correspon
legend({'c200op','c100op','c50op','c200p','c100p','c50p'})
ylabel('Number of Requests')
xlabel('Train Departure Time')

this_instance = 'p_op_distribution_50_100_200';
print(gcf, [saveData_Figures,this_instance,'.jpg'], '-djpeg', '-r300');
savefig(gcf, [saveData_Figures,this_instance,'.fig']);
