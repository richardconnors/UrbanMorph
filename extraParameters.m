% Below are some parameters still missing:
% 	Charging speed: 0.83 kw/minutes
% 	Average walking speed: 0.085 km/minute
% 	Average bus speed: 0.83 km/minutes
% 	Number of dummies for each chargers: 3
% 	Maximum walking distance: 1.5 km
% 	Service time at bus stops: 0.5 minutes
% 	Maximum waiting time at transit stops: 15 minutes

% 	Bus type 1
% 	Maximum capacity Q_max_1 = 10 passengers
% 	Maximum battery capacity B_max_1 = 35.8 kw
% 	Minimum state-of-charge: 10% of B_max_1
% 	Maximum state-of-charge: 80% of B_max_1
% 	Energy consumption rate: 0.24 kw/km
% 	Number of bus type 1: 
% 	((number_of_customers/2)*0.7)/Q_max_1 => round up the number
% 	Half customers are served by this bus type with the maximum ridership 70%
% 	Initial battery level: 
% 	If there are n buses of this type, the first one is 20% of the capacity, the second one is 30%.... but no more than 80% of the maximum
% 	So for the n^th bus within this type, the initial battery level is: B_max_2×min(0.2 + 0.1(n-1),  0.8)

% 	Bus type 2
% 	Maximum capacity Q_max_2 = 20 passengers
% 	Maximum battery capacity B_max_2 = 53.7 kw
% 	Minimum state-of-charge: 10% of B_max_2
% 	Maximum state-of-charge: 80% of B_max_2
% 	Energy consumption rate: 0.29 kw/km
% 	Number of bus type 1: ((number_of_customers/2)*0.7)/Q_max_2 => round up the number
% 	Initial battery level: B_max_2×min(0.2 + 0.1(n-1),  0.8)

% It would be nice if you could add somewhere the 
% number of passengers 
% number of transit stops
% number of bus stops 
% number of chargers
% otherwise the numbers can be gotten from the length of the lines in
% .txt/.csv file but sometimes the last line might be empty. Probably there
% are other ways to retrieve these numbers, let me know if you have a
% better solution! Also, can you please change the coordinates  unit to km
% just to be consistent with the code.
