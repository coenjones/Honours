source('module.r')  # Loading the implementation

AFL = read.csv('Data/AFL0919.csv')
head(AFL)

data = Data2Mat(AFL)  # Converting list of results into Win and Loss matrices
WLTable(data)
data$w[1:4, 1:4]  # Win matrix


vanilla = BTC(data, mod = 0, ref = 'last')
vanilla$Theta
vanilla$Std.Errors
BT_plot(vanilla)


teamspec = BTC(data, mod = 3)
teamspec$Theta
teamspec$Alpha
teamspec$Theta.Errs
teamspec$Alpha.Errs
BT_plot(teamspec)
BT_test(teamspec, vanilla, data)
BT_test(teamspec, BTC(data, 1), data) # Fit a common-homeground and compare



# Creating a Victoria vs Interstate Hierarchy
Vics = c("Richmond", "Western Bulldogs", "Geelong", "North Melbourne", 
         "Collingwood", "Hawthorn", "Carlton", "Melbourne", "St Kilda", 
         "Essendon")
Inters = c("Port Adelaide", "West Coast", "Gold Coast", "GWS Giants",
           "Fremantle", "Adelaide", "Brisbane", "Sydney")
Groups = rbind(Vics, Inters)
R = SimpleHierarchy(Groups, data)
R[1:5, 1:5]

minih = BTC(data, mod = 2, rel = R)
minih$Alpha
minih$Std.Errors
BT_test(minih, vanilla, data, rel = R)
BT_test(minih, teamspec, data, rel = R)
BT_test(minih, BTC(data, 1), data, rel = R)



# Some other nice things
BT_bootstrap(minih, AFL, runs = 5, rel = R)
BT_predict(minih, 'Brisbane', 'Sydney', rel = R)
for(k in 1:30){print(MatchSim(minih, 'Brisbane', 'Sydney', rel = R))}
head(SeasonSim(minih, AFL, rel = R))
