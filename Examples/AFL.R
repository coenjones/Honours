source('module.R')
data = read.csv('Data/AFL0919.csv')
mat = Data2Mat(data)


Vics = c("Richmond", "Western Bulldogs", "Geelong", "North Melbourne", 
         "Collingwood", "Hawthorn", "Carlton", "Melbourne", "St Kilda", 
         "Essendon")
Inters = c("Port Adelaide", "West Coast", "Gold Coast", "GWS Giants",
           "Fremantle", "Adelaide", "Brisbane", "Sydney")
Groups = rbind(Vics, Inters)

R = SimpleHierarchy(Groups, mat)

vanilla = BTC(mat, mod = 0)
common = BTC(mat, mod = 1)
minih = BTC(mat, mod = 2, rel = R)
single = BTC(mat, mod = 4)
hierarchical = BTC(mat, mod = 5, rel = R)


BT_plot(vanilla)
BT_plot(single)


BT_summary(vanilla)
BT_summary(common)
BT_summary(minih)
BT_summary(single)
BT_summary(hierarchical)


WLTable(mat)


TOspec = BTC(mat, mod = 6)


BT_test(vanilla, common, mat)
BT_test(common, minih, mat, rel = R)
BT_test(minih, single, mat, rel = R)
BT_test(common, single, mat)
BT_test(single, TOspec, mat)
BT_test(single, hierarchical, mat, rel = R)



theta = H2Vec(hierarchical)[1:18]
alpha = H2Vec(hierarchical)[19:54]
HH = Hess5(theta, alpha, mat, R)



## Number of Games between teams
tot_hom = mat$h + t(mat$h)

hgames = c()
for(t1 in 1:18){
  for(t2 in 1:18){
    if(t1 < t2){
      new_row = c(mat$teams[t1], mat$teams[t2], tot_hom[t1,t2])
      hgames = rbind(hgames, new_row)
    }
  }
}
hgames = as.data.frame(hgames)
rownames(hgames) = NULL
hgames[order(hgames$V3),]