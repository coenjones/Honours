source('module.r')
nfl16 = read.csv('Data/nfl-2016.csv')
nfl17 = read.csv('Data/nfl-2017.csv')
nfl18 = read.csv('Data/nfl-2018.csv')
nfl19 = read.csv('Data/nfl-2019.csv')
data = rbind(nfl16, nfl17, nfl18, nfl19)
dat = ResultSplitter(data)
mat = Data2Mat(dat)

unique = BTC(mat, mod = 4)

BT_plot(unique)
