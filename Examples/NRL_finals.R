source('module.r')
set.seed(2021)

training = read.csv('Data/nrl-2021.csv')
mat = Data2Mat(ResultSplitter(training))
model = BTC(mat, 1)

premiers = c()
granny= c()
nruns = 10000


for(run in 1:nruns){
  # Week 1
  QF1 = MatchSim(model, 'Storm', 'Sea Eagles')
  QF2 = MatchSim(model, 'Panthers', 'Rabbitohs')
  EF1 = MatchSim(model, 'Roosters', 'Titans')
  EF2 = MatchSim(model, 'Eels', 'Knights')
  
  
  # Week 2
  if(QF1 == 'Storm'){
    SF1 = MatchSim(model, 'Sea Eagles', EF1)
  }
  if(QF1 == 'Sea Eagles'){
    SF1 = MatchSim(model, 'Storm', EF1)
  }
  if(QF2 == 'Panthers'){
    SF2 = MatchSim(model, 'Rabbitohs', EF2)
  }
  if(QF2 == 'Rabbitohs'){
    SF2 = MatchSim(model, 'Panthers', EF2)
  }
  
  
  # Week 3
  PF1 = MatchSim(model, QF1, SF2)
  PF2 = MatchSim(model, QF2, SF1)
  
  
  # Week 4
  if(runif(1) > 0.5){
    GF = MatchSim(model, PF1, PF2)
    premiers = cbind(premiers, GF)
    granny = rbind(granny, c(PF1, PF2, 1, 0))
  }
  else{
    GF = MatchSim(model, PF2, PF1)
    premiers = cbind(premiers, GF)
    granny = rbind(granny, c(PF2, PF1, 1, 0))
  }
  
}

# How many times does each team become premier?
table(premiers)/nruns
table(premiers)

# How many times does each pairing occur?
granny = data.frame(granny)
colnames(granny) = c('home.team', 'away.team', 'home_scores', 'away_scores')
data = Data2Mat(granny)
upper.tri(data$w) * data$w + upper.tri(data$w) * t(data$w)
(upper.tri(data$w) * data$w + upper.tri(data$w) * t(data$w)) / nruns



# Plots
ptable = table(premiers)
ordered = ptable[order(ptable)]
df = data.frame(ordered)
ggplot(df, aes(x=premiers, y=Freq)) + 
  geom_bar(stat = "identity") + 
  coord_flip()


# Comparing to BET365 odds
odds = data.frame(rbind(c('Storm',2.2), c('Panthers', 3.8), c('Rabbitohs', 5.5), c('Sea Eagles', 10), c('Roosters', 19)))
colnames(odds) = c('team', 'odds')
odds$probs = (1 / as.numeric(odds$odds)) / (sum(1 / as.numeric(odds$odds)))
prtable = ptable/nruns
odds$ests = c(prtable['Storm'], prtable['Panthers'], prtable['Rabbitohs'], prtable['Sea Eagles'], prtable['Roosters'])


ggplot(odds, aes(x = reorder(team, probs - ests), y = probs - ests)) + 
  geom_bar(stat = 'identity') +
  xlab('Team') + 
  ylab('Difference between Bet365 and BT Implied Probabilities') +
  coord_flip()
