source('module.r')
set.seed(2021)

training = read.csv('Data/afl-2021-UTC.csv')
mat = Data2Mat(ResultSplitter(training[1:198,]))
model = BTC(mat, 1)

premiers = c()
granny= c()
nruns = 10000


for(run in 1:nruns){
  # Week 1
  QF1 = MatchSim(model, 'Melbourne', 'Brisbane Lions')
  QF2 = MatchSim(model, 'Port Adelaide', 'Geelong Cats')
  EF1 = MatchSim(model, 'Western Bulldogs', 'Essendon')
  EF2 = MatchSim(model, 'Sydney Swans', 'GWS Giants')
  
  
  # Week 2
  if(QF1 == 'Melbourne'){
    SF1 = MatchSim(model, 'Brisbane Lions', EF1)
  }
  if(QF1 == 'Brisbane Lions'){
    SF1 = MatchSim(model, 'Melbourne', EF1)
  }
  if(QF2 == 'Port Adelaide'){
    SF2 = MatchSim(model, 'Geelong Cats', EF2)
  }
  if(QF2 == 'Geelong Cats'){
    SF2 = MatchSim(model, 'Port Adelaide', EF2)
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



# Comparing to Sportsbet odds
odds = data.frame(rbind(c('Port Adelaide',3.5), c('Melbourne', 3.75), c('Brisbane Lions', 5.5), 
                        c('Geelong Cats', 7), c('Western Bulldogs', 10), c('Sydney Swans', 15),
                        c('GWS Giants', 21), c('Essendon', 26)))
colnames(odds) = c('team', 'odds')
odds$probs = (1 / as.numeric(odds$odds)) / (sum(1 / as.numeric(odds$odds)))
prtable = ptable/nruns
odds$ests = c(prtable['Port Adelaide'], prtable['Melbourne'], prtable['Brisbane Lions'], 
              prtable['Geelong Cats'], prtable['Western Bulldogs'], prtable['Sydney Swans'],
              prtable['GWS Giants'], prtable['Essendon'])


ggplot(odds, aes(x = reorder(team, probs - ests), y = probs - ests)) + 
  geom_bar(stat = 'identity') +
  xlab('Team') + 
  ylab('Difference between Sportsbet and BT Implied Probabilities') +
  coord_flip()
