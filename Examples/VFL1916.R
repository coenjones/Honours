## 1916 VFL Season example - for illustrating properties of the package
#
## Fixtures and Results sourced from Wikipedia - Further sourced from 
## Rogers, S. & Brown, A., Every Game Ever Played: VFL/AFL Results 1897–1997 (Sixth Edition), 
## Viking Books, (Ringwood), 1998. ISBN 0-670-90809-6
source('module.r')

data = read.csv('Data/VFL1916.csv')
mat = Data2Mat(data)

model = BTC(mat, 0)


# Showing the stability of our stochastic optimisation.
thetas = model$Theta
ses = model$Std.Errors

for(k in 1:99){
  thetas = cbind(thetas, BTC(mat, 0)$Theta)
  ses = cbind(ses, BTC(mat, 0)$Std.Errors)
}

max(thetas['Carlton',]) - min(thetas['Carlton',])
max(thetas['Richmond',]) - min(thetas['Richmond',])
max(thetas['Collingwood',]) - min(thetas['Collingwood',])

apply(thetas[,], 1, sd)
