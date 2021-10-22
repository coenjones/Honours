source('module.R')
data = read.csv('Data/NBA16-20.csv')
mat = Data2Mat(data)
R = NBA_relat()


vanilla = BTC(mat, mod = 0)
common = BTC(mat, mod = 1)
minih = BTC(mat, mod = 2, rel = R)
single = BTC(mat, mod = 4)
hierarchical = BTC(mat, mod = 5, rel = R)


BT_plot(vanilla)
BT_plot(common)
BT_plot(minih)
BT_plot(single)
BT_plot(hierarchical)


BT_summary(vanilla)
BT_summary(common)
BT_summary(minih)
BT_summary(single)
BT_summary(hierarchical)


BT_test(vanilla, common, mat)
BT_test(common, minih, mat, R)
BT_test(minih, single, mat, R)
BT_test(single, hierarchical, mat, R)



BT_predict(vanilla, 'New York Knicks', 'Brooklyn Nets')
BT_predict(common, 'New York Knicks', 'Brooklyn Nets')
BT_predict(minih, 'New York Knicks', 'Brooklyn Nets', R)
BT_predict(single, 'New York Knicks', 'Brooklyn Nets')
BT_predict(hierarchical, 'New York Knicks', 'Brooklyn Nets', R)


theta = H2Vec(hierarchical)[1:30]
alpha = H2Vec(hierarchical)[31:120]

HH = Hess5(theta, alpha, mat, R)
isSymmetric(HH)
se = Inverse(HH)
