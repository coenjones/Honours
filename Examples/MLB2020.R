## Example using data in the 'xx - yy' form, such as from fixturedownload.com

source('module.R')
data = read.csv('Data/mlb-2020-UTC.csv')
data = ResultSplitter(data)

mat = Data2Mat(data)
model = BTC(mat, 0)
BT_plot(model)
