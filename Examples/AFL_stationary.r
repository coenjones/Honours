source('module.R')
data = read.csv('Data/AFL0919.csv')
data = data[1:(8*198 - 1),]
mat = Data2Mat(data)

overall = BTC(mat, mod = 1, ref = 'Gold Coast')
theta_avg = overall$Theta - mean(overall$Theta[,])



Rolling = function(data, season, model_type, window_size = 1, rel = NULL, ref = 'last', rollref = NULL){
  # Let season either be a vector of indices for the season, or a single value which dictates how many games are in a season
  nteams = length(unique(data$home.team))
  if(length(season) == 1){
    season = seq(1, dim(data)[1], season)
  }
  
  theta_matrix = matrix(rep(0, nteams*(length(season) - window_size + 1)), nrow = nteams)
  terr_matrix = matrix(rep(0, (nteams - 1)*(length(season) - window_size + 1)), nrow = nteams - 1)

  if(model_type > 0){
    alpha_matrix = matrix(rep(0, (nparams(Data2Mat(data), model_type, rel) - nteams)*(length(season) - window_size + 1)),
                          nrow = nparams(Data2Mat(data), model_type, rel) - nteams)
  }
  row.names(theta_matrix) = unique(data$home.team)
  row.names(terr_matrix) = unique(data$home.team)[-which(unique(data$home.team) == ref)]

  for(k in 1:(length(season) - window_size)){
    new_mat = Data2Mat(data[season[k]:(season[k + window_size] - 1), ])
    new_model = BTC(new_mat, model_type, rel = rel, ref = ref)
    for(team in new_mat$teams){
      theta_matrix[team, k] = new_model$Theta[team,]
      if(team != ref){
        terr_matrix[team, k] = new_model$Std.Errors[team,]
      }
    }
    alpha_matrix[1, k] = new_model$Alpha
  }
  
  # Handling final season, as we cannot guarantee it's length
  new_mat = Data2Mat(data[season[length(season) - window_size + 1]:dim(data)[1], ])
  new_model = BTC(new_mat, model_type, rel = rel, ref = ref)
  for(team in mat$teams){
    theta_matrix[team, length(season) - window_size + 1] = new_model$Theta[team,]
    if(team != ref){
      terr_matrix[team, length(season) - window_size + 1] = new_model$Std.Errors[team,]
    }
  }
  alpha_matrix[1, length(season) - window_size + 1] = new_model$Alpha
  
  # Handling the different ways of rebasing the rolls
  if(rollref == 'mean'){
    for(k in 1:dim(theta_matrix)[2]){
      theta_matrix[,k] = theta_matrix[,k] - mean(theta_matrix[,k])
    }
  }
  
  return(list(thetas = theta_matrix, alphas = alpha_matrix,
              terrs = terr_matrix))
}


afl_seasons = c(1, 199, 397, 595, 793, 990, 1188, 1386)

output = Rolling(data, afl_seasons, 1, window_size = 1, ref = 'Gold Coast', rollref = 'mean')
ages = read.csv('Data/AFL_average_age.csv')


team = 'Hawthorn'
plotting = data.frame(
  seasons = c(2019, 2018, 2017, 2016, 2015, 2014, 2013, 2012),
  estimates = output$thetas[team,],
  errs = 1.65*output$terrs[team,],
  overall = theta_avg[team,],
  ov_err = 1.65*overall$Std.Errors[team,],
  av_age = t(ages[ages$Team == team, -1])
)
colnames(plotting) = c('seasons','estimates','errs','overall','ov_err','av_age')
rownames(plotting) = NULL

# Plotting estimates for the windows
ggplot(plotting) + 
  geom_errorbar(aes(x = seasons, ymin = estimates - errs, ymax = estimates + errs), width = 0.4, colour = 'orange', alpha = 0.9, size = 0.9) + 
  geom_point(aes(x = seasons, y = estimates), size = 3, fill='skyblue', alpha = 0.7, shape = 23) + 
  geom_line(aes(x = seasons, y = overall), colour = 'red') + 
  geom_ribbon(aes(x = seasons, y = overall, ymax = overall + ov_err, ymin = overall - ov_err), 
              alpha=0.2, fill = 'pink')
  

# Plotting estimates versus average age
xyplot(estimates ~ av_age, data = plotting, type = c('p','r'), main = team)




# Plotting with ALL data
ests = output$thetas
age_and_est = data.frame(
  estimates = c(t(ests[order(row.names(ests)),])),
  age = c(t(data.matrix(ages[,-1])))
)
xyplot(estimates ~ age, data = age_and_est, type = c('p'), main = 'AFL')
