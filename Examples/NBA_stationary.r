source('module.R')
data = read.csv('Data/NBA16-20.csv')
mat = Data2Mat(data)

overall = BTC(mat, mod = 1)
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
    mat = Data2Mat(data[season[k]:(season[k + window_size] - 1), ])
    model = BTC(mat, model_type, rel = rel, ref = ref)
    for(team in mat$teams){
      theta_matrix[team, k] = model$Theta[team,]
      if(team != ref){
        terr_matrix[team, k] = model$Std.Errors[team,]
      }
    }
    alpha_matrix[1, k] = model$Alpha
  }
  
  # Handling final season, as we cannot guarantee it's length
  mat = Data2Mat(data[season[length(season) - window_size + 1]:(dim(data)[1]), ])
  model = BTC(mat, model_type, rel = rel, ref = ref)
  for(team in mat$teams){
    theta_matrix[team, length(season) - window_size + 1] = model$Theta[team,]
    if(team != ref){
      terr_matrix[team, length(season) - window_size + 1] = model$Std.Errors[team,]
    }
  }
  alpha_matrix[1, length(season) - window_size + 1] = model$Alpha
  
  # Handling the different ways of rebasing the rolls
  if(rollref == 'mean'){
    for(k in 1:dim(theta_matrix)[2]){
      theta_matrix[,k] = theta_matrix[,k] - mean(theta_matrix[,k])
    }
  }
  
  return(list(thetas = theta_matrix, alphas = alpha_matrix,
              terrs = terr_matrix))
}


nba_seasons = c(1, 1310, 2622, 3934)

output = Rolling(data, nba_seasons, 1, window_size = 1, ref = 'New York Knicks', rollref = 'mean')

team = 'Milwaukee Bucks'
plotting = data.frame(
  seasons = c(17, 18, 19, 20),
  estimates = output$thetas[team,],
  errs = 1.65*output$terrs[team,],
  overall = theta_avg[team,],
  ov_err = 1.65*overall$Std.Errors[team,]
)

ggplot(plotting) + 
  geom_errorbar(aes(x = seasons, ymin = estimates - errs, ymax = estimates + errs), width = 0.4, colour = 'orange', alpha = 0.9, size = 0.9) + 
  geom_point(aes(x = seasons, y = estimates), size = 3, fill='skyblue', alpha = 0.7, shape = 23) + 
  geom_line(aes(x = seasons, y = overall), colour = 'red') + 
  geom_ribbon(aes(x = seasons, y = overall, ymax = overall + ov_err, ymin = overall - ov_err), 
              alpha=0.2, fill = 'pink')
