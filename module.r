source("functions.R")
library(ggplot2)
library(matlib)


BTC = function(data, mod, rel = NULL, ref = 'last'){
  # Coen's Bradley Terry Model 
  # Returns the 'thetas' with the first parameter set to 0
  
  data$h = data$w + data$l
  data$teams = colnames(data$w)
  
  teams = data$teams
  if(mod == 0){
    theta = runif(length(teams)) # Initial guess for optimisation
    model = optim(par = theta, fn = LL, gr = LLdash, mat = data, method = 'BFGS') # Fitting model
    if(ref == 'last'){
      ref = data$teams[which.min(model$par)]
    }
    Theta = model$par - model$par[which(data$teams == ref)] # Rebasing so that min = 0
    
    # Excluding the min and finding standard errors
    HH = Hess(Theta, data)
    HH = HH[-which(data$teams == ref), -which(data$teams == ref)] # Removing fixed theta
    se = Inverse(HH)
    
    # Taking only the diagonals
    d = col(se) - row(se)
    errors = sqrt(split(se, d)$'0')
    
    return(list(Theta = data.frame(Theta, row.names = teams),
                Std.Errors = data.frame(errors, row.names = data$teams[-which(data$teams == ref)]),
                teams = teams, type = mod, ref = ref))
  }
  if(mod == 1){
    params = runif(length(teams) + 1)
    model = optim(par = params, fn = LL1, gr = LL1dash, mat = data, method = 'BFGS')
    
    if(ref == 'last'){
      ref = teams[which.min(model$par[1:length(teams)])]
    }    
    
    theta = model$par[1:length(teams)] - model$par[which(teams == ref)]
    alpha = model$par[(1 + length(teams))]
    
    # Standard Errors
    HH = Hess1(theta, alpha, data)
    HH = HH[-which(teams == ref), -which(teams == ref)]
    se = Inverse(HH)
    d = col(se) - row(se)
    errors = sqrt(split(se, d)$'0')
    
    return(list(Theta = data.frame(Theta = theta, row.names = teams),
                Alpha = alpha,
                Std.Errors = data.frame(errors, row.names = c(teams[-which(teams == ref)], 'Alpha')),
                teams = teams, type = mod, ref = ref))
  }
  if(mod == 2){
    params = runif(length(teams) + max(rel))
    model = optim(par = params, fn = LL2, gr = LL2dash, mat = data, R=rel, method = 'BFGS')
    
    if(ref == 'last'){
      ref = teams[which.min(model$par[1:length(teams)])]
    }
    
    theta = model$par[1:length(teams)] - model$par[which(teams == ref)]
    alpha = model$par[(1 + length(teams)):(max(rel) + length(teams))]
    
    # Finding standard errors
    HH = Hess2(theta, alpha, data, rel)
    HH = HH[-which(teams == ref), -which(teams == ref)]
    se = Inverse(HH)
    
    # Taking only the diagonals
    d = col(se) - row(se)
    errors = sqrt(split(se, d)$'0')
    
    return(list(Theta = data.frame(Theta = theta, row.names = teams), 
                Alpha = data.frame(Alpha = alpha, row.names = 1:max(rel)),
                Std.Errors = data.frame(errors, row.names = c(teams[-which(teams == ref)], 1:max(rel))),
                teams = teams, type = mod, ref = ref))
  }
  if(mod == 3){
    params = runif(2*length(teams))
    model = optim(par = params, fn = LL4, gr = LL4dash, mat = data, method = 'BFGS')
    
    if(ref == 'last'){
      ref = teams[which.min(model$par[1:length(teams)])]
    } 
    
    theta = model$par[1:length(teams)] - model$par[which(teams == ref)]
    alpha = model$par[(1 + length(teams)):(2*length(teams))]
    
    # Std Errs
    HH = Hess4(theta, alpha, data)
    HH = HH[-which(teams == ref), -which(teams == ref)]
    se = Inverse(HH)
    d = col(se) - row(se)
    errors = sqrt(split(se, d)$'0')
    theta.errors = errors[1:(length(teams)-1)]
    alpha.errors = errors[length(teams):dim(HH)[1]]
    
    return(list(Theta = data.frame(Theta = theta, row.names = teams), 
                Alpha = data.frame(Alpha = alpha, row.names = teams),
                Theta.Errs = data.frame(theta.errors, row.names = teams[-which(teams == ref)]),
                Alpha.Errs = data.frame(alpha.errors, row.names = teams),
                teams = teams, type = mod, ref = ref))
  }
  if(mod == 4){
    params = runif((1+max(rel))*length(teams))
    model = optim(par = params, fn = LL5, gr = LL5dash, mat = data, R=rel, method = 'BFGS')
    
    if(ref == 'last'){
      ref = teams[which.min(model$par[1:length(teams)])]
    } 
    
    theta = model$par[1:length(teams)] - model$par[which(teams == ref)]
    Theta = data.frame(theta, row.names = teams)
    level = 1
    Alpha = NA
    while(level <= max(rel)){
      Level = model$par[(level*length(teams) + 1) : ((level+1)*length(teams))]
      Alpha = data.frame(Alpha, Level, row.names = teams)
      level = level + 1
    }
    return(list(Theta = Theta, 
                Alpha = Alpha[,-1],
                teams = teams, type = mod, ref = ref))
  }
  if(mod == 5){
    params = runif((length(teams)+1)*length(teams))
    model = optim(par = params, fn = LL6, gr = LL6dash, mat = data, method = 'BFGS')
    
    if(ref == 'last'){
      ref = teams[which.min(model$par[1:length(teams)])]
    } 
    
    model$par[1:length(teams)] = model$par[1:length(teams)] - model$par[which(teams == ref)]
    return(list(par = model$par, teams = teams, type = mod, ref = ref))
  }
  
}

BT_plot = function(model){
  ## Takes in a Bradley-Terry model from the BTC function and produces a plot
  ## of the parameters by team.
  
  ## For hierarchical models, it prints each plot every ten seconds, starting
  ## with H_1, up to H_n
  
  if(model$type == 0){
    # Vanilla Model
    vanilla = model$Theta
    params = vanilla[order(-vanilla$Theta),, drop = FALSE]
    print(ggplot(params, aes(x=reorder(rownames(params), Theta), y=Theta)) +
      geom_segment(aes(x=reorder(rownames(params), Theta), xend=reorder(rownames(params), Theta), y=0, yend=Theta)) +
      geom_point(color=rgb(0.2,0.7,0.1,0.5), size=3, alpha=0.6) +
      theme_light() +
      xlab('Teams') +
      coord_flip() +
      theme(
        panel.grid.major.y = element_blank(),
        panel.border = element_blank(),
        axis.ticks.y = element_blank()
      )
    )
  }
  if(model$type == 1){
    sing = data.frame(Theta = model$Theta[,], Home = model$Theta[,] + model$Alpha,
                      row.names = model$teams)
    print(
      ggplot(sing, aes(x=reorder(rownames(sing), Theta), y=Theta)) + 
        geom_segment(aes(x=reorder(rownames(sing), Theta), xend = reorder(rownames(sing), Theta), y = Theta, yend = Home)) +
        geom_point( aes(x=reorder(rownames(sing), Theta), y=Home), color=rgb(0.2,0.7,0.1,0.5), size=3 ) +
        geom_point( aes(x=reorder(rownames(sing), Theta), y=Theta), color=rgb(0.7,0.2,0.1,0.5), size=3 ) +
        coord_flip()+ 
        xlab('Teams') +
        theme_light() +
        theme(
          panel.grid.major.y = element_blank(),
          panel.border = element_blank(),
          axis.ticks.y = element_blank()
        )
    )
  }
  if(model$type == 2){
    ### TODO 
    ### - try and make it so user interrupts, rather than a 10 second sleep
    ### - insert titles that specify what level of the hierarchy has been plotted
    for(k in model$Alpha[,]){
      sing = data.frame(Theta = model$Theta[,], Home = model$Theta[,] + k,
                        row.names = model$teams)
      print(
        ggplot(sing, aes(x=reorder(rownames(sing), Theta), y=Theta)) + 
          geom_segment(aes(x=reorder(rownames(sing), Theta), xend = reorder(rownames(sing), Theta), y = Theta, yend = Home)) +
          geom_point( aes(x=reorder(rownames(sing), Theta), y=Home), color=rgb(0.2,0.7,0.1,0.5), size=3 ) +
          geom_point( aes(x=reorder(rownames(sing), Theta), y=Theta), color=rgb(0.7,0.2,0.1,0.5), size=3 ) +
          coord_flip()+ 
          xlab('Teams') +
          theme_light() +
          theme(
            panel.grid.major.y = element_blank(),
            panel.border = element_blank(),
            axis.ticks.y = element_blank()
          )
      )
      Sys.sleep(10)
    }
  }
  if(model$type == 3){
    # Non-Hierarchical Models
    ### TODO - fix this so that sing can work
    sing = data.frame(Theta = model$Theta[,], Home = model$Theta[,] + model$Alpha[,],
                      row.names = model$teams)
    print(
      ggplot(sing, aes(x=reorder(rownames(sing), Theta), y=Theta)) + 
      geom_segment(aes(x=reorder(rownames(sing), Theta), xend = reorder(rownames(sing), Theta), y = Theta, yend = Home)) +
      geom_point( aes(x=reorder(rownames(sing), Theta), y=Home), color=rgb(0.2,0.7,0.1,0.5), size=3 ) +
      geom_point( aes(x=reorder(rownames(sing), Theta), y=Theta), color=rgb(0.7,0.2,0.1,0.5), size=3 ) +
      coord_flip()+ 
      xlab('Teams') +
      theme_light() +
      theme(
        panel.grid.major.y = element_blank(),
        panel.border = element_blank(),
        axis.ticks.y = element_blank()
        )
    )
  }
  if(model$type == 4){
    ### TODO 
    ### - try and make it so user interrupts, rather than a 10 second sleep
    ### - insert titles that specify what level of the hierarchy has been plotted
    for(k in 1:length(model$Alpha)){
      sing = data.frame(Theta = model$Theta[,], 
                        Home = model$Theta[,] + model$Alpha[,k],
                        row.names = model$teams)
      print(
        ggplot(sing, aes(x=reorder(rownames(sing), Theta), y=Theta)) + 
          geom_segment(aes(x=reorder(rownames(sing), Theta), xend = reorder(rownames(sing), Theta), y = Theta, yend = Home)) +
          geom_point( aes(x=reorder(rownames(sing), Theta), y=Home), color=rgb(0.2,0.7,0.1,0.5), size=3 ) +
          geom_point( aes(x=reorder(rownames(sing), Theta), y=Theta), color=rgb(0.7,0.2,0.1,0.5), size=3 ) +
          coord_flip()+ 
          xlab('Teams') +
          theme_light() +
          theme(
            panel.grid.major.y = element_blank(),
            panel.border = element_blank(),
            axis.ticks.y = element_blank()
          )
      )
      Sys.sleep(10)
    }
  }
}

WLTable = function(mat){
  # Just a function that creates a Win/Loss table from the data given
  # Like the model, it handles draws by adding 0.5 to both columns
  WLtab = c()
  for(team in mat$teams){
    WLtab = rbind(WLtab, c(sum(mat$w[team,]) + sum(mat$l[, team]),
                           sum(mat$l[team,]) + sum(mat$w[, team]),
                           (sum(mat$w[team,]) + sum(mat$l[, team])) / (sum(mat$w[team,]) + sum(mat$l[, team]) + sum(mat$l[team,]) + sum(mat$w[, team]))))
  }
  WLtab = data.frame(WLtab, row.names = mat$teams)
  colnames(WLtab) = c('Wins','Losses','Win Rate')
  print(WLtab[order(-WLtab$`Win Rate`), ])
}

BT_summary = function(model){
  # Function for providing a summary of a Bradley-Terry model for our package
  rf = 3 # Rounding factor
  
  if(model$type == 0){
    # Vanilla Model - Just prints out all Thetas in descending order, with the
    # lowest team set to 0
    print(round(model$Theta[order(-model$Theta),, drop = FALSE], rf))
  }
  
  if(model$type == 1){
    print(round(model$Theta[order(-model$Theta),, drop = FALSE], rf))
    cat("With Home-Ground advantage equal to: ", round(model$Alpha, rf))
  }
  
  if(model$type == 2){
    print(round(model$Theta[order(-model$Theta),, drop = FALSE], rf))
    cat('With Hierarchical Home-Ground advantages of: ', '\n')
    print(round(model$Alpha, rf))
  }
  
  
  if(model$type == 3 | model$type == 4){
    # Unique Homeground Model - Prints out Theta and Alpha, ordered by Theta 
    print(cbind(round(model$Theta[order(-model$Theta), ,drop = FALSE], rf),
                round(model$Alpha[order(-model$Theta),,drop=FALSE], rf)))
  }
}

BT_predict = function(model, home_team, away_team, rel = NULL){
  # Predicts the probability of either side winning a match based off the fitted
  # BT model
  if(home_team %in% model$teams & away_team %in% model$teams){
    if(model$type == 0){
      home_win = exp(model$Theta[home_team,])/(exp(model$Theta[home_team,]) + exp(model$Theta[away_team,]))
    }
    if(model$type == 1){
      home_win = exp(model$Theta[home_team,] + model$Alpha)/(exp(model$Theta[home_team,]+ model$Alpha) + exp(model$Theta[away_team,]))
    }
    if(model$type == 2){
      home_win = exp(model$Theta[home_team,] + model$Alpha[rel[home_team, away_team],])/(exp(model$Theta[home_team,]+ model$Alpha[rel[home_team, away_team],]) + exp(model$Theta[away_team,]))
    }
    if(model$type == 3){
      home_win = exp(model$Theta[home_team,] + model$Alpha[home_team,])/(exp(model$Theta[home_team,]+ model$Alpha[home_team,]) + exp(model$Theta[away_team,]))
    }
    if(model$type == 4){
      home_win = exp(model$Theta[home_team,] + model$Alpha[home_team, rel[home_team, away_team]])/(exp(model$Theta[home_team,]+ model$Alpha[home_team, rel[home_team, away_team]]) + exp(model$Theta[away_team,]))
    }
    #cat(home_team, 'win with probability', round(home_win, 4), '\n')
    #cat(away_team, 'win with probability', round(1 - home_win, 4))
  return(home_win)
  }
  #else{
    #cat('These teams are not a part of the given model', '\n')
    #cat('Check for spelling errors, or make sure you have provided the correct model')
  #}
}

BT_test = function(model1, model2, mat, rel = NULL){
  ## Function for performing the hypothesis tests outlined in the thesis
  
  m1 = model1; m2 = model2
  if(model1$type > model2$type){
    m1 = model2; m2 = model1
  }
  if(model1$type == model2$type){
    return()
  }
  if(m1$type == 0 & m2$type == 1){
    stat = -2*(LL1(c(m2$Theta[,], m2$Alpha), mat) - LL(m1$Theta[,], mat))
    cat('P-Value is ', 1 - pchisq(stat, 1), '\n')
    cat('with a statistic of ', stat, 'and 1 degree of freedom')
  }
  if(m1$type == 0 & m2$type == 2){
    stat = -2*(LL2(H2Vec(m2), mat, R) - LL(m1$Theta[,], mat))
    df = length(H2Vec(m2)) - length(m1$teams)
    cat('P-Value is ', 1 - pchisq(stat, df), '\n')
    cat('with a statistic of ', stat, 'and', df, 'degrees of freedom')
  }
  if(m1$type == 0 & m2$type == 3){
    stat = -2*(LL4(c(m2$Theta[,], m2$Alpha[,]), mat) - LL(m1$Theta[,], mat))
    cat('P-Value is ', 1 - pchisq(stat, length(m1$teams)), '\n')
    cat('with a statistic of ', stat, 'and', length(m1$teams), 'degrees of freedom')
    }
  if(m1$type == 0 & m2$type == 4){
    stat = -2*(LL5(H2Vec(m2), mat, rel) - LL(m1$Theta[,], mat))
    deg = dim(m2$Alpha)[1] * dim(m2$Alpha)[2]
    cat('P-Value is ', (1 - pchisq(stat, deg)), '\n')
    cat('with a statistic of ', stat, 'and', deg, 'degrees of freedom')
  }
  if(m1$type == 0 & m2$type == 5){
    stat = -2*(LL6(m2$par, mat) - LL(m1$Theta[,], mat))
    df = (length(m1$teams))*(length(m1$teams) - 1)
    cat('P-Value is ', 1 - pchisq(stat, df), '\n')
    cat('with a statistic of ', stat, 'and', df, 'degrees of freedom.')
    
  }
  if(m1$type == 1 & m2$type == 2){
    stat = -2*(LL2(H2Vec(m2), mat, rel) - LL1(c(m1$Theta[,], m1$Alpha), mat))
    deg = length(H2Vec(m2)) - length(m1$teams) - 1
    cat('P-Value is ', (1 - pchisq(stat, deg)), '\n')
    cat('with a statistic of ', stat, 'and', deg, 'degrees of freedom')
  }
  if(m1$type == 1 & m2$type == 3){
    stat = -2*(LL4(c(m2$Theta[,], m2$Alpha[,]), mat) - LL1(c(m1$Theta[,], m1$Alpha), mat))
    cat('P-Value is ', (1 - pchisq(stat, length(m1$teams) - 1)), '\n')
    cat('with a statistic of ', stat, 'and', length(m1$teams) - 1, 'degrees of freedom')
  }
  if(m1$type == 1 & m2$type == 4){
    stat = -2*(LL5(H2Vec(m2), mat, rel) - LL1(c(m1$Theta[,], m1$Alpha), mat))
    deg = (dim(m2$Alpha)[1] * dim(m2$Alpha)[2]) - 1
    cat('P-Value is ', (1 - pchisq(stat, deg)), '\n')
    cat('with a statistic of ', stat, 'and', deg, 'degrees of freedom')
  }
  if(m1$type == 1 & m2$type == 5){
    stat = -2*(LL6(m2$par, mat) - LL1(c(m1$Theta[,], m1$Alpha), mat))
    deg = length(m2$par) - length(m1$teams) - 1
    cat('P-Value is ', (1 - pchisq(stat, deg)), '\n')
    cat('with a statistic of ', stat, 'and', deg, 'degrees of freedom')
  }
  if(m1$type == 2 & m2$type == 3){
    stat = -2*(LL4(H2Vec(m2), mat) - LL2(H2Vec(m1), mat, rel))
    deg = length(H2Vec(m2)) - length(H2Vec(m1))
    cat('P-Value is ', (1 - pchisq(stat, deg)), '\n')
    cat('with a statistic of ', stat, 'and', deg, 'degrees of freedom')
  }
  if(m1$type == 2 & m2$type == 4){
    stat = -2*(LL5(H2Vec(m2), mat, rel) - LL2(H2Vec(m1), mat, rel))
    deg = length(H2Vec(m2)) - length(H2Vec(m1))
    cat('P-Value is ', (1 - pchisq(stat, deg)), '\n')
    cat('with a statistic of ', stat, 'and', deg, 'degrees of freedom')
  }
  if(m1$type == 2 & m2$type == 5){
    stat = -2*(LL6(m2$par, mat) - LL2(H2Vec(m1), mat, rel))
    deg = length(pwise$par) - length(H2Vec(m1))
    cat('P-Value is ', (1 - pchisq(stat, deg)), '\n')
    cat('with a statistic of ', stat, 'and', deg, 'degrees of freedom')
  }
  if(m1$type == 3 & m2$type == 4){
    stat = -2*(LL5(H2Vec(m2), mat, rel) - LL4(c(m1$Theta[,], m1$Alpha[,]), mat))
    deg = dim(m2$Alpha)[1] * (dim(m2$Alpha)[2] - 1)
    cat('P-Value is ', (1 - pchisq(stat, deg)), '\n')
    cat('with a statistic of ', stat, 'and', deg, 'degrees of freedom')
    }
  if(m1$type == 3 & m2$type == 5){
    stat = -2*(LL6(m2$par, mat) - LL4(c(m1$Theta[,], m1$Alpha[,]), mat))
    cat('P-Value is ', 1 - pchisq(stat, (length(m1$teams) - 1)^2), '\n')
    cat('with a statistic of ', stat, 'and', (length(m1$teams) - 1)^2, 'degrees of freedom.')
  }
  if(m1$type == 4 & m2$type == 5){
    stat = -2*(LL6(m2$par, mat) - LL5(H2Vec(m1), mat, rel))
    df = length(m2$par) - length(H2Vec(m1))
    cat('P-Value is ', (1 - pchisq(stat, df)), '\n')
    cat('with a statistic of ', stat, 'and', df, 'degrees of freedom')
  }
}

BT_bootstrap = function(model, fixture, runs = 100, rel = NULL){
  # Run once, then run n-1 times
  print('Warning: runtime for BT_bootstrap() may be long for large datasets and/or number of runs')
  terrs = NULL; aerrs = NULL
  badteam = which.min(model$Theta[[1]])
  m = BTC(Data2Mat(SeasonSim(model, fixture, rel)), model$type, rel)
  theta = (m$Theta - m$Theta[[1]][badteam])
  if(model$type > 0){
    alpha = m$Alpha
  }
  for(run in 1:(runs - 1)){
    m = BTC(Data2Mat(SeasonSim(model, fixture, rel)), model$type, rel)
    theta = cbind(theta, m$Theta - m$Theta[[1]][badteam])
    if(model$type > 0){
      alpha = cbind(alpha, m$Alpha)
    }
  }
  terrs = apply(theta[,], 1, sd)
  if(model$type > 1){
    aerrs = apply(alpha[,], 1, sd)
  }
  if(model$type == 1){
    aerrs = sd(alpha)
  }
  return(c(terrs, aerrs))
}

Data2Mat = function(df){
  # Dataframe should have the columns 'home.team', 'away.team', 'home_scores' and 'away_scores'
  # This function then returns the h, w and l matrices required for our BT package
  
  teams = unique(df$home.team)
  
  df[,'home_scores'] = as.numeric(df[,'home_scores'])
  df[,'away_scores'] = as.numeric(df[,'away_scores'])
  
  
  w = matrix(data = 0, nrow = length(teams), ncol = length(teams), dimnames = list(teams, teams))
  l = matrix(data = 0, nrow = length(teams), ncol = length(teams), dimnames = list(teams, teams))
  
  for(i in 1:nrow(df)){
    row = df[i,]
    w[row$home.team, row$away.team] = w[row$home.team, row$away.team] + (row$home_scores > row$away_scores) + 0.5*(row$home_scores == row$away_scores)
    l[row$home.team, row$away.team] = l[row$home.team, row$away.team] + (row$home_scores < row$away_scores) + 0.5*(row$home_scores == row$away_scores)
  }
  
  h = w + l
  return(list(h=h, w=w, l=l, teams = teams))
}
ResultSplitter = function(data){
  ## Takes data with a 'Results' column only, of the form '1 - 0' when the Home
  ## team wins 1 - 0, and converts it into a compatable dataframe.
  ##
  ## Data from fixturedownload.com, such as mlb-2020-UTC.csv, come in this form
  ## (see MLB2020.R)
  
  home.score = c()
  away.score = c()
  for(k in 1:nrow(data)){
    home.score = c(home.score, strsplit(data[k, ]$Result, ' - ')[[1]][1])
    away.score = c(away.score, strsplit(data[k, ]$Result, ' - ')[[1]][2])
  }
  
  data = data.frame(home.team = data$Home.Team, away.team = data$Away.Team, 
                    home_scores = home.score, away_scores = away.score)
  
  return(data)
}
