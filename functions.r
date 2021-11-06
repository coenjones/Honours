SimpleHierarchy = function(groups, mat){
  # Groups should be a matrix where each ROW represents a group of teams
  # This function will then return a relationship matrix, with a 2 allocated
  # to all pairs where the two items are in the same group, a 1 when in 
  # seperate groups, and a 0 when both items are the same.
  # Any pairs where one or more teams are not in ANY group given will also 
  # be marked with a 1.
  
  n.teams = length(mat$teams)
  
  R = matrix(rep(1, n.teams^2), nrow = n.teams, dimnames = c(list(mat$teams), list(mat$teams)))
  for(grp in 1:nrow(groups)){
    for(t1 in 1:length(groups[grp,])){
      for(t2 in 1:length(groups[grp,])){
        if(groups[grp, t1] != groups[grp, t2])
          R[groups[grp, t1], groups[grp, t2]] = 2
      }
    }
  }
  for(t1 in 1:n.teams){
    R[t1, t1] = 0
  }
  return(R)
}
SeasonSim = function(model, fixture, rel = NULL){
  # Takes a set of fixtures and simulates results for each game using a previously fitted BT Model
  newfix = data.frame(fixture)
  for(match in 1:nrow(newfix)){
    newfix[match, 'home_scores'] = (BT_predict(model, 
                                               newfix[match, 'home.team'],
                                               newfix[match, 'away.team'], rel) > runif(1))
    newfix[match, 'away_scores'] = 1 - newfix[match, 'home_scores']
  }
  return(newfix)
}
MatchSim = function(model, home_team, away_team, rel = NULL){
  p = BT_predict(model, home_team, away_team, rel = rel)
  if(p > runif(1)){
    return(home_team)
  }
  else{
    return(away_team)
  }
}


# Model 0 - No Homeground
LL = function(params, mat){
  # mat contains mat$h, mat$w, mat$l
  # params is the vector of thetas
  teams = colnames(mat$h)
  
  log_like = 0 
  
  for(t1 in 1:length(teams)){
    for(t2 in 1:length(teams)){
      log_like = log_like + mat$w[teams[t1], teams[t2]]*params[t1] + mat$l[teams[t1],teams[t2]]*params[t2] - mat$h[teams[t1],teams[t2]]*log(exp(params[t1]) + exp(params[t2]))
    }
  }
  return(-1*log_like)
}
LLdash = function(params, mat){
  teams = colnames(mat$h)
  output = rep(0, length(teams))
  for(t1 in 1:length(teams)){
    lastterm = 0
    for(t2 in 1:length(teams)){
      lastterm = lastterm - (mat$h[t1,t2] + mat$h[t2,t1])*(exp(params[t1])/(exp(params[t1]) + exp(params[t2])))
    }
    output[t1] = sum(mat$w[t1,]) + sum(mat$l[,t1]) + lastterm
  }
  return(-1*output)
}
Hess = function(params, mat){
  t = length(mat$teams)
  output = matrix(rep(0, t*t), t, t)
  for(i in 1:t){
    for(j in 1:t){
      if(i == j){
        for(k in 1:t){
          output[i, j] = output[i, j] + ((mat$h[i, k] + mat$h[k, i]) * 
                                           ((exp(params[i] + params[k])) / 
                                           ((exp(params[i]) + exp(params[k]))^2)))
        }
      }
      else{
        output[i, j] = -1*((mat$h[i, j] + mat$h[j, i]) * 
                        (exp(params[i] + params[j])) / 
                        ((exp(params[i]) + exp(params[j]))^2))
      }
    }
  }
  return(output)
}


# Model 1 - Common Homeground
LL1 = function(params, mat){
  # mat contains mat$h, mat$w, mat$l
  # params is the vector of thetas, followed by the vector of deltas.
  teams = colnames(mat$h)
  TM = length(teams)
  
  log_like = 0 
  
  for(t1 in 1:length(teams)){
    for(t2 in 1:length(teams)){
      log_like = log_like + mat$w[teams[t1], teams[t2]]*(params[t1] + params[1 + TM]) + 
        mat$l[teams[t1],teams[t2]]*params[t2] - 
        mat$h[teams[t1],teams[t2]]*log(exp(params[t1] + params[1 + TM]) + exp(params[t2]))
    }
  }
  return(-1*log_like)
}
LL1dash = function(params, mat){
  teams = colnames(mat$h)
  TM = length(teams)
  output = params
  
  # dL/dTheta
  for(t1 in 1:TM){
    lastterm = 0
    for(t2 in 1:TM){
      lastterm = lastterm - 
        (mat$h[t1,t2])*(exp(params[t1] + params[1+TM])/(exp(params[t1] + params[1+TM]) + exp(params[t2]))) -
        (mat$h[t2,t1])*(exp(params[t1])/(exp(params[t2] + params[1+TM]) + exp(params[t1])))
    }
    output[t1] = sum(mat$w[t1,]) + sum(mat$l[,t1]) + lastterm
  }
  
  # dL/dAlpha
  output[t1 + 1] = 0
  for(t1 in 1:TM){
    lastterm = 0
    for(t2 in 1:TM){
      lastterm = lastterm - (mat$h[t1,t2])*(exp(params[t1] + params[1+TM])/(exp(params[t1] + params[1+TM]) + exp(params[t2])))
    }
    output[TM + 1] = output[TM + 1] + sum(mat$w[t1,]) + lastterm
  }
  return(-1*output)
}
Hess1 = function(theta, alpha, mat){
  t = length(mat$teams)
  output = matrix(rep(0, (t+1)^2), t+1, t+1)
  for(i in 1:t){
    for(j in 1:t){
      if(i == j){
        for(k in 1:t){
          output[i, j] = output[i, j] + exp(theta[i] + theta[k] + alpha) * 
                            (mat$h[i,k] / (exp(theta[i] + alpha) + exp(theta[k]))^2 +
                             mat$h[k,i] / (exp(theta[k] + alpha) + exp(theta[i]))^2)
        }
      }
      else{
        output[i, j] = -1*exp(theta[i] + theta[j] + alpha) * 
                          ((mat$h[i,j] / (exp(theta[i] + alpha) + exp(theta[j]))^2) +
                          (mat$h[j,i] / (exp(theta[j] + alpha) + exp(theta[i]))^2))
      }
    }
  }
  for(i in 1:t){
    for(j in 1:t){
      output[i, t+1] = output[i, t+1] + (exp(theta[i] + theta[j] + alpha) * 
                                           ((mat$h[i,j] / (exp(theta[i] + alpha) + exp(theta[j]))^2) -
                                            (mat$h[j,i] / (exp(theta[j] + alpha) + exp(theta[i]))^2)))
      output[t+1, t+1] = output[t+1, t+1] + (exp(theta[i] + theta[j] + alpha) * mat$h[i,j] / 
                                               (exp(theta[i] + alpha) + exp(theta[j]))^2)
    }
    output[t+1, i] = output[i, t+1]
  }
  return(output)
}


# Model 2 - Common Hierarchical
LL2 = function(params, mat, R){
  # mat contains mat$h, mat$w, mat$l
  # params is the vector of thetas, followed by the vector of deltas.
  teams = colnames(mat$h)
  TM = length(teams)
  
  log_like = 0 
  
  for(t1 in 1:length(teams)){
    for(t2 in 1:length(teams)){
      log_like = log_like + mat$w[teams[t1], teams[t2]]*(params[t1] + params[R[t1, t2] + TM]) + 
        mat$l[teams[t1],teams[t2]]*params[t2] - 
        mat$h[teams[t1],teams[t2]]*log(exp(params[t1] + params[R[t1, t2] + TM]) + exp(params[t2]))
    }
  }
  return(-1*log_like)
}
LL2dash = function(params, mat, R){
  teams = colnames(mat$h)
  TM = length(teams)
  output = rep(0, TM+max(R))
  
  # dL/dTheta
  for(t1 in 1:TM){
    lastterm = 0
    for(t2 in 1:TM){
      lastterm = lastterm - 
        (mat$h[t1,t2])*(exp(params[t1] + params[R[t1, t2]+TM])/(exp(params[t1] + params[R[t1, t2]+TM]) + exp(params[t2]))) -
        (mat$h[t2,t1])*(exp(params[t1])/(exp(params[t2] + params[R[t1, t2]+TM]) + exp(params[t1])))
    }
    output[t1] = sum(mat$w[t1,]) + sum(mat$l[,t1]) + lastterm
  }
  
  # dL/dAlpha
  for(k in 1:max(R)){
    output[TM + k] = 0
  }
  
  for(t1 in 1:TM){
    for(t2 in 1:TM){
      output[TM + R[t1, t2]] = output[TM + R[t1, t2]] + mat$w[t1, t2] -
        (mat$h[t1,t2])*(exp(params[t1] + params[R[t1, t2] +TM])/(exp(params[t1] + params[R[t1, t2]+TM]) + exp(params[t2])))
    }
  }
  return(-1*output)
}
Hess2 = function(theta, alpha, mat, R){
  t = length(mat$teams)
  n = length(alpha)
  output = matrix(rep(0, (t + n)^2), t + n, t + n)
  for(i in 1:t){
    for(j in 1:t){
      if(i == j){
        for(k in 1:t){
          output[i, i] = output[i, i] + exp(theta[i] + theta[k] + alpha[max(R[i,k], 1)]) * 
            (mat$h[i,k] / (exp(theta[i] + alpha[max(R[i,k], 1)]) + exp(theta[k]))^2 +
             mat$h[k,i] / (exp(theta[k] + alpha[max(R[i,k], 1)]) + exp(theta[i]))^2)
        }
      }
      else{
        output[i, j] = -1*exp(theta[i] + theta[j] + alpha[max(R[i,k], 1)]) * 
          ((mat$h[i,j] / (exp(theta[i] + alpha[max(R[i,k], 1)]) + exp(theta[j]))^2) +
             (mat$h[j,i]/ (exp(theta[j] + alpha[max(R[i,k], 1)]) + exp(theta[i]))^2))
      }
      for(lev in 1:n){
        output[lev + t, lev + t] = output[lev + t, lev + t] + 
          ((mat$h[i,j]*exp(theta[i] + theta[j] + alpha[lev]) / 
              (exp(theta[i] + alpha[lev]) + exp(theta[j]))^2) * 
             (lev == R[i,j]))
        output[i, lev + t] = output[i, lev + t] + 
          ((exp(theta[i] + theta[j] + alpha[lev]) * 
              (mat$h[i,j] / (exp(theta[i] + alpha[lev]) + exp(theta[j]))^2 -
                 mat$h[j,i] / (exp(theta[j] + alpha[lev]) + exp(theta[i]))^2)) *
             (lev == R[i,j]))
      }
    }
    for(lev in 1:n){
      output[lev + t, i] = output[i, lev + t]
    }
  }
  return(output)
}


# Model 3 - Unique Homeground
LL4 = function(params, mat){
  # mat contains mat$h, mat$w, mat$l
  # params is the vector of thetas, followed by the vector of deltas.
  teams = colnames(mat$h)
  TM = length(teams)
  
  log_like = 0 
  
  for(t1 in 1:length(teams)){
    for(t2 in 1:length(teams)){
      log_like = log_like + mat$w[teams[t1], teams[t2]]*(params[t1] + params[t1 + TM]) + 
        mat$l[teams[t1],teams[t2]]*params[t2] - 
        mat$h[teams[t1],teams[t2]]*log(exp(params[t1] + params[t1 + TM]) + exp(params[t2]))
    }
  }
  return(-1*log_like)
}
LL4dash = function(params, mat){
  teams = colnames(mat$h)
  TM = length(teams)
  output = rep(0, 2*TM)
  
  # dL/dTheta
  for(t1 in 1:TM){
    lastterm = 0
    for(t2 in 1:TM){
      lastterm = lastterm - 
        (mat$h[t1,t2])*(exp(params[t1] + params[t1+TM])/(exp(params[t1] + params[t1+TM]) + exp(params[t2]))) -
        (mat$h[t2,t1])*(exp(params[t1])/(exp(params[t2] + params[t2+TM]) + exp(params[t1])))
    }
    output[t1] = sum(mat$w[t1,]) + sum(mat$l[,t1]) + lastterm
  }
  
  # dL/dAlpha
  for(t1 in 1:TM){
    lastterm = 0
    for(t2 in 1:TM){
      lastterm = lastterm - (mat$h[t1,t2])*(exp(params[t1] + params[t1+TM])/(exp(params[t1] + params[t1+TM]) + exp(params[t2])))
    }
    output[t1 + TM] = sum(mat$w[t1,]) + lastterm
  }
  return(-1*output)
}
Hess4 = function(theta, alpha, mat){
  t = length(mat$teams)
  output = matrix(rep(0, (2*t)^2), 2*t, 2*t)
  for(i in 1:t){
    for(j in 1:t){
      if(i == j){
        for(k in 1:t){
          output[i, i] = output[i, i] + exp(theta[i] + theta[k]) * 
                        (mat$h[i,k]*exp(alpha[i]) / (exp(theta[i] + alpha[i]) + exp(theta[k]))^2 +
                         mat$h[k,i]*exp(alpha[k]) / (exp(theta[k] + alpha[k]) + exp(theta[i]))^2)
          output[t+i, t+i] = output[t+i, t+i] + exp(theta[i] + theta[k] + alpha[i]) * 
                              (mat$h[i,k] / (exp(theta[i] + alpha[i]) + exp(theta[k]))^2)
        }
      }
      else{
        output[i, j] = -1*exp(theta[i] + theta[j]) * 
                       ((mat$h[i,j]*exp(alpha[i]) / (exp(theta[i] + alpha[i]) + exp(theta[j]))^2) +
                        (mat$h[j,i]*exp(alpha[j]) / (exp(theta[j] + alpha[j]) + exp(theta[i]))^2))
        output[i, j+t] = mat$h[j,i]*exp(theta[i] + theta[j] + alpha[j])/((exp(theta[j]+alpha[j]) + exp(theta[i]))^2) 
        output[j+t, i] = output[i, j+t]
      }
    }
    output[i, t+i] = -1*output[t+i, t+i]
    output[t+i, i] = output[i, t+i]
  }
  return(output)
}


# Model 4 - Hierarchical 
LL5 = function(params, mat, R){
  # mat contains mat$h, mat$w, mat$l
  # params is the vector of thetas, followed by the vector of deltas.
  teams = colnames(mat$h)
  TM = length(teams)
  
  log_like = 0 
  
  for(t1 in 1:length(teams)){
    for(t2 in (1:length(teams))[-t1]){
      log_like = log_like + mat$w[teams[t1], teams[t2]]*(params[t1] + params[R[t1, t2]*TM + t1]) + 
        mat$l[teams[t1],teams[t2]]*params[t2] - 
        mat$h[teams[t1],teams[t2]]*log(exp(params[t1] + params[R[t1, t2]*TM + t1]) + exp(params[t2]))
    }
  }
  return(-1*log_like)
}
LL5dash = function(params, mat, R){
  teams = colnames(mat$h)
  TM = length(teams)
  output = rep(0, length(params))
  index = 1:TM
  
  
  # dL/dTheta
  for(t1 in 1:TM){
    dldt = 0
    for(t2 in (1:TM)[-t1]){
      dldt = dldt + mat$w[teams[t1],teams[t2]] + mat$l[teams[t2],teams[t1]] - 
        (mat$h[teams[t1],teams[t2]] * exp(params[t1] + params[t1 + TM*R[t1,t2]])/(exp(params[t2]) + exp(params[t1] + params[t1 + TM*R[t1,t2]]))) -
        (mat$h[teams[t2],teams[t1]] * exp(params[t1])/(exp(params[t2] + params[t2 + TM*R[t1,t2]]) + exp(params[t1])))
        }
    output[t1] = dldt
  }
  
  # dL/dAlpha
  for(t1 in 1:TM){
    for(t2 in 1:TM){
      output[t1 + TM*R[t1, t2]] = output[t1 + TM*R[t1, t2]] + mat$w[t1, t2] -
        (mat$h[t1,t2])*(exp(params[t1] + params[t1 + TM*R[t1, t2]])/(exp(params[t1] + params[t1 + TM*R[t1, t2]]) + exp(params[t2])))
    }
  }
  return(-1*output)
}
Hess5 = function(theta, alpha, mat, R){
  t = length(mat$teams)
  n = max(R)
  output = matrix(rep(0, ((n+1)*t)^2), (n + 1)*t, (n + 1)*t)
  for(i in 1:t){
    for(j in 1:t){
      if(i == j){
        for(k in 1:t){
          # dThetai^2
          output[i, i] = output[i, i] + exp(theta[i] + theta[k]) * 
            ((mat$h[i,k] * exp(alpha[i + t*(max(R[i,k], 1) - 1)]) /
             (exp(theta[i] + alpha[i + t*(max(R[i,k], 1) - 1)]) + exp(theta[k]))^2) + 
             (mat$h[k,i] * exp(alpha[k + t*(max(R[i,k], 1) - 1)]) /
             (exp(theta[k] + alpha[k + t*(max(R[i,k], 1) - 1)]) + exp(theta[i]))^2))
          # dThetai dAlphai
          output[i, i + t*max(R[i,k], 1)] = output[i, i + t*max(R[i,k], 1)] + 
            mat$h[i, k] * exp(theta[i] + theta[k] + alpha[i + t*(max(R[i,k], 1) - 1)]) / 
            ((exp(theta[i] + alpha[i + t*(max(R[i,k], 1) - 1)]) + exp(theta[k]))^2)
          output[i + t*max(R[i,k], 1), i] = output[i + t*max(R[i,k], 1), i] + 
            mat$h[i, k] * exp(theta[i] + theta[k] + alpha[i + t*(max(R[i,k], 1) - 1)]) / 
            ((exp(theta[i] + alpha[i + t*(max(R[i,k], 1) - 1)]) + exp(theta[k]))^2)
        }
      }
      else{
        # dThetai dThetaj
        output[i, j] = -1*exp(theta[i] + theta[j]) * 
          ((mat$h[i,j] * exp(alpha[i + t*(max(R[i,j], 1) - 1)]) /
           (exp(theta[i] + alpha[i + t*(max(R[i,j], 1) - 1)]) + exp(theta[j]))^2) + 
           (mat$h[j,i] * exp(alpha[j + t*(max(R[i,j], 1) - 1)]) /
           (exp(theta[j] + alpha[j + t*(max(R[i,j], 1) - 1)]) + exp(theta[i]))^2))
        # dThetai dAlphaj
        output[i, j + t*max(R[i,j], 1)] = -1 * mat$h[j,i] * 
          exp(theta[i] + theta[j] + alpha[j + t*(max(R[i,j] ,1) - 1)]) / 
          ((exp(theta[j] + alpha[j + t*(max(R[i,j] ,1) - 1)]) + exp(theta[i]))^2)
        output[j + t*(max(R[i,j], 1)), i] = output[i, j + t*(max(R[i,j] ,1))]
        # dAlphai^2
        output[i + t*R[i,j], i + t*R[i,j]] = output[i + t*R[i,j], i + t*R[i,j]] + 
          mat$h[i, j] * exp(theta[i] + theta[j] + alpha[i + t*(R[i,j] - 1)]) / 
          ((exp(theta[i] + alpha[i + t*(R[i,j] - 1)]) + exp(theta[j]))^2)
      }
    }
  }
  return(output)
}

# Model 5 - Pairwise Homegrounds 
LL6 = function(params, mat){
  # mat contains mat$h, mat$w, mat$l
  # params is the vector of thetas, followed by the vector of deltas.
  # this is the unique pairwise home advantage model
  teams = colnames(mat$h)
  TM = length(teams)
  
  log_like = 0 
  
  for(t1 in 1:length(teams)){
    for(t2 in 1:length(teams)){
      log_like = log_like + mat$w[teams[t1], teams[t2]]*(params[t1] + params[t1*TM + t2]) + 
        mat$l[teams[t1],teams[t2]]*params[t2] - 
        mat$h[teams[t1],teams[t2]]*log(exp(params[t1] + params[t1*TM + t2]) + exp(params[t2]))
    }
  }
  return(-1*log_like)
}
LL6dash = function(params, mat){
  teams = colnames(mat$h)
  TM = length(teams)
  output = rep(0, (TM+1)*TM)
  
  # dL/dTheta
  for(t1 in 1:TM){
    lastterm = 0
    for(t2 in 1:TM){
      lastterm = lastterm - 
        (mat$h[t1,t2])*(exp(params[t1] + params[t1*TM + t2])/(exp(params[t1] + params[t1*TM + t2]) + exp(params[t2]))) -
        (mat$h[t2,t1])*(exp(params[t1])/(exp(params[t2] + params[t2*TM + t1]) + exp(params[t1])))
    }
    output[t1] = sum(mat$w[t1,]) + sum(mat$l[,t1]) + lastterm
  }
  
  # dL/dAlpha
  for(t1 in 1:TM){
    for(t2 in 1:TM){
      output[t1*TM + t2] = mat$w[t1, t2] - mat$h[t1,t2]*(exp(params[t1] + params[t1*TM + t2])/(exp(params[t1] + params[t1*TM + t2]) + exp(params[t2])))
    }
  }
  return(-1*output)
}



NBA_relat = function(){
  # Function for creating the NBA relationship matrix
  data = read.csv('Data/NBA16-20.csv')
  teams = unique(data$home.team)
  ATL = c('Boston Celtics',
          'Brooklyn Nets',
          'New York Knicks',
          'Philadelphia 76ers',
          'Toronto Raptors')
  CNT = c('Chicago Bulls',
          'Cleveland Cavaliers',
          'Detroit Pistons',
          'Indiana Pacers',
          'Milwaukee Bucks')
  SE = c('Atlanta Hawks',
         'Charlotte Hornets',
         'Miami Heat',
         'Orlando Magic',
         'Washington Wizards')
  NW = c('Denver Nuggets',
         'Minnesota Timberwolves',
         'Oklahoma City Thunder',
         'Portland Trail Blazers',
         'Utah Jazz')
  PAC = c('Golden State Warriors',
          'Los Angeles Clippers',
          'Los Angeles Lakers',
          'Phoenix Suns',
          'Sacramento Kings')
  SW = c('Dallas Mavericks',
         'Houston Rockets',
         'Memphis Grizzlies',
         'New Orleans Pelicans',
         'San Antonio Spurs')
  
  D = matrix(data = 0, nrow = length(teams), ncol = length(teams), dimnames = list(teams, teams))
  C = matrix(data = 0, nrow = length(teams), ncol = length(teams), dimnames = list(teams, teams))
  L = matrix(data = 0, nrow = length(teams), ncol = length(teams), dimnames = list(teams, teams))
  
  east = data.frame(ATL, SE, CNT)
  west = data.frame(PAC, SW, NW)
  nba = data.frame(east, west)
  
  # Designating Divisions
  for(div in 1:ncol(nba)){
    for(t1 in 1:nrow(nba[div])){
      for(t2 in 1:nrow(nba[div])){
        if(nba[t1, div] != nba[t2, div])
          D[nba[t1,div], nba[t2, div]] = 1
      }
    }
  }
  
  # Counting Conferences
  for(div1 in 1:ncol(east)){
    for(div2 in 1:ncol(east)){
      for(t1 in 1:nrow(east[div1])){
        for(t2 in 1:nrow(east[div2])){
          if(nba[t1, div1] != nba[t2, div2]){
            C[nba[t1,div1], nba[t2, div2]] = 1
          }
        }
      }
    }
  }
  for(div1 in 1:ncol(west)){
    for(div2 in 1:ncol(west)){
      for(t1 in 1:nrow(west[div1])){
        for(t2 in 1:nrow(west[div2])){
          if(nba[t1, div1 + 3] != nba[t2, div2 + 3]){
            C[nba[t1,div1 + 3], nba[t2, div2 + 3]] = 1
          }  
        }
      }
    }
  }
  
  for(t1 in teams){
    for(t2 in teams){
      if(t1 != t2){
        L[t1, t2] = 1
      }
    }
  }
  
  R = C + D + L
  
  return(R)
}

H2Vec = function(model){
  # Takes the output of a hierarchical model and converts the parameters into
  # a vector that can be used for future inputs
  num.h = dim(model$Alpha)[2]
  output = c(model$Theta[,])
  for(k in 1:num.h){
    output = c(output, model$Alpha[,k])
  }
  return(output)
}

nparams = function(mat, model_type, rel = NULL){
  if(model_type == 0){return(length(mat$teams))}
  if(model_type == 1){return(length(mat$teams) + 1)}
  if(model_type == 2){return(length(mat$teams) + max(rel))}
  if(model_type == 3){return(2*length(mat$teams))}
  if(model_type == 4){return((1 + max(rel))*length(mat$teams))}
  if(model_type == 5){return(length(mat$teams)*(length(mat$teams) + 1))}
}
