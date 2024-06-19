
Stab_Hier_new <- function (data, mat, order.q = c(1, 2), Alltime = TRUE, start_T = NULL, end_T = NULL){
  if (Alltime == FALSE) {
    data <- data[, start_T:end_T]
  }
  data <- as.matrix(data)
  data[which(data==0)] <- 10^(-15)
  data <- as.data.frame(data)
  
  TT <- ncol(data)
  gamma_zz <- apply(data, 2, sum)
  gamma <- sapply(order.q, function(qq) {
    if (qq == 1) {
      gamma_zz <- gamma_zz[gamma_zz != 0]
      gamma <- (-1) * sum((gamma_zz/sum(gamma_zz)) * log((gamma_zz/sum(gamma_zz))))/log(TT)
    }
    else {
      gamma <- (1 - sum((gamma_zz/sum(gamma_zz))^qq))/(1 - TT^(1 - qq))
    }
  })
  gamma <- data.frame(hier = ncol(mat) + 1, order_q = order.q, gamma_value = gamma)
  
  hier_i_alphaanddiff <- function(data, mat, q, i) {
    TT <- ncol(data)
    index <- unique(mat[, 1:i])
    if (i == 1) {
      take <- sapply(1:length(index), function(rr) which(mat[,1:i] == index[rr]))
    }else {
      take <- sapply(1:nrow(index), function(rr) which(apply(mat[,1:i], 1, function(mm) sum(mm == index[rr, ])) == i))
    }
    value <- lapply(take, function(tak) {
      zz <- data[tak, ]
      zz <- apply(zz, 2, sum)
      outout <- sapply(q, function(qq) {
        if (length(zz[zz != 0]) != 0) {
          if (qq == 1) {
            zz <- zz[zz != 0]
            out <- (-1) * sum((zz/sum(zz)) * log((zz/sum(zz))))/log(TT)
          }
          else {
            out <- (1 - sum((zz/sum(zz))^qq))/(1 - TT^(1 - qq))
          }
          ww <- sum(zz)/sum(data)
          vv <- ww * out
        }
        else {
          out <- 0
          ww <- 0
          vv <- ww * out
        }
        return(vv)
      })
      return(outout)
    })
    value <- do.call(rbind, value)
    if (length(q) == 1) {
      result <- sum(value)
    }else {
      result <- apply(value, 2, sum)
    }
    
    if (i >= 2) {
      index2 <- unique(mat[, 1:(i - 1)])
      if (i == 2) {
        take2 <- sapply(1:length(index2), function(rr2) {
          which(index[, 1:(i - 1)] == index2[rr2])
        })
      }else{
        take2 <- sapply(1:nrow(index2), function(rr2) {
          which(apply(index[, 1:(i - 1)], 1, function(mm) sum(mm == index2[rr2, ])) == (i - 1))
        })
      }
      
      subdat <- lapply(take, function(tak) {
        zz <- data[tak, ]
        zz <- apply(zz, 2, sum)
        return(zz)
      })
      
      if(i == 2 & length(index2)<=2){
        subdat2 <- apply(take2, 2, function(tak2) {
          zz <- subdat[tak2]
          zz <- do.call(rbind, zz)
          
          vv <- sapply(q, function(qq) {
            if (qq == 1) {
              zz[which(zz == 0)] <- 10^(-15)
              z_ik <- apply(zz,1,sum)
              w_ik <- z_ik/sum(data)
              w_plusk <- sum(zz)/sum(data)
              value <- sum((-1)*w_ik*log(w_ik/w_plusk)/log(TT))
              
            }else {
              back <- t(apply(zz, 1, function(z)(z/sum(z))^qq))
              z_ik <- apply(zz,1,sum)
              mid <- (z_ik/sum(z_ik))-(z_ik/sum(z_ik))^qq
              front <- sum(zz)/sum(data)
              value <- sum(apply(back,2,function(z)z*mid*front)/(1 - TT^(1 - qq)))
              
            }
            return(value)
          })
          return(vv)
        })
        subdat2 <- t(subdat2)
        
      }else{
        subdat2 <- lapply(take2, function(tak2) {
          zz <- subdat[tak2]
          zz <- do.call(rbind, zz)
          vv <- sapply(q, function(qq) {
            if (qq == 1) {
              zz[which(zz == 0)] <- 10^(-15)
              z_ik <- apply(zz,1,sum)
              w_ik <- z_ik/sum(data)
              w_plusk <- sum(zz)/sum(data)
              value <- sum((-1)*w_ik*log(w_ik/w_plusk)/log(TT))
              
            }else {
              back <- t(apply(zz, 1, function(z)(z/sum(z))^qq))
              z_ik <- apply(zz,1,sum)
              mid <- (z_ik/sum(z_ik))-(z_ik/sum(z_ik))^qq
              front <- sum(zz)/sum(data)
              value <- sum(apply(back,2,function(z)z*mid*front)/(1 - TT^(1 - qq)))
              
            }
            return(value)
          })
          return(vv)
        })
        subdat2 <- do.call(rbind, subdat2)
      }
      
      if (length(q) == 1) {
        subdat2 <- sum(subdat2)
      }
      else {
        subdat2 <- apply(subdat2, 2, sum)
      }
    }else {
      subdat2 <- lapply(take, function(tak) {
        zz <- data[tak, ]
        zz <- apply(zz, 2, sum)
        vv <- sapply(q, function(qq) {
          if (length(zz[zz != 0]) != 0) {
            if (qq == 1) {
              w_plusk <- (sum(zz)/sum(data))
              out <- (-1)*w_plusk*log(w_plusk)/log(TT)
            }
            else {
              back <- (zz/sum(zz))^qq
              mid <- (sum(zz)/sum(data))-(sum(zz)/sum(data))^qq
              out <- sum(mid*back/(1-TT^(1-qq)))
            }
          }
          else {
            out <- 0
          }
          return(out)
        })
        return(vv)
      })
      subdat2 <- do.call(rbind, subdat2)
      subdat2 <- apply(subdat2,2,sum)
    }
    return(data.frame(hier = rep((ncol(mat) + 1 - i), length(q)), 
                      order_q = q, gamma_alpha = result, beta_max = subdat2))
  }
  stab <- c()
  for (ii in 1:ncol(mat)) {
    stab <- rbind(stab, hier_i_alphaanddiff(data, mat, order.q,ii))
  }
  gamma_value <- rbind(gamma, gamma, data.frame(stab[, 1:2],gamma_value = stab[, 3]))
  gamma_value <- filter(gamma_value, hier != 1)
  gamma_value$hier <- rep(c((ncol(mat) + 1):1), each = length(order.q))
  alpha_value <- c(rep(NA, length(order.q)), stab[, 3])
  beta_max <- c(rep(NA, length(order.q)), stab[, 4])
  alldata <- cbind(gamma_value, alpha_value, beta_max)
  alldata <- cbind(alldata[, 1:4], beta_value = (alldata[,3] - alldata[, 4])/alldata[, 5])
  alldata <- as.data.frame(alldata)
  colnames(alldata) <- c("Hier", "Order_q", "Gamma", "Alpha", 
                         "Beta (normalized)")
  return(alldata)
}



OUT2 <- Stab_Hier_new(Jena_hierarchical_data, Jena_hierarchical_mat, order.q = seq(0.1,2,0.1), Alltime = TRUE)
OUT1 <- Stab_Hier(Jena_hierarchical_data, Jena_hierarchical_mat, order.q = seq(0.1,2,0.1), Alltime = TRUE)
