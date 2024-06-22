#' Calculate stability of the time series data for single assemblage.
#'
#' \code{Stab_Single} is a function that calculate stability of the time series data (like biomass, productivity, etc.) for single assemblage.
#'
#' @param data can be input as a \code{vector} of time series data, or \code{data.frame} (assemblages by times).
#' @param order.q a numerical vector specifying the orders of stability. Default is c(1,2).
#' @param Alltime \code{TRUE} or \code{FALSE}, to decide whether to use all the times in the data.
#' @param start_T (argument only for \code{Alltime = FALSE}) a positive integer specifying the start column of time in the data.
#' @param end_T (argument only for \code{Alltime = FALSE}) a positive integer specifying the end column of time in the data.
#'
#'
#'
#' @return a dataframe with columns "Plot/Community", "Order_q" and "Stability".
#'
#' @examples
#' # Stability of each single plot
#' data("Jena_plot_biomass")
#' single_plot <- do.call(rbind, Jena_plot_biomass)
#' output_single_plot <- Stab_Single(data=single_plot, order.q=c(1,2), Alltime=TRUE)
#' output_single_plot
#'
#' # Stability of each single species in each plot
#' data("Jena_species_biomass")
#' single_species <- do.call(rbind, Jena_species_biomass)
#' output_single_species <- Stab_Single(data=single_species, order.q=c(1,2), Alltime=TRUE)
#' output_single_species
#'
#' @export

Stab_Single <- function(data, order.q=c(1,2), Alltime=TRUE, start_T=NULL, end_T=NULL){

  if(is.vector(data)){
    data <- matrix(data, nrow=1)
  }

  NA_num <- sum(is.na(data))
  if(NA_num!=0){
    stop('There are some NA in the data.')
  }
  if(sum(which(order.q<=0))){
    stop('Order q need larger than 0.')
  }
  # if(sum(apply(data,1,sum)==0)!=0){
  #   stop('There is at least one plot/community whose time series data are all zeros.')
  # }

  if(Alltime!=TRUE){
    if(is.null(start_T) | is.null(end_T)){
       stop('Need to set the length of time series for calculating.')
    }
    if((start_T>=end_T) | length(start_T)!=1 | length(end_T)!=1){
       stop('Starting and ending time need to be a number, and ending time needs larger than starting time.')
    }
  }

  stability <- function(vector,q){
    if(sum(vector!=0)==0){
      out <- NA
    }else{
      K <- length(vector)
      vector[which(vector==0)] <- 10^(-15)
      if(q==1){
        H <- sum((vector/sum(vector))*log(vector/sum(vector)))*(-1)
        out <- H/log(K)
      }else{
        up <- 1-sum((vector/sum(vector))^q)
        out <- up/(1-K^(1-q))
      }
    }
    return(out)
  }

  if(Alltime==TRUE){
    stab <- as.matrix(apply(data, 1, function(vec) sapply(order.q, function(qq) stability(vec, q=qq))))
  }else{
    subdata <- data[,c(start_T:end_T)]
    stab <- as.matrix(apply(subdata, 1, function(vec) sapply(order.q, function(qq) stability(vec, q=qq))))
  }
  result <- data.frame(Assemblage=rep(rownames(as.data.frame(data)), length(order.q)),
                       Order_q=rep(order.q, each=nrow(data)),
                       Stability=as.vector(t(stab)))
  colnames(result)[1] <- c("Plot/Community")
  return(result)
}





#' Calculate stability and synchrony of the time series data for multiple assemblages.
#'
#' \code{Stab_Syn_Multiple} is a function that calculate (Gamma, Alpha and Beta) stability and synchrony of the time series data (like biomass, productivity, etc.) for multiple assemblages.
#'
#' @param data can be input as a \code{data.frame/matrix} (assemblages by times), or a \code{list} of \code{data.frames} with each dataframe representing a assemblages-by-times data.
#' @param order.q a numerical vector specifying the orders of stability and synchrony. Default is c(1,2).
#' @param Alltime \code{TRUE} or \code{FALSE}, to decide whether to use all the times in (every) dataframe.
#' @param start_T (argument only for \code{Alltime = FALSE}) a positive integer specifying the start column of time in (every) dataframe.
#' @param end_T (argument only for \code{Alltime = FALSE}) a positive integer specifying the end column of time in (every) dataframe.
#'
#'
#' @return a dataframe with columns "Place", "Order_q", "Gamma", "Alpha", "Beta" and "Synchrony".
#'
#' @examples
#' # Stability of multiple plots
#' data("Jena_plot_biomass")
#' multiple_plot <- Jena_plot_biomass
#' output_multi_plot <- Stab_Syn_Multiple(data=multiple_plot, order.q=c(1,2), Alltime=TRUE)
#' output_multi_plot
#'
#' # Stability of multiple species in each plot
#' data("Jena_species_biomass")
#' multiple_species <- Jena_species_biomass
#' output_multi_species <- Stab_Syn_Multiple(data=multiple_species, order.q=c(1,2), Alltime=TRUE)
#' output_multi_species
#'
#'
#' @export

Stab_Syn_Multiple <- function(data, order.q=c(1,2), Alltime=TRUE, start_T=NULL, end_T=NULL){

  NA_num <- sum(is.na(data))
  if(NA_num!=0){
    stop('There are some NA in the data.')
  }
  if(sum(which(order.q<=0))){
    stop('Order q need larger than 0.')
  }

  # if(is.data.frame(data) | is.matrix(data)){
  #   if(nrow(data)<=1){
  #     stop('There is at least one list that has only a single assemblage.')
  #   }
  #   if(sum(which(apply(data,1,sum)==0))!=0){
  #     stop('There is at least one assemblage whose time series data are all zeros.')
  #   }
  # }else{
  #   num_row <- do.call(rbind, lapply(data, function(ZZ) nrow(ZZ)))
  #   if(sum(which(as.vector(num_row)<=1))!=0){
  #     stop('There is at least one list that has only a single assemblage.')
  #   }
  #   zero <- do.call(rbind, lapply(data, function(ZZ) sum(apply(ZZ,1,sum)==0)))
  #   if(sum(which(as.vector(zero)!=0))!=0){
  #     stop('There is at least one assemblage whose time series data are all zeros.')
  #   }
  # }

  if(Alltime!=TRUE){
    if(is.null(start_T) | is.null(end_T)){
      stop('Need to set the length of time series for calculating.')
    }
    if((start_T>=end_T) | length(start_T)!=1 | length(end_T)!=1){
      stop('Starting and ending time need to be a number, and ending time needs larger than starting time.')
    }
  }

  Stabillity_Multiple <- function(ZZ,q){
    K <- ncol(ZZ)
    #ZZ[is.na(ZZ)] <- 0

    z_iplus <- apply(ZZ,1,sum)
    if(length(which(z_iplus==0))!=0){
      ZZ <- ZZ[-which(z_iplus==0),]
      #z_iplus <- z_iplus[z_iplus!=0]
    }
    ZZ <- as.matrix(ZZ)
    ZZ[which(ZZ==0)] <- 10^(-15)
    z_iplus <- apply(ZZ,1,sum)

    z_plusk <- apply(ZZ,2,sum)
    z_plusplus <- sum(ZZ)
    # 0622_revise
    p_i <- as.data.frame(apply(ZZ,2,function(w)w/z_iplus))
    p_pool <- z_plusk/z_plusplus
    ww <- z_iplus/z_plusplus

    if(q==1){
      p_i_new <- p_i
      p_i_new[which(p_i_new==0)] <- 10^(-15)

      p_pool_new <- p_pool
      p_pool_new[which(p_pool_new==0)] <- 10^(-15)

      alpha <- (-1/log(K))*sum((ZZ/z_plusplus)*log(p_i_new))
      gamma <- (-1/log(K))*sum(p_pool*log(p_pool_new))
    }else{
      alpha <- (1-sum(apply(p_i,2,function(w)(w^q)*ww)))/(1-K^(1-q))
      gamma <- (1-sum(p_pool^q))/(1-K^(1-q))
    }

    return(c(gamma, alpha, (1-alpha)/(1-gamma), (1-alpha)-(1-gamma)))
  }

  #type="C"or"U"
  Synchrony <- function(ZZ,q,type="C"){
    M <- nrow(ZZ)
    K <- ncol(ZZ)
    #ZZ[is.na(ZZ)] <- 0

    if(M==1){
      value <- 0
    }else{
      z_iplus <- apply(ZZ,1,sum)
      if(length(which(z_iplus==0))!=0){
        ZZ <- ZZ[-which(z_iplus==0),]
        #z_iplus <- z_iplus[z_iplus!=0]
      }
      ZZ <- as.matrix(ZZ)
      ZZ[which(ZZ==0)] <- 10^(-15)
      z_iplus <- apply(ZZ,1,sum)

      if(length(z_iplus)<=1){
        value <- NA
      }else{
        z_plusk <- apply(ZZ,2,sum)
        z_plusplus <- sum(ZZ)
        p_i <- apply(ZZ,2,function(w)w/z_iplus)
        ww <- z_iplus/z_plusplus

        if(q==1){
          pool <- z_plusk/z_plusplus
          pool[which(pool==0)] <- 10^(-15)

          A <- sum(apply(p_i,2,function(w){
            nonzero <- which(w!=0)
            value <- sum(ww[nonzero]*w[nonzero]*log(w[nonzero]))
            return(value)
          })
          )
          G <- sum((z_plusk/z_plusplus)*log(z_plusk/z_plusplus))
          DOWN <- sum(ww*log(ww))*(-1)
          value <- (A-G)/DOWN

        }else{
          A <- sum(apply(p_i^q,2,function(w)w*ww))
          G <- sum((z_plusk/z_plusplus)^q)
          DOWN <- sum(apply(p_i^q,2,function(w)w*(ww-ww^q)))
          SUB <- sum((ZZ/z_plusplus)^q)/sum((z_plusk/z_plusplus)^q)

          if(type=="C"){
            value <- (A-G)/DOWN
          }else{
            value <- ((A-G)/DOWN)*SUB
          }
        }
      }
    }
    return(1-value)
  }


  if(is.data.frame(data) | is.matrix(data)){

    if(Alltime==TRUE){
      subdata <- data
    }else{
      subdata <- data[,c(start_T:end_T)]
    }
    out <- as.matrix(sapply(order.q, function(qq) c(Stabillity_Multiple(subdata,q=qq), Synchrony(subdata,q=qq))))
    result <- data.frame(Place=rep(1, length(order.q)),
                         Order_q=order.q, t(out))
    colnames(result)[3:7] <- c("Stab_Gamma", "Stab_Alpha", "Stab_Beta_multiplicative", "Stab_Beta_additive", "Synchrony")
    # revise_0605
    result <- result[,-5]

  }else{

    out <- lapply(order.q, function(qq){
                  cal <- lapply(data, function(ZZ){
                            if(Alltime==T){
                              subZZ <- ZZ
                            }else{
                              subZZ <- ZZ[,c(start_T:end_T)]
                            }
                            outout <- c(Stabillity_Multiple(subZZ,q=qq), Synchrony(subZZ,q=qq))
                            result <- data.frame(Order_q=qq, t(outout))
                            colnames(result)[2:6] <- c("Stab_Gamma", "Stab_Alpha", "Stab_Beta_multiplicative", "Stab_Beta_additive", "Synchrony")
                            # revise_0605
                            result <- result[,-4]
                            return(result)
                          })
                  cal2 <- do.call(rbind, cal)
                  calcal <- data.frame(Place=names(data), cal2)
                  return(calcal)
           })
    result <- do.call(rbind, out)
  }
  colnames(result) <- c("Place", "Order_q", "Gamma", "Alpha", "Beta", "Synchrony")
  return(result)

}


#' Calculate stability of the time series data for hierarchical structure.
#'
#' \code{Stab_Hier} is a function that calculate stability of the time series data (like biomass, productivity, etc.) for hierarchical structure.
#'
#' @param data can be input as \code{data.frame} (assemblages by times).
#' @param mat hierarchical structure of data.
#' @param order.q a numerical vector specifying the orders of stability. Default is c(1,2).
#' @param Alltime \code{TRUE} or \code{FALSE}, to decide whether to use all the times in the data.
#' @param start_T (argument only for \code{Alltime = FALSE}) a positive integer specifying the start column of time in the data.
#' @param end_T (argument only for \code{Alltime = FALSE}) a positive integer specifying the end column of time in the data.
#'
#' @import dplyr
#'
#' @return a dataframe with columns "Hier", "Order_q", and stability "Gamma", "Alpha" and "Beta (normalized)".
#'
#' @examples
#'
#' data("Jena_hierarchical_data")
#' data("Jena_hierarchical_mat")
#' output_hier <- Stab_Hier(data=Jena_hierarchical_data, mat=Jena_hierarchical_mat,
#'                          order.q=c(1,2), Alltime=TRUE)
#' output_hier
#'
#'
#' @export

Stab_Hier <- function(data, mat, order.q=c(1,2), Alltime=TRUE, start_T=NULL, end_T=NULL){
  if(Alltime==FALSE){
    data <- data[,start_T:end_T]
  }
  TT <- ncol(data)

  del <- which(apply(data,1,sum)==0)
  if(length(del)!=0){
    data <- data[-del,]
    mat <- mat[-del,]
  }

  data <- as.matrix(data)
  data[which(data==0)] <- 10^(-15)
  data <- as.data.frame(data)

  gamma_zz <- apply(data, 2, sum)
  gamma <- sapply(order.q, function(qq){
    if(qq==1){
      gamma_zz <- gamma_zz[gamma_zz!=0]
      gamma <- (-1)*sum((gamma_zz/sum(gamma_zz))*log((gamma_zz/sum(gamma_zz))))/log(TT)
    }else{
      gamma <- (1-sum((gamma_zz/sum(gamma_zz))^qq))/(1-TT^(1-qq))
    }
  })
  gamma <- data.frame(hier=ncol(mat)+1, order_q=order.q, gamma_value=gamma)


  hier_i_alphaanddiff <- function(data, mat, q, i){
    TT <- ncol(data)
    index <- unique(mat[,1:i])
    if(i==1){
      take <- sapply(1:length(index), function(rr)which(mat[,1:i]==index[rr]))
    }else{
      take <- sapply(1:nrow(index), function(rr)which(apply(mat[,1:i], 1, function(mm)sum(mm==index[rr,]))==i))
    }
    value <- lapply(take, function(tak){
      zz <- data[tak, ]
      zz <- apply(zz, 2, sum)
      outout <- sapply(q, function(qq){
        if(length(zz[zz!=0])!=0){
          if(qq==1){
            zz <- zz[zz!=0]
            out <- (-1)*sum((zz/sum(zz))*log((zz/sum(zz))))/log(TT)
          }else{
            out <- (1-sum((zz/sum(zz))^qq))/(1-TT^(1-qq))
          }
          ww <- sum(zz)/sum(data)
          vv <- ww*out
        }else{
          out <- 0
          ww <- 0
          vv <- ww*out
        }
        return(vv)
      })
      return(outout)
    })

    value <- do.call(rbind, value)
    if(length(q)==1){
      result <- sum(value)
    }else{
      result <- apply(value,2,sum)
    }

    ##diff
    if(i>=2){
      index2 <- unique(mat[,1:(i-1)])

      if(i==2){
        take2 <- sapply(1:length(index2), function(rr2){which(index[,1:(i-1)]==index2[rr2])})
      }else{
        take2 <- sapply(1:nrow(index2), function(rr2){which(apply(index[,1:(i-1)], 1, function(mm)sum(mm==index2[rr2,]))==(i-1))})
      }

      subdat <- lapply(take, function(tak){
        zz <- data[tak, ]
        zz <- apply(zz, 2, sum)
        return(zz)
      })

      if(i==2 & length(index2)<=2){

        subdat2 <- apply(take2, 2, function(tak2) {
          zz <- subdat[tak2]
          zz <- do.call(rbind, zz)
          vv <- sapply(q, function(qq) {
            if (qq == 1) {
              zz[which(zz == 0)] <- 10^(-15)
              value <- (-1)*sum((zz/sum(zz))*log(zz/sum(zz)))*(sum(zz)/sum(data))/log(TT)
            }else {
              value <- ((1 - sum((zz/sum(zz))^qq)) * (sum(zz)/sum(data)))/(1 - TT^(1 - qq))
            }
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
              value <- (-1)*sum((zz/sum(zz))*log(zz/sum(zz)))*(sum(zz)/sum(data))/log(TT)
            }else {
              value <- ((1 - sum((zz/sum(zz))^qq)) * (sum(zz)/sum(data)))/(1 - TT^(1 - qq))
            }
          })
          return(vv)
        })
        subdat2 <- do.call(rbind, subdat2)
      }

      if(length(q)==1){
        subdat2 <- sum(subdat2)
      }else{
        subdat2 <- apply(subdat2,2,sum)
      }

    }else{
      subdat2 <- lapply(take, function(tak){
        zz <- data[tak, ]
        zz <- apply(zz, 2, sum)
        vv <- sapply(q, function(qq){
          if(length(zz[zz!=0])!=0){
            if(qq==1){
              zz <- zz[zz!=0]
              out <- (-1)*sum((zz/sum(data))*log((zz/sum(data))))/log(TT)
            }else{
              out <- sum((zz/sum(data))^qq)
            }
          }else{
            out <- 0
          }
          return(out)
        })
        return(vv)
      })

      subdat2 <- do.call(rbind, subdat2)
      subdat2 <- sapply(1:length(q), function(ll){
                    if(q[ll]==1){
                      out <- sum(subdat2[,ll])
                    }else{
                      out <- (1-sum(subdat2[,ll]))/(1-TT^(1-q[ll]))
                    }
                  })
    }
    return(data.frame(hier=rep((ncol(mat)+1-i),length(q)), order_q=q, gamma_alpha=result, beta_max=subdat2-result))
  }

  stab <- c()
  for(ii in 1:ncol(mat)){
    stab <- rbind(stab, hier_i_alphaanddiff(data, mat, order.q, ii))
  }

  gamma_value <- rbind(gamma, gamma, data.frame(stab[,1:2], gamma_value=stab[,3]))
  gamma_value <- filter(gamma_value, hier!=1)
  gamma_value$hier <- rep(c((ncol(mat)+1):1), each=length(order.q))

  alpha_value <- c(rep(NA,length(order.q)), stab[,3])
  beta_max <- c(rep(NA,length(order.q)), stab[,4])
  alldata <- cbind(gamma_value, alpha_value, beta_max)

  alldata <- cbind(alldata[,1:4], beta_value=(alldata[,3]-alldata[,4])/alldata[,5])
  alldata <- as.data.frame(alldata)
  colnames(alldata) <- c("Hier","Order_q","Gamma","Alpha","Beta (normalized)")

  return(alldata)
}



#' ggplot2 extension for a Stab_Single, Stab_Syn_Multiple or Stab_Hier object with q-profile.
#'
#' \code{ggStab_Syn_qprofile} is a graphical function that based on the output from the function \code{Stab_Single}, \code{Stab_Syn_Multiple} or \code{Stab_Hier}. It provides to graph the q-profile of stability (and synchrony if is multiple assemblages).
#'
#' @param output the output obtained from \code{Stab_Single}, \code{Stab_Syn_Multiple} or \code{Stab_Hier}.
#'
#' @import ggpubr
#'
#' @return For a \code{Stab_Single} object, this function return a figure of q-profile for stability .
#' For a \code{Stab_Syn_Multiple} object, this function return a figure that contains q-profile for (Gamma, Alpha, Beta) stability and synchrony.
#' For a \code{Stab_Hier} object, this function return a figure that contains q-profile for gamma stability of highest hierarchical level and alpha stability of other hirarchical level, and also q-profile for normalized beta stability.
#'
#'
#' @examples
#' data("Jena_plot_biomass")
#' data("Jena_species_biomass")
#' data("Jena_hierarchical_data")
#' data("Jena_hierarchical_mat")
#'
#' ## Single assemblage
#' # Stability of each single plot
#' single_plot <- do.call(rbind, Jena_plot_biomass)
#' output_single_plot_q <- Stab_Single(data=single_plot[c(12,38),],
#'                                     order.q=seq(0.1,2,0.1), Alltime=TRUE)
#' ggStab_Syn_qprofile(output=output_single_plot_q)
#'
#' # Stability of each single species
#' single_species <- do.call(rbind, Jena_species_biomass)
#' output_single_species_q <- Stab_Single(data=single_species[c(40,49),],
#'                                        order.q=seq(0.1,2,0.1), Alltime=TRUE)
#' ggStab_Syn_qprofile(output=output_single_species_q)
#'
#'
#' ## Multiple assemblages
#' # Stability of multiple plots
#' multiple_plot <- Jena_plot_biomass
#' output_multi_plot_q <- Stab_Syn_Multiple(data=multiple_plot[c(9,11)],
#'                                          order.q=seq(0.1,2,0.1), Alltime=TRUE)
#' ggStab_Syn_qprofile(output=output_multi_plot_q)
#'
#' # Stability of multiple species in plot
#' multiple_species <- Jena_species_biomass
#' output_multi_species_q <- Stab_Syn_Multiple(data=multiple_species[c(62,70)],
#'                                             order.q=seq(0.1,2,0.1), Alltime=TRUE)
#' ggStab_Syn_qprofile(output=output_multi_species_q)
#'
#'
#' ## Hierarchies
#' output_hier_q <- Stab_Hier(data=Jena_hierarchical_data, mat=Jena_hierarchical_mat,
#'                            order.q=seq(0.1,2,0.1), Alltime=TRUE)
#' ggStab_Syn_qprofile(output=output_hier_q)
#'
#' @export

ggStab_Syn_qprofile <- function(output){

  if(length(which(colnames(output)=="Stability"))!=0){
    if(length(which(colnames(output)=="Plot/Community"))==0 | length(which(colnames(output)=="Order_q"))==0){
      stop('Please put the complete output of "Stab_Single", "Stab_Syn_Multiple" or "Stab_Hier" function.')
    }else{
      outtype <- "single"
    }
  }else if(length(which(colnames(output)=="Hier"))!=0){
    if(length(which(colnames(output)=="Order_q"))==0 | length(which(colnames(output)=="Gamma"))==0 | length(which(colnames(output)=="Alpha"))==0
       | length(which(colnames(output)=="Beta (normalized)"))==0){
      stop('Please put the complete output of "Stab_Single", "Stab_Syn_Multiple" or "Stab_Hier" function.')
    }else{
      outtype <- "hier"
    }
  }else{
    if(length(which(colnames(output)=="Place"))==0 | length(which(colnames(output)=="Order_q"))==0 | length(which(colnames(output)=="Gamma"))==0
       | length(which(colnames(output)=="Alpha"))==0 | length(which(colnames(output)=="Beta"))==0　
       | length(which(colnames(output)=="Synchrony"))==0){
      stop('Please put the complete output of "Stab_Single", "Stab_Syn_Multiple" or "Stab_Hier" function.')
    }else{
      outtype <- "multiple"
    }
  }

  if(outtype=="single"){
    output$`Plot/Community` <- factor(output$`Plot/Community`, levels=unique(output$`Plot/Community`))

    plotout <- ggplot(data=output, aes(x=Order_q, y=Stability, color=`Plot/Community`))+
                  geom_line(linewidth=1.2)+
                  ylab(label="Stability")+
                  xlab(label="Order of q")+
                  labs(color="Plot/Community")+ theme_bw()+
                  theme(legend.title = element_text(size=13), legend.text = element_text(size=12),
                        legend.key.size = unit(0.8, 'cm'), axis.title = element_text(size=16))


  }else if(outtype=="hier"){
    maxhier <- max(output$Hier)
    hier_num <- unique(output$Hier)
    qq <- unique(output$Order_q)
    type_name <- paste("S_alpha","(",hier_num[-1],")", sep="")
    type_diff <- paste("hier",hier_num[-1]+1,"-",hier_num[-1], sep="")

    plotdat1 <- data.frame(Order_q = rep(qq, length(hier_num)),
                           Stability = c(filter(output, Hier==maxhier)$Gamma, filter(output, Hier!=maxhier)$Alpha),
                           type = rep(c("S_gamma",type_name), each=length(qq)))
    plotdat1$type <- factor(plotdat1$type, levels=c("S_gamma",type_name))

    plotdat2 <- data.frame(Order_q = rep(qq, (length(hier_num)-1)),
                           Diff = filter(output, Hier!=maxhier)$`Beta (normalized)`,
                           type = rep(type_diff, each=length(qq)))
    plotdat2$type <- factor(plotdat2$type, levels=type_diff)

    plotout1 <- ggplot(data=plotdat1, aes(x=Order_q, y=Stability, color=type))+
                  geom_line(linewidth=1.2)+
                  ylab("Hierarchical Stability")+
                  labs(color="")+
                  theme_bw()+
                  theme(axis.text=element_text(size=10), axis.title=element_text(size=16),
                        plot.margin = unit(c(1,1,1,1), "cm"),
                        legend.key.size = unit(0.8, 'cm'),
                        legend.text = element_text(size=12),legend.position="bottom")

    if(length(type_name)+1>=4){
      plotout1 <- plotout1 + guides(color=guide_legend(nrow=2,byrow=TRUE))
    }

    plotout2 <- ggplot(data=plotdat2, aes(x=Order_q, y=Diff, color=type))+
                  geom_line(linewidth=1.2)+
                  ylab(label=expression(paste("Difference of Stability\n  (Normalized Beta)")))+
                  labs(color="")+
                  theme_bw()+
                  theme(axis.text=element_text(size=10), axis.title=element_text(size=16),
                        plot.margin = unit(c(1,1,1,1), "cm"),
                        legend.key.size = unit(0.8, 'cm'),
                        legend.text = element_text(size=12),legend.position="bottom")

    if(length(type_diff)>=4){
      plotout2 <- plotout2 + guides(color=guide_legend(nrow=2,byrow=TRUE))
    }


    plotout <- ggarrange(plotout1, plotout2, ncol = 2)


  }else{
    output$Place <- factor(output$Place, levels=unique(output$Place))

    # plotout <- list()
    stab_plotdat <- data.frame(Place=rep(output$Place,4),
                               Order_q=rep(output$Order_q,4),
                               value=c(output$Gamma,output$Alpha,output$Beta,output$Synchrony),
                               type=rep(c("Gamma","Alpha","Beta","Synchrony"),each=nrow(output)))
    stab_plotdat$type <- factor(stab_plotdat$type, levels=c("Gamma","Alpha","Beta","Synchrony"))

    plotout <- ggplot(data=stab_plotdat, aes(x=Order_q, y=value, color=Place))+
                    geom_line(linewidth=1.2)+
                    facet_wrap(.~type, nrow=2, scales = "free")+
                    ylab(label="Stability and Synchrony")+
                    xlab(label="Order of q")+
                    labs(color="Place")+ theme_bw()+
                    theme(strip.text = element_text(size=13), legend.title = element_text(size=13),
                          legend.text = element_text(size=12), legend.key.size = unit(0.8, 'cm'),
                          axis.title = element_text(size=16))
#
#     plotout[[2]] <- ggplot(data=output, aes(x=Order_q, y=Synchrony, color=Place))+
#                       geom_line(linewidth=1.2)+
#                       ylab(label="Synchrony")+
#                       xlab(label="Order of q")+
#                       labs(color="Place")+ theme_bw()+
#                       theme(legend.title = element_text(size=13), legend.text = element_text(size=12),
#                             legend.key.size = unit(0.8, 'cm'), axis.title = element_text(size=16))

  }
  return(plotout)
}



#' ggplot2 extension for a Stab_Single or Stab_Syn_Multiple object to analysis with an diversity (or other) variable.
#'
#' \code{ggStab_Syn_analysis} is a graphical function that based on the output from the function \code{Stab_Single} or \code{Stab_Syn_Multiple}. It provides to graph relationships between stability (and synchrony if is multiple assemblages) and an additional diversity (or other) variable .
#'
#' @param output the output obtained from \code{Stab_Single} or \code{Stab_Syn_Multiple} and needs to combine with a column that sets as \code{x_variable}. Also, if \code{by_group} is not \code{NULL}, the output also need to combine with the column that sets as \code{by_group}.
#' @param x_variable name of the column of diversity (or other) variable, that will use as the x-axis in the plot.
#' @param by_group name of the column that is a categorical variable for plotting points with different color. And it is required if \code{model = "LMM"}, model uses it as random effect for intercept and slope. Default is \code{NULL}.
#' @param model specifying the fitting model, \code{model = "lm"} for linear model; \code{model = "LMM"} for linear mixed model with random effects for intercept and slope. Default is \code{model = "LMM"}.
#'
#'
#' @import ggplot2
#' @import stringr
#' @import lme4
#' @import lmerTest
#' @import dplyr
#'
#'
#' @return For an \code{Stab_Single} object, this function return a figure of diversity (or other) variable vs. stability.
#' For an \code{Stab_Syn_Multiple} object, this function return a figure that is about diversity (or other) variable vs. (Gamma, Alpha, Beta) stability and synchrony.
#'
#' @examples
#' data("Jena_plot_biomass")
#' data("Jena_species_biomass")
#' data("Jena_hierarchical_data")
#' data("Jena_hierarchical_mat")
#'
#' ## Single assemblage
#' # Stability of each single plot
#' single_plot <- do.call(rbind, Jena_plot_biomass)
#' output_single_plot_div <- Stab_Single(data=single_plot, order.q=c(1,2), Alltime=TRUE)
#' output_single_plot_div <- data.frame(output_single_plot_div,
#'                                      sowndiv=as.numeric(do.call(rbind,
#'                                        strsplit(output_single_plot_div[,1],"[._]+"))[,2]),
#'                                      block=do.call(rbind,
#'                                        strsplit(output_single_plot_div[,1],"[._]+"))[,1])
#'
#' ggStab_Syn_analysis(output=output_single_plot_div, x_variable="sowndiv",
#'                     by_group="block", model="LMM")
#'
#' # Stability of each single species
#' single_species <- do.call(rbind, Jena_species_biomass)
#' output_single_species_div <- Stab_Single(data=single_species,
#'                                          order.q=c(1,2), Alltime=TRUE)
#' output_single_species_div <- data.frame(output_single_species_div,
#'                               sowndiv=as.numeric(do.call(rbind,
#'                                       strsplit(output_single_species_div[,1],"[._]+"))[,3]),
#'                               block=do.call(rbind,
#'                                     strsplit(output_single_species_div[,1],"[._]+"))[,2])
#'
#' ggStab_Syn_analysis(output=output_single_species_div, x_variable="sowndiv",
#'                     by_group="block", model="LMM")
#'
#'
#' ## Multiple assemblages
#' # Stability of multiple plots
#' multiple_plot <- Jena_plot_biomass
#' output_multi_plot_div <- Stab_Syn_Multiple(data=multiple_plot, order.q=c(1,2), Alltime=TRUE)
#' output_multi_plot_div <- data.frame(output_multi_plot_div, sowndiv=rep(c(16,8,4,2,1),8),
#'                                     block=rep(rep(c("B1","B2","B3","B4"),each=5),2))
#'
#' ggStab_Syn_analysis(output=output_multi_plot_div, x_variable="sowndiv",
#'                     by_group="block", model="LMM")
#'
#' # Stability of multiple species in plot
#' multiple_species <- Jena_species_biomass
#' output_multi_species_div <- Stab_Syn_Multiple(data=multiple_species,
#'                                               order.q=c(1,2), Alltime=TRUE)
#' output_multi_species_div <- data.frame(output_multi_species_div,
#'                              sowndiv=as.numeric(do.call(rbind,
#'                                      strsplit(output_multi_species_div[,1],"_"))[,3]),
#'                              block=do.call(rbind,
#'                                    strsplit(output_multi_species_div[,1],"_"))[,2])
#'
#' ggStab_Syn_analysis(output=output_multi_species_div, x_variable="sowndiv",
#'                     by_group="block", model="LMM")
#'
#' @export

ggStab_Syn_analysis <- function(output, x_variable, by_group=NULL, model="LMM"){

  # x_variable_ori <- x_variable
  # if(is.null(by_group)==FALSE){by_group_ori <- by_group}

  if(length(which(colnames(output)=="Stability"))!=0){
    if((length(which(colnames(output)=="Plot/Community"))==0 & length(which(colnames(output)=="Plot.Community"))==0) | length(which(colnames(output)=="Order_q"))==0){
      stop('Please put the complete output of "Stab_Single" or "Stab_Syn_Multiple" function.')
    }
  }else{
    if(length(which(colnames(output)=="Place"))==0 | length(which(colnames(output)=="Order_q"))==0 | length(which(colnames(output)=="Gamma"))==0
       | length(which(colnames(output)=="Alpha"))==0 | length(which(colnames(output)=="Beta"))==0　
       | length(which(colnames(output)=="Synchrony"))==0){
      stop('Please put the complete output of "Stab_Single" or "Stab_Syn_Multiple" function.')
    }
  }
  if(length(which(colnames(output)==x_variable))==0){
    stop('The output data need to combine a column setting as x_variable.')
  }else{
    colnames(output)[which(colnames(output)==x_variable)] <- c("Xvariable")
  }
  if(is.null(by_group)==FALSE){
    if(length(which(colnames(output)==by_group))==0){
      stop('The output data need to combine a column setting as by_group.')
    }else{
      colnames(output)[which(colnames(output)==by_group)] <- c("Gvariable")
    }
  }
  if(model=="LMM"){
    if(is.null(by_group)==TRUE){
      stop('For linear mixed model, you need to set by_group variable as random effect.')
    }
  }

  if(is.null(by_group)==FALSE){
    output$Gvariable <- as.factor(output$Gvariable)
  }

  if(length(which(colnames(output)=="Stability"))!=0){

    # fit model
    lm_sign <-c()
    lm_slope <- c()
    plotdata <- c()
    for(qq in unique(output$Order_q)){
      subdata <- filter(output, Order_q==qq)
      if(model=="lm"){
        MODEL <- lm(Stability ~ Xvariable, subdata)
        summary <- summary(MODEL)
        sign <- summary$coefficients[2,4]
        pred_value <- predict(MODEL, newdata=subdata)
      }else{
        MODEL <- lmer(Stability ~ 1 + Xvariable + (1 + Xvariable | Gvariable), subdata)
        summary <- summary(MODEL)
        sign <- summary$coefficients[2,5]
        pred_value <- predict(MODEL, newdata=subdata, re.form=NA)
      }

      lm_sign <- rbind(lm_sign, c(sign=ifelse(sign<0.05, "significant", "non-significant")))
      lm_slope <- rbind(lm_slope, c(slope=summary$coefficients[2,1]))
      plotdata <- rbind(plotdata, data.frame(subdata, pred=pred_value))
    }
    lm_sign <- as.data.frame(lm_sign)
    lm_slope <- as.data.frame(lm_slope)
    rownames(lm_sign) <- paste("q = ", unique(output$Order_q), sep="")
    rownames(lm_slope) <- paste("q = ", unique(output$Order_q), sep="")

    plotdata$Order_q <- paste("q = ", plotdata$Order_q, sep="")

    plotdata$sign <- sapply(plotdata$Order_q, function(yy){
                        sign <- lm_sign[which(rownames(lm_sign)==yy),1]
                        return(sign)
                     })
    plotdata$sign <- factor(plotdata$sign, levels=c("significant", "non-significant"))

    slope_text <- data.frame(slope = paste("slope = ",round(lm_slope[,1],4),sep=""),
                             Order_q = rownames(lm_slope))

    if(is.null(by_group)==FALSE){

      plotout <- ggplot(plotdata, aes(x=Xvariable, y=Stability))+
                    geom_point(aes(color=Gvariable), size=2.7)+
                    geom_line(aes(x=Xvariable, y=pred, linetype = sign), linewidth=1.2, color="black")+
                    scale_linetype_manual(values=c("solid","dashed"), drop = FALSE)+
                    facet_wrap(.~Order_q, scales="fixed", ncol = min(5, length(unique(plotdata$Order_q))))+
                    labs(linetype="", color=by_group)+
                    xlab(label=x_variable)+
                    ylab(label="Stability")+ theme_bw()+
                    geom_text(data=slope_text, aes(x = -Inf, y = -Inf, label = slope),
                              x=max(plotdata$Xvariable, na.rm = TRUE),
                              y=min(plotdata$Stability, na.rm = TRUE), size=4.5, hjust=1.1, vjust=0.1)+
                    theme(strip.text = element_text(size=13), legend.title = element_text(size=13),
                          legend.text = element_text(size=12), legend.key.size = unit(0.8, 'cm'),
                          axis.title = element_text(size=16))

    }else{
      plotout <- ggplot(plotdata, aes(x=Xvariable, y=Stability))+
                    geom_point(size=2.7)+
                    geom_line(aes(x=Xvariable, y=pred, linetype = sign), linewidth=1.2, color="black")+
                    scale_linetype_manual(values=c("solid","dashed"), drop = FALSE)+
                    facet_wrap(.~Order_q, scales="fixed", ncol = min(5, length(unique(plotdata$Order_q))))+
                    labs(linetype="")+
                    xlab(label=x_variable)+
                    ylab(label="Stability")+ theme_bw()+
                    geom_text(data=slope_text, aes(x = -Inf, y = -Inf, label = slope),
                              x=max(plotdata$Xvariable, na.rm = TRUE),
                              y=min(plotdata$Stability, na.rm = TRUE), size=5, hjust=1, vjust=0.1)+
                    theme(strip.text = element_text(size=13), legend.title = element_text(size=13),
                          legend.text = element_text(size=12), legend.key.size = unit(0.8, 'cm'),
                          axis.title = element_text(size=16))
    }

  }else{

    # fit model
    lm_sign <-c()
    lm_slope <- c()
    plotdata <- c()
    for(qq in unique(output$Order_q)){
      subdata <- filter(output, Order_q==qq)
      if(model=="lm"){
        sign_num <- 4

        MODEL_G <- lm(Gamma ~ Xvariable, subdata)
        sum_gamma <- summary(MODEL_G)
        pred_G <- predict(MODEL_G, newdata=subdata)
        MODEL_A <- lm(Alpha ~ Xvariable, subdata)
        sum_alpha <- summary(MODEL_A)
        pred_A <- predict(MODEL_A, newdata=subdata)
        # MODEL_BM <- lm(Stab_Beta_multiplicative ~ Xvariable, subdata)
        # sum_beta_multi <- summary(MODEL_BM)
        # pred_BM <- predict(MODEL_BM, newdata=subdata)
        MODEL_BA <- lm(Beta ~ Xvariable, subdata)
        sum_beta_add <- summary(MODEL_BA)
        pred_BA <- predict(MODEL_BA, newdata=subdata)
        MODEL_S <- lm(Synchrony ~ Xvariable, subdata)
        sum_syn <- summary(MODEL_S)
        pred_S <- predict(MODEL_S, newdata=subdata)

      }else{
        sign_num <- 5

        MODEL_G <- lmer(Gamma ~ 1 + Xvariable + (1 + Xvariable | Gvariable), subdata)
        sum_gamma <- summary(MODEL_G)
        pred_G <- predict(MODEL_G, newdata=subdata, re.form=NA)
        MODEL_A <- lmer(Alpha ~ 1 + Xvariable + (1 + Xvariable | Gvariable), subdata)
        sum_alpha <- summary(MODEL_A)
        pred_A <- predict(MODEL_A, newdata=subdata, re.form=NA)
        # MODEL_BM <- lmer(Stab_Beta_multiplicative ~ 1 + Xvariable + (1 + Xvariable | Gvariable), subdata)
        # sum_beta_multi <- summary(MODEL_BM)
        # pred_BM <- predict(MODEL_BM, newdata=subdata, re.form=NA)
        MODEL_BA <- lmer(Beta ~ 1 + Xvariable + (1 + Xvariable | Gvariable), subdata)
        sum_beta_add <- summary(MODEL_BA)
        pred_BA <- predict(MODEL_BA, newdata=subdata, re.form=NA)
        MODEL_S <- lmer(Synchrony ~ 1 + Xvariable + (1 + Xvariable | Gvariable), subdata)
        sum_syn <- summary(MODEL_S)
        pred_S <- predict(MODEL_S, newdata=subdata, re.form=NA)
      }

      lm_sign <- rbind(lm_sign,
                       c(gamma=ifelse(sum_gamma$coefficients[2,sign_num]<0.05, "significant", "non-significant"),
                         alpha=ifelse(sum_alpha$coefficients[2,sign_num]<0.05, "significant", "non-significant"),
                         beta_add=ifelse(sum_beta_add$coefficients[2,sign_num]<0.05, "significant", "non-significant"),
                         synchrony=ifelse(sum_syn$coefficients[2,sign_num]<0.05, "significant", "non-significant")))
      lm_slope <- rbind(lm_slope,
                        c(gamma=sum_gamma$coefficients[2,1],
                          alpha=sum_alpha$coefficients[2,1],
                          beta_add=sum_beta_add$coefficients[2,1],
                          synchrony=sum_syn$coefficients[2,1]))
      plotdata <- rbind(plotdata, data.frame(subdata, pred_G=pred_G, pred_A=pred_A, pred_BA=pred_BA, pred_S=pred_S))
    }
    lm_sign <- as.data.frame(lm_sign)
    lm_slope <- as.data.frame(lm_slope)
    rownames(lm_sign) <- paste("q = ", unique(output$Order_q), sep="")
    rownames(lm_slope) <- paste("q = ", unique(output$Order_q), sep="")

    plotdata$Order_q <- paste("q = ", plotdata$Order_q, sep="")

    plotdata$sign_G <- sapply(plotdata$Order_q, function(yy){
      lm_sign[which(rownames(lm_sign)==yy),1]
    })
    plotdata$sign_A <- sapply(plotdata$Order_q, function(yy){
      lm_sign[which(rownames(lm_sign)==yy),2]
    })
    plotdata$sign_BA <- sapply(plotdata$Order_q, function(yy){
      lm_sign[which(rownames(lm_sign)==yy),3]
    })
    plotdata$sign_S <- sapply(plotdata$Order_q, function(yy){
      lm_sign[which(rownames(lm_sign)==yy),4]
    })

    plotdata_Stab <- data.frame(Place = rep(plotdata$Place,4),
                                Order_q = rep(plotdata$Order_q,4),
                                Stability = c(plotdata$Gamma, plotdata$Alpha, plotdata$Beta, plotdata$Synchrony),
                                pred = c(plotdata$pred_G, plotdata$pred_A, plotdata$pred_BA, plotdata$pred_S),
                                sign = c(plotdata$sign_G, plotdata$sign_A, plotdata$sign_BA, plotdata$sign_S),
                                type = rep(c("Gamma","Alpha","Beta","Synchrony"), each = nrow(plotdata)),
                                Xvariable = rep(plotdata$Xvariable,4))

    if(is.null(by_group)==FALSE){
      plotdata_Stab$Gvariable <- rep(plotdata$Gvariable,4)
      plotdata_Stab$Gvariable <- as.factor(plotdata_Stab$Gvariable)
      plotdata$Gvariable <- as.factor(plotdata$Gvariable)
    }

    plotdata_Stab$sign <- factor(plotdata_Stab$sign, levels=c("significant", "non-significant"))
    plotdata_Stab$type <- factor(plotdata_Stab$type, levels = c("Gamma","Alpha","Beta","Synchrony"))
    plotdata_Stab$Order_q <- as.factor(plotdata_Stab$Order_q)


    # plotdata$sign_S <- factor(plotdata$sign_S, levels=c("significant", "non-significant"))
    # plotdata$Order_q <- as.factor(plotdata$Order_q)


    slope_text_Stab <- data.frame(slope = paste("slope = ",round(as.vector(as.matrix(lm_slope)),4),sep=""),
                                  Order_q = rep(rownames(lm_slope),4),
                                  type = rep(c("Gamma","Alpha","Beta","Synchrony"), each = nrow(lm_slope)))
    slope_text_Stab$type <- factor(slope_text_Stab$type, levels = c("Gamma","Alpha","Beta","Synchrony"))
    slope_text_Stab$Order_q <- as.factor(slope_text_Stab$Order_q)

    # slope_text_Syn <- data.frame(slope = paste("slope = ",round(lm_slope[,5],4),sep=""),
    #                              Order_q = rownames(lm_slope))
    # slope_text_Syn$Order_q <- as.factor(slope_text_Syn$Order_q)

    tyG <- min(filter(plotdata_Stab, type=="Gamma")$Stability, na.rm = TRUE)
    tyA <- min(filter(plotdata_Stab, type=="Alpha")$Stability, na.rm = TRUE)
    # tyBM <- max(filter(plotdata_Stab, type=="Beta (multiplicative)")$Stability)
    tyBA <- max(filter(plotdata_Stab, type=="Beta")$Stability, na.rm = TRUE)
    tyS <- min(filter(plotdata_Stab, type=="Synchrony")$Stability, na.rm = TRUE)

    if(is.null(by_group)==FALSE){

      plotout1 <- ggplot(plotdata_Stab, aes(x=Xvariable, y=Stability))+
                    geom_point(aes(color=Gvariable), size=2.7)+
                    geom_line(aes(x=Xvariable, y=pred, linetype = sign), linewidth=1.2, color="black")+
                    scale_linetype_manual(values=c("solid","dashed"), drop = FALSE)+
                    facet_grid(type~Order_q, scales = "free_y")+
                    labs(linetype="", color=by_group)+
                    xlab(label=x_variable)+
                    ylab(label="Stability and Synchrony")+ theme_bw()+
                    geom_text(data=slope_text_Stab, aes(x = -Inf, y = -Inf, label = slope),
                              x=max(plotdata_Stab$Xvariable, na.rm = TRUE),
                              y=rep(c(tyG, tyA, tyBA, tyS),each=2), size=5,
                              hjust=1, vjust=rep(c(0.1,0.1,1.1,0.1),each=2))+
                    theme(strip.text = element_text(size=13), legend.title = element_text(size=13),
                          legend.text = element_text(size=12), legend.key.size = unit(0.8, 'cm'),
                          axis.title = element_text(size=16))

      # plotout2 <- ggplot(plotdata, aes(x=Xvariable, y=Synchrony))+
      #               geom_point(aes(color=Gvariable), size=2.7)+
      #               geom_line(aes(x=Xvariable, y=pred_S, linetype = sign_S), linewidth=1.2, color="black")+
      #               scale_linetype_manual(values=c("solid","dashed"), drop = FALSE)+
      #               facet_wrap(.~Order_q, scales = "fixed")+
      #               labs(linetype="", color=by_group)+
      #               xlab(label=x_variable)+
      #               ylab(label="Synchrony")+ theme_bw()+
      #               geom_text(data=slope_text_Syn, aes(x = -Inf, y = -Inf, label = slope),
      #                         x=max(plotdata$Xvariable),
      #                         y=min(plotdata$Synchrony), size=5,
      #                         hjust=1, vjust=0.1)+
      #               theme(strip.text = element_text(size=13), legend.title = element_text(size=13),
      #                     legend.text = element_text(size=12), legend.key.size = unit(0.8, 'cm'),
      #                     axis.title = element_text(size=16))


    }else{

      plotout1 <- ggplot(plotdata_Stab, aes(x=Xvariable, y=Stability))+
                      geom_point(size=2.7)+
                      geom_line(aes(x=Xvariable, y=pred, linetype = sign), linewidth=1.2, color="black")+
                      scale_linetype_manual(values=c("solid","dashed"), drop = FALSE)+
                      facet_grid(type~Order_q, scales = "free_y")+
                      labs(linetype="")+
                      xlab(label=x_variable)+
                      ylab(label="Stability and Synchrony")+ theme_bw()+
                      geom_text(data=slope_text_Stab, aes(x = -Inf, y = -Inf, label = slope),
                                x=max(plotdata_Stab$Xvariable, na.rm = TRUE),
                                y=rep(c(tyG, tyA, tyBA, tyS),each=2), size=5,
                                hjust=1, vjust=rep(c(0.1,1.1),each=4))+
                      theme(strip.text = element_text(size=13), legend.title = element_text(size=13),
                            legend.text = element_text(size=12), legend.key.size = unit(0.8, 'cm'),
                            axis.title = element_text(size=16))

      # plotout2 <- ggplot(plotdata, aes(x=Xvariable, y=Synchrony))+
      #                 geom_point(size=2.7)+
      #                 geom_line(aes(x=Xvariable, y=pred_S, linetype = sign_S), linewidth=1.2, color="black")+
      #                 scale_linetype_manual(values=c("solid","dashed"), drop = FALSE)+
      #                 facet_wrap(.~Order_q, scales = "fixed")+
      #                 labs(linetype="")+
      #                 xlab(label=x_variable)+
      #                 ylab(label="Synchrony")+ theme_bw()+
      #                 geom_text(data=slope_text_Syn, aes(x = -Inf, y = -Inf, label = slope),
      #                           x=max(plotdata$Xvariable),
      #                           y=min(plotdata$Synchrony), size=5,
      #                           hjust=1, vjust=0.1)+
      #                 theme(strip.text = element_text(size=13), legend.title = element_text(size=13),
      #                       legend.text = element_text(size=12), legend.key.size = unit(0.8, 'cm'),
      #                       axis.title = element_text(size=16))

    }
    plotout <- plotout1
    # plotout[[1]] <- plotout1
    # plotout[[2]] <- plotout2
  }

  return(plotout)
}






## ========== no visible global function definition for R CMD check ========== ##
utils::globalVariables(c("COLOR", "Order_q", "Place", "Stability", "Synchrony"
))




