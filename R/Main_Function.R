#' Calculate stability of the time series data for single assemblage.
#'
#' \code{Stay_Single} is a function that calculate stability of the time series data (like biomass, productivity, etc.) for single assemblage.
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
#' data("Jena_plot_biomass_data")
#' single_plot <- do.call(rbind, Jena_plot_biomass_data)
#' output_single_plot <- Stay_Single(data=single_plot, order.q=c(1,2), Alltime=TRUE)
#' output_single_plot
#'
#' # Stability of each single species in each plot
#' data("Jena_species_biomass_data")
#' single_species <- do.call(rbind, Jena_species_biomass_data)
#' output_single_species <- Stay_Single(data=single_species, order.q=c(1,2), Alltime=TRUE)
#' output_single_species
#'
#' @export

Stay_Single <- function(data, order.q=c(1,2), Alltime=TRUE, start_T=NULL, end_T=NULL){

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
#' \code{Stay_Multiple} is a function that calculate (Gamma, Alpha and Beta) stability and synchrony of the time series data (like biomass, productivity, etc.) for multiple assemblages.
#'
#' @param data can be input as a \code{data.frame/matrix} (assemblages by times), or a \code{list} of \code{data.frames} with each dataframe representing a assemblages-by-times data.
#' @param order.q a numerical vector specifying the orders of stability and synchrony. Default is c(1,2).
#' @param Alltime \code{TRUE} or \code{FALSE}, to decide whether to use all the times in (every) dataframe.
#' @param start_T (argument only for \code{Alltime = FALSE}) a positive integer specifying the start column of time in (every) dataframe.
#' @param end_T (argument only for \code{Alltime = FALSE}) a positive integer specifying the end column of time in (every) dataframe.
#'
#'
#' @return a dataframe with columns "Site", "Order_q", "Gamma", "Alpha", "Beta" and "Synchrony".
#'
#' @examples
#' # Stability of multiple plots
#' data("Jena_plot_biomass_data")
#' multiple_plot <- Jena_plot_biomass_data
#' output_multi_plot <- Stay_Multiple(data=multiple_plot, order.q=c(1,2), Alltime=TRUE)
#' output_multi_plot
#'
#' # Stability of multiple species in each plot
#' data("Jena_species_biomass_data")
#' multiple_species <- Jena_species_biomass_data
#' output_multi_species <- Stay_Multiple(data=multiple_species, order.q=c(1,2), Alltime=TRUE)
#' output_multi_species
#'
#'
#' @export

Stay_Multiple <- function(data, order.q=c(1,2), Alltime=TRUE, start_T=NULL, end_T=NULL){

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
    result <- data.frame(Site=rep(1, length(order.q)),
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
                  calcal <- data.frame(Site=names(data), cal2)
                  return(calcal)
           })
    result <- do.call(rbind, out)
  }
  colnames(result) <- c("Site", "Order_q", "Gamma", "Alpha", "Beta", "Synchrony")
  return(result)

}


#' Calculate stability of the time series data for hierarchical structure.
#'
#' \code{Stay_Hier} is a function that calculate stability of the time series data (like biomass, productivity, etc.) for hierarchical structure.
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
#' @return a dataframe with columns "Hier", "Order_q", and stability "Gamma", "Alpha", "Beta" and "Synchrony".
#'
#' @examples
#'
#' data("Jena_hierarchical_biomass_data")
#' data("Jena_hierarchical_mat")
#' output_hier <- Stay_Hier(data=Jena_hierarchical_biomass_data, mat=Jena_hierarchical_mat,
#'                          order.q=c(1,2), Alltime=TRUE)
#' output_hier
#'
#'
#' @export

Stay_Hier <- function(data, mat, order.q=c(1,2), Alltime=TRUE, start_T=NULL, end_T=NULL){
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

  alldata <- cbind(alldata[,1:4], beta=alldata[,3]-alldata[,4],
                   synchrony=1-((alldata[,3]-alldata[,4])/alldata[,5]))
  alldata <- as.data.frame(alldata)
  colnames(alldata) <- c("Hier","Order_q","Gamma","Alpha","Beta","Synchrony")

  return(alldata)
}



#' ggplot2 extension for a Stay_Single, Stay_Multiple or Stay_Hier object with q-profile.
#'
#' \code{ggStay_qprofile} is a graphical function that based on the output from the function \code{Stay_Single}, \code{Stay_Multiple} or \code{Stay_Hier}. It provides to graph the q-profile of stability (and synchrony if is multiple assemblages).
#'
#' @param output the output obtained from \code{Stay_Single}, \code{Stay_Multiple} or \code{Stay_Hier}.
#'
#' @import ggpubr
#'
#'
#' @return For a \code{Stay_Single} object, this function return a figure of q-profile for stability .
#' For a \code{Stay_Multiple} object, this function return a figure that contains q-profile for (Gamma, Alpha, Beta) stability and synchrony.
#' For a \code{Stay_Hier} object, this function return a figure that contains q-profile for gamma stability of highest hierarchical level and alpha stability of other hierarchical level. And it also provides a figure about the relationship between each decomposition of overall stability and the synchrony of each hierarchical level and order.q.
#'
#'
#' @examples
#' data("Jena_plot_biomass_data")
#' data("Jena_species_biomass_data")
#' data("Jena_hierarchical_biomass_data")
#' data("Jena_hierarchical_mat")
#'
#' ## Single assemblage
#' # Stability of each single plot
#' single_plot <- do.call(rbind, Jena_plot_biomass_data)
#' output_single_plot_q <- Stay_Single(data=single_plot[c(12,38),],
#'                                     order.q=seq(0.1,2,0.1), Alltime=TRUE)
#' ggStay_qprofile(output=output_single_plot_q)
#'
#' # Stability of each single species
#' single_species <- do.call(rbind, Jena_species_biomass_data)
#' output_single_species_q <- Stay_Single(data=single_species[c(40,49),],
#'                                        order.q=seq(0.1,2,0.1), Alltime=TRUE)
#' ggStay_qprofile(output=output_single_species_q)
#'
#'
#' ## Multiple assemblages
#' # Stability of multiple plots
#' multiple_plot <- Jena_plot_biomass_data
#' output_multi_plot_q <- Stay_Multiple(data=multiple_plot[c(9,11)],
#'                                          order.q=seq(0.1,2,0.1), Alltime=TRUE)
#' ggStay_qprofile(output=output_multi_plot_q)
#'
#' # Stability of multiple species in plot
#' multiple_species <- Jena_species_biomass_data
#' output_multi_species_q <- Stay_Multiple(data=multiple_species[c(62,70)],
#'                                             order.q=seq(0.1,2,0.1), Alltime=TRUE)
#' ggStay_qprofile(output=output_multi_species_q)
#'
#'
#' ## Hierarchies
#' output_hier_q <- Stay_Hier(data=Jena_hierarchical_biomass_data,
#'                            mat=Jena_hierarchical_mat,
#'                            order.q=seq(0.1,2,0.1), Alltime=TRUE)
#' ggStay_qprofile(output=output_hier_q)
#'
#' @export

ggStay_qprofile <- function(output){

  if(length(which(colnames(output)=="Stability"))!=0){
    if(length(which(colnames(output)=="Plot/Community"))==0 | length(which(colnames(output)=="Order_q"))==0){
      stop('Please put the complete output of "Stay_Single", "Stay_Multiple" or "Stay_Hier" function.')
    }else{
      outtype <- "single"
    }
  }else if(length(which(colnames(output)=="Hier"))!=0){
    if(length(which(colnames(output)=="Order_q"))==0 | length(which(colnames(output)=="Gamma"))==0 | length(which(colnames(output)=="Alpha"))==0
       | length(which(colnames(output)=="Beta"))==0 | length(which(colnames(output)=="Synchrony"))==0){
      stop('Please put the complete output of "Stay_Single", "Stay_Multiple" or "Stay_Hier" function.')
    }else{
      outtype <- "hier"
    }
  }else{
    if(length(which(colnames(output)=="Site"))==0 | length(which(colnames(output)=="Order_q"))==0 | length(which(colnames(output)=="Gamma"))==0
       | length(which(colnames(output)=="Alpha"))==0 | length(which(colnames(output)=="Beta"))==0　
       | length(which(colnames(output)=="Synchrony"))==0){
      stop('Please put the complete output of "Stay_Single", "Stay_Multiple" or "Stay_Hier" function.')
    }else{
      outtype <- "multiple"
    }
  }


  if(outtype=="single"){

    # Check if the number of unique 'Assemblage' is 4 or less
    if (length(unique(output$`Plot/Community`)) <= 4){
      cbPalette <- c("#EA0000","#0066CC","#64A600","#FFA042")
    }else{
      # If there are more than 4 assemblages, start with the same predefined color palette
      # Then extend the palette by generating additional colors using the 'ggplotColors' function
      cbPalette <- c(c("#EA0000","#0066CC","#64A600","#FFA042"),
                     ggplotColors(length(unique(output$`Plot/Community`))-4))
    }

    output$`Plot/Community` <- factor(output$`Plot/Community`, levels=unique(output$`Plot/Community`))

    plotout <- ggplot(data=output, aes(x=Order_q, y=Stability, color=`Plot/Community`))+
      geom_line(linewidth=1.2)+
      ylab(label="Stability")+
      xlab(label="Order of q")+
      labs(color="Plot/Community")+ theme_bw()+
      scale_colour_manual(values = cbPalette) +
      theme(legend.title = element_text(size=13), legend.text = element_text(size=12),
            legend.key.size = unit(0.8, 'cm'), axis.title = element_text(size=16))


  }else if(outtype=="hier"){
    maxhier <- max(output$Hier)
    hier_num <- unique(output$Hier)
    qq <- unique(output$Order_q)
    type_name <- c(paste("Gamma(",maxhier,")",sep=""),paste("Alpha","(",hier_num[-1],")", sep=""))
    type_name2_2 <- c(paste("Beta","(",hier_num[-1],")", sep=""), "Alpha(1)")
    type_diff <- paste("Hier",hier_num[-1], sep="")

    plotdat1 <- data.frame(Order_q = rep(qq, length(hier_num)),
                           Stability = c(filter(output, Hier==maxhier)$Gamma, filter(output, Hier!=maxhier)$Alpha),
                           type = rep(type_name, each=length(qq)))
    plotdat1$type <- factor(plotdat1$type, levels=type_name)

    plotdat2 <- data.frame(Order_q = rep(qq, (length(hier_num)-1)),
                           Synchrony = filter(output, Hier!=maxhier)$Synchrony,
                           type = rep(type_diff, each=length(qq)))
    plotdat2$type <- factor(plotdat2$type, levels=type_diff)

    plotdat2_2 <- data.frame(Order_q = rep(qq, length(hier_num)),
                             Stability = c(filter(output, Hier!=maxhier)$Beta, filter(output, Hier==1)$Alpha),
                             type = rep(type_name2_2, each=length(qq)))
    plotdat2_2$type <- factor(plotdat2_2$type, levels=type_name2_2)


    # Check if the number of unique 'Assemblage' is 4 or less
    if (length(unique(plotdat1$type)) <= 4){
      cbPalette <- c("#EA0000","#0066CC","#64A600","#FFA042")
    }else{
      cbPalette <- c(c("#EA0000","#0066CC","#64A600","#FFA042"),
                     ggplotColors(length(unique(plotdat1$type))-4))
    }

    if (length(unique(plotdat2$type)) <= 4){
      cbPalette_2 <- c("#EA0000","#0066CC","#64A600","#FFA042")
    }else{
      cbPalette_2 <- c(c("#EA0000","#0066CC","#64A600","#FFA042"),
                     ggplotColors(length(unique(plotdat2$type))-4))
    }

    plotout1 <- ggplot(data=plotdat1, aes(x=Order_q, y=Stability, color=type))+
                  geom_line(linewidth=1.2)+
                  ylab("Hierarchical Stability")+
                  labs(color="")+
                  scale_colour_manual(values = cbPalette) +
                  theme_bw()+
                  theme(axis.text=element_text(size=10), axis.title=element_text(size=16),
                        plot.margin = unit(c(1,1,1,1), "cm"),
                        legend.key.size = unit(0.8, 'cm'),
                        legend.text = element_text(size=12))

    plotout2 <- ggplot(data=plotdat2, aes(x=Order_q, y=Synchrony, color=type))+
                  geom_line(linewidth=1.2)+
                  ylab(label=expression(paste("Hierarchical Synchrony")))+
                  labs(color="Hierarchical level")+
                  scale_colour_manual(values = cbPalette_2) +
                  theme_bw()+
                  theme(axis.text=element_text(size=10), axis.title=element_text(size=16),
                        plot.margin = unit(c(1,1,1,1), "cm"),
                        legend.key.size = unit(0.8, 'cm'),
                        legend.text = element_text(size=12))

    plotout2_2 <- ggplot(data=plotdat2_2, aes(x=Order_q, y=Stability, color=type))+
                    geom_line(linewidth=1.2)+
                    ylab(label=expression(paste("Hierarchical Stability")))+
                    labs(color="")+
                    scale_colour_manual(values = cbPalette) +
                    theme_bw()+
                    theme(axis.text=element_text(size=10), axis.title=element_text(size=16),
                          plot.margin = unit(c(1,1,1,1), "cm"),
                          legend.key.size = unit(0.8, 'cm'),
                          legend.text = element_text(size=12))

    # if legend in the bottom can use
    # if(length(type_diff)>=4){
    #   plotout2 <- plotout2 + guides(color=guide_legend(nrow=2,byrow=TRUE))
    # }


    # if(sum(output$Order_q==1)!=0 | sum(output$Order_q==2)!=0){
    #   output_q1 <- filter(output, Order_q==1)
    #   output_q1 <- output_q1[order(output_q1$Hier,decreasing=TRUE),]
    #   output_q2 <- filter(output, Order_q==2)
    #   output_q2 <- output_q2[order(output_q2$Hier,decreasing=TRUE),]
    #
    #   label_names <- c(paste("Beta(",(output_q1$Hier-1)[-length(output_q1$Hier)],")", sep=""), "Alpha(1)")
    #   if(nrow(output_q1)==0){
    #     q_stab_q1 <- c()
    #   }else{
    #     q_stab_q1 <- data.frame(Hier=label_names,
    #                             Order_q=rep(1,length(output_q1$Hier)),
    #                             beta_alpha=c((output_q1$Gamma-output_q1$Alpha)[2:length(output_q1$Hier)], output_q1$Alpha[length(output_q1$Hier)]))
    #   }
    #
    #   if(nrow(output_q2)==0){
    #     q_stab_q2 <- c()
    #   }else{
    #     q_stab_q2 <- data.frame(Hier=label_names,
    #                             Order_q=rep(2,length(output_q2$Hier)),
    #                             beta_alpha=c((output_q2$Gamma-output_q2$Alpha)[2:length(output_q2$Hier)], output_q2$Alpha[length(output_q2$Hier)]))
    #   }
    #
    #   q_stab_q1$percent <- ((q_stab_q1$beta_alpha)/sum(q_stab_q1$beta_alpha))*100
    #   q_stab_q2$percent <- ((q_stab_q2$beta_alpha)/sum(q_stab_q2$beta_alpha))*100
    #
    #   plotdat_q12 <- rbind(q_stab_q1, q_stab_q2)
    #   plotdat_q12$Hier <- factor(plotdat_q12$Hier, levels=label_names)
    #   plotdat_q12$Order_q <- as.factor(plotdat_q12$Order_q)
    #
    #   # Check if the number of unique 'Assemblage' is 4 or less
    #   if (length(unique(plotdat_q12$Hier)) <= 4){
    #     cbPalette2 <- c("#EA0000","#0066CC","#73BF00","#FFAF60")
    #   }else{
    #     # If there are more than 4 assemblages, start with the same predefined color palette
    #     # Then extend the palette by generating additional colors using the 'ggplotColors' function
    #     cbPalette2 <- c(c("#EA0000","#0066CC","#73BF00","#FFAF60"),
    #                    ggplotColors(length(unique(plotdat_q12$Hier))-4))
    #   }
    #
    #   plotout3 <- ggplot(data=plotdat_q12, aes(x=Order_q, y=percent, fill=Hier))+
    #     geom_bar(stat="identity", width = 0.5)+
    #     ylab(label="Stability (% of total)")+
    #     labs(fill = "Hierarchical level")+theme_bw()+
    #     scale_fill_manual(values = cbPalette2) +
    #     theme(axis.text=element_text(size=10), axis.title=element_text(size=13))
    #
    #
    #   label_names2 <- paste("Hier_",(output_q1$Hier-1)[-length(output_q1$Hier)], sep="")
    #   if(nrow(output_q1)==0){
    #     q_stab_q1_2 <- c()
    #   }else{
    #     q_stab_q1_2 <- data.frame(Hier=label_names2, Order_q=rep(1,length(output_q1$Hier)-1),
    #                               beta_alpha=c(output_q1$Synchrony[2:length(output_q1$Hier)]))
    #   }
    #
    #   if(nrow(output_q2)==0){
    #     q_stab_q2_2 <- c()
    #   }else{
    #     q_stab_q2_2 <- data.frame(Hier=label_names2, Order_q=rep(2,length(output_q2$Hier)-1),
    #                               beta_alpha=c(output_q2$Synchrony[2:length(output_q2$Hier)]))
    #   }
    #   q_stab_q1_2[,4] <- q_stab_q1_2[,3]/sum(q_stab_q1_2[,3])
    #   q_stab_q2_2[,4] <- q_stab_q2_2[,3]/sum(q_stab_q2_2[,3])
    #
    #   plotdat_q12_2 <- rbind(q_stab_q1_2, q_stab_q2_2)
    #   colnames(plotdat_q12_2)[c(3,4)] <- c("Synchrony","percent")
    #   plotdat_q12_2$Hier <- factor(plotdat_q12_2$Hier, levels=label_names2)
    #   plotdat_q12_2$Order_q <- as.factor(plotdat_q12_2$Order_q)
    #
    #   # Check if the number of unique 'Assemblage' is 4 or less
    #   if (length(unique(plotdat_q12_2$Hier)) <= 4){
    #     cbPalette3 <- c("#EA0000","#0066CC","#73BF00","#FFAF60")
    #   }else{
    #     # If there are more than 4 assemblages, start with the same predefined color palette
    #     # Then extend the palette by generating additional colors using the 'ggplotColors' function
    #     cbPalette3 <- c(c("#EA0000","#0066CC","#73BF00","#FFAF60"),
    #                     ggplotColors(length(unique(plotdat_q12_2$Hier))-4))
    #   }
    #
    #   plotout4 <- ggplot(data=plotdat_q12_2, aes(x=Order_q, y=percent, fill=Hier))+
    #     geom_bar(stat="identity", width = 0.5)+
    #     ylab(label="Synchrony (% of total)")+
    #     labs(fill = "Hierarchical level")+theme_bw()+
    #     scale_fill_manual(values = cbPalette3) +
    #     theme(axis.text=element_text(size=10), axis.title=element_text(size=13))
    # }

    # if(sum(output$Order_q==1)!=0 | sum(output$Order_q==2)!=0){
    #   plotout <- list()
    #   # plotout[[1]] <- ggarrange(plotout1, plotout2, ncol = 2)
    #   plotout[[1]] <- plotout1
    #   plotout[[2]] <- ggarrange(plotout3, plotout4, ncol = 2)
    # }else{
    #   # plotout <- ggarrange(plotout1, plotout2, ncol = 2)
    #   plotout <- list()
    #   plotout[[1]] <- plotout1
    # }

    plotout <- list()
    plotout[[1]] <- plotout1
    plotout[[2]] <- ggarrange(plotout2_2, plotout2, ncol = 2)

  }else{
    output$Site <- factor(output$Site, levels=unique(output$Site))

    # plotout <- list()
    stab_plotdat <- data.frame(Site=rep(output$Site,4),
                               Order_q=rep(output$Order_q,4),
                               value=c(output$Gamma,output$Alpha,output$Beta,output$Synchrony),
                               type=rep(c("Gamma","Alpha","Beta","Synchrony"),each=nrow(output)))
    stab_plotdat$type <- factor(stab_plotdat$type, levels=c("Gamma","Alpha","Beta","Synchrony"))

    # Check if the number of unique 'Assemblage' is 4 or less
    if (length(unique(stab_plotdat$Site)) <= 4){
      cbPalette <- c("#EA0000","#0066CC","#64A600","#FFA042")
    }else{
      # If there are more than 4 assemblages, start with the same predefined color palette
      # Then extend the palette by generating additional colors using the 'ggplotColors' function
      cbPalette <- c(c("#EA0000","#0066CC","#64A600","#FFA042"),
                      ggplotColors(length(unique(stab_plotdat$Site))-4))
    }

    plotout <- ggplot(data=stab_plotdat, aes(x=Order_q, y=value, color=Site))+
      geom_line(linewidth=1.2)+
      facet_wrap(.~type, nrow=2, scales = "free")+
      ylab(label="Stability and Synchrony")+
      xlab(label="Order of q")+
      labs(color="Site")+ theme_bw()+
      scale_color_manual(values = cbPalette) +
      theme(strip.text = element_text(size=13), legend.title = element_text(size=13),
            legend.text = element_text(size=12), legend.key.size = unit(0.8, 'cm'),
            axis.title = element_text(size=16))
    #
    #     plotout[[2]] <- ggplot(data=output, aes(x=Order_q, y=Synchrony, color=Site))+
    #                       geom_line(linewidth=1.2)+
    #                       ylab(label="Synchrony")+
    #                       xlab(label="Order of q")+
    #                       labs(color="Site")+ theme_bw()+
    #                       theme(legend.title = element_text(size=13), legend.text = element_text(size=12),
    #                             legend.key.size = unit(0.8, 'cm'), axis.title = element_text(size=16))
  }

  return(plotout)
}



#' ggplot2 extension for a Stay_Single or Stay_Multiple object to analysis with an diversity (or other) variable.
#'
#' \code{ggStay_analysis} is a graphical function that based on the output from the function \code{Stay_Single} or \code{Stay_Multiple}. It provides to graph relationships between stability (and synchrony if is multiple assemblages) and an additional diversity (or other) variable .
#'
#' @param output the output obtained from \code{Stay_Single} or \code{Stay_Multiple} and needs to combine with a column that sets as \code{x_variable}. Also, if \code{by_group} is not \code{NULL}, the output also need to combine with the column that sets as \code{by_group}.
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
#' @return For an \code{Stay_Single} object, this function return a figure of diversity (or other) variable vs. stability.
#' For an \code{Stay_Multiple} object, this function return a figure that is about diversity (or other) variable vs. (Gamma, Alpha, Beta) stability and synchrony.
#'
#' @examples
#' data("Jena_plot_biomass_data")
#' data("Jena_species_biomass_data")
#' data("Jena_hierarchical_biomass_data")
#' data("Jena_hierarchical_mat")
#'
#' ## Single assemblage
#' # Stability of each single plot
#' single_plot <- do.call(rbind, Jena_plot_biomass_data)
#' output_single_plot_div <- Stay_Single(data=single_plot, order.q=c(1,2), Alltime=TRUE)
#' output_single_plot_div <- data.frame(output_single_plot_div,
#'                                      sowndiv=as.numeric(do.call(rbind,
#'                                        strsplit(output_single_plot_div[,1],"[._]+"))[,2]),
#'                                      block=do.call(rbind,
#'                                        strsplit(output_single_plot_div[,1],"[._]+"))[,1])
#'
#' ggStay_analysis(output=output_single_plot_div, x_variable="sowndiv",
#'                     by_group="block", model="LMM")
#'
#' # Stability of each single species
#' single_species <- do.call(rbind, Jena_species_biomass_data)
#' output_single_species_div <- Stay_Single(data=single_species,
#'                                          order.q=c(1,2), Alltime=TRUE)
#' output_single_species_div <- data.frame(output_single_species_div,
#'                               sowndiv=as.numeric(do.call(rbind,
#'                                       strsplit(output_single_species_div[,1],"[._]+"))[,3]),
#'                               block=do.call(rbind,
#'                                     strsplit(output_single_species_div[,1],"[._]+"))[,2])
#'
#' ggStay_analysis(output=output_single_species_div, x_variable="sowndiv",
#'                     by_group="block", model="LMM")
#'
#'
#' ## Multiple assemblages
#' # Stability of multiple plots
#' multiple_plot <- Jena_plot_biomass_data
#' output_multi_plot_div <- Stay_Multiple(data=multiple_plot, order.q=c(1,2), Alltime=TRUE)
#' output_multi_plot_div <- data.frame(output_multi_plot_div, sowndiv=rep(c(16,8,4,2,1),8),
#'                                     block=rep(rep(c("B1","B2","B3","B4"),each=5),2))
#'
#' ggStay_analysis(output=output_multi_plot_div, x_variable="sowndiv",
#'                     by_group="block", model="LMM")
#'
#' # Stability of multiple species in plot
#' multiple_species <- Jena_species_biomass_data
#' output_multi_species_div <- Stay_Multiple(data=multiple_species,
#'                                               order.q=c(1,2), Alltime=TRUE)
#' output_multi_species_div <- data.frame(output_multi_species_div,
#'                              sowndiv=as.numeric(do.call(rbind,
#'                                      strsplit(output_multi_species_div[,1],"_"))[,3]),
#'                              block=do.call(rbind,
#'                                    strsplit(output_multi_species_div[,1],"_"))[,2])
#'
#' ggStay_analysis(output=output_multi_species_div, x_variable="sowndiv",
#'                     by_group="block", model="LMM")
#'
#' @export

ggStay_analysis <- function(output, x_variable, by_group=NULL, model="LMM"){

  # x_variable_ori <- x_variable
  # if(is.null(by_group)==FALSE){by_group_ori <- by_group}

  if(length(which(colnames(output)=="Stability"))!=0){
    if((length(which(colnames(output)=="Plot/Community"))==0 & length(which(colnames(output)=="Plot.Community"))==0) | length(which(colnames(output)=="Order_q"))==0){
      stop('Please put the complete output of "Stay_Single" or "Stay_Multiple" function.')
    }
  }else{
    if(length(which(colnames(output)=="Site"))==0 | length(which(colnames(output)=="Order_q"))==0 | length(which(colnames(output)=="Gamma"))==0
       | length(which(colnames(output)=="Alpha"))==0 | length(which(colnames(output)=="Beta"))==0　
       | length(which(colnames(output)=="Synchrony"))==0){
      stop('Please put the complete output of "Stay_Single" or "Stay_Multiple" function.')
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
    part_fit <- c()
    plotdata_text_part <- c()
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

        ran <- c()
        for(bb in rownames(coef(MODEL)$Gvariable)){
          ran <- rbind(ran, range(filter(subdata, Gvariable==bb)$Xvariable))
        }
        slope_intercept <- data.frame(Order_q=rep(paste("q = ",qq,sep=""),nrow(coef(MODEL)$Gvariable)),
                                      Gvariable=rownames(coef(MODEL)$Gvariable),
                                      intercept=coef(MODEL)$Gvariable[,1],
                                      slope=coef(MODEL)$Gvariable[,2],
                                      ran)

        colnames(slope_intercept)[5:6] <- c("x_min","x_max")
        part_fit <- rbind(part_fit, slope_intercept)

        plotdata_text_part <- rbind(plotdata_text_part,
                                    data.frame(slope = paste("slope = ", format(round(coef(MODEL)$Gvariable[,2],4), nsmall=4), sep=""),
                                               Order_q = rep(paste("q = ",qq,sep=""), nrow(coef(MODEL)$Gvariable)),
                                               Gvariable=rownames(coef(MODEL)$Gvariable)))
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

    if(model=="LMM"){
      plotdata_text_part$Order_q <- as.factor(plotdata_text_part$Order_q)
      plotdata_text_part$Gvariable <- as.factor(plotdata_text_part$Gvariable)
      part_fit$Order_q <- as.factor(part_fit$Order_q)
      part_fit$Gvariable <- as.factor(part_fit$Gvariable)
    }


    if(is.null(by_group)==FALSE){

      # Check if the number of unique 'Assemblage' is 4 or less
      if (length(unique(plotdata$Gvariable)) <= 4){
        cbPalette <- c("#EA0000","#64A600","#0066CC","#FFA042")
      }else{
        # If there are more than 4 assemblages, start with the same predefined color palette
        # Then extend the palette by generating additional colors using the 'ggplotColors' function
        cbPalette <- c(c("#EA0000","#64A600","#0066CC","#FFA042"),
                       ggplotColors(length(unique(plotdata$Gvariable))-4))
      }

      ### 0717 revise
      if(length(unique(plotdata$Gvariable)) == 4){
        hhjust1 <- 1
        vvjust1 <- -3.5
        hhjust2 <- rep(c(2.1,1),4)
        vvjust2 <- rep(c(-1.7,-1.7,0.1,0.1),2)
      }else{
        gnum <- length(unique(plotdata$Gvariable))
        if(gnum%%2==0){
          create <- gnum/2
        }else{
          create <- (gnum + 1)/2
        }
        hhjust1 <- 1
        vcreate <- 0.1+(-1.8)*c(create:0)
        vvjust1 <- vcreate[1]
        if(gnum%%2==1){
          hhjust2 <- rep(c(2.1,1),create)[-length(rep(c(2.1,1),create))]
          hhjust2 <- rep(hhjust2,2)
          vvjust2 <- rep(vcreate[-1], each=2)[-length(rep(vcreate[-1], each=2))]
          vvjust2 <- rep(vvjust2, 2)
        }else{
          hhjust2 <- rep(rep(c(2.1,1),create),2)
          vvjust2 <- rep(rep(vcreate[-1], each=2),2)
        }
      }

      if(model=="LMM"){
        plotout <- ggplot()+
          geom_point(data = plotdata, aes(x=Xvariable, y=Stability, color=Gvariable), size=2.7)+
          geom_segment(data = part_fit,
                       aes(x=x_min, xend=x_max, y=intercept+slope*x_min, yend=intercept+slope*x_max, color=Gvariable))+
          geom_line(data = plotdata, aes(x=Xvariable, y=pred, linetype = sign), linewidth=1.2, color="black")+
          scale_linetype_manual(values=c("solid","dashed"), drop = FALSE)+
          scale_color_manual(values = cbPalette)+
          facet_wrap(.~Order_q, scales="fixed", ncol = min(5, length(unique(plotdata$Order_q))))+
          labs(linetype="", color=by_group)+
          xlab(label=x_variable)+
          ylab(label="Stability")+ theme_bw()+
          geom_text(data=slope_text, aes(x = -Inf, y = -Inf, label = slope),
                    x=max(plotdata$Xvariable, na.rm = TRUE),
                    y=min(plotdata$Stability, na.rm = TRUE), size=3.5, hjust=hhjust1, vjust=vvjust1)+
          geom_text(data=plotdata_text_part,
                    aes(x = -Inf, y = -Inf, label=slope, color=Gvariable),
                    x=max(plotdata$Xvariable, na.rm = TRUE),
                    y=min(plotdata$Stability, na.rm = TRUE), size=3.5,
                    hjust=hhjust2, vjust=vvjust2, key_glyph=draw_key_path)+
          theme(strip.text = element_text(size=13), legend.title = element_text(size=13),
                legend.text = element_text(size=12), legend.key.size = unit(0.8, 'cm'),
                axis.title = element_text(size=16))+
          guides(linetype = guide_legend(keywidth = 2.5))

      }else{
        plotout <- ggplot()+
          geom_point(data = plotdata, aes(x=Xvariable, y=Stability, color=Gvariable), size=2.7)+
          geom_line(data = plotdata, aes(x=Xvariable, y=pred, linetype = sign), linewidth=1.2, color="black")+
          scale_linetype_manual(values=c("solid","dashed"), drop = FALSE)+
          scale_color_manual(values = cbPalette)+
          facet_wrap(.~Order_q, scales="fixed", ncol = min(5, length(unique(plotdata$Order_q))))+
          labs(linetype="", color=by_group)+
          xlab(label=x_variable)+
          ylab(label="Stability")+ theme_bw()+
          geom_text(data=slope_text, aes(x = -Inf, y = -Inf, label = slope),
                    x=max(plotdata$Xvariable, na.rm = TRUE),
                    y=min(plotdata$Stability, na.rm = TRUE), size=5, hjust=1, vjust=0.1)+
          theme(strip.text = element_text(size=13), legend.title = element_text(size=13),
                legend.text = element_text(size=12), legend.key.size = unit(0.8, 'cm'),
                axis.title = element_text(size=16))+
          guides(linetype = guide_legend(keywidth = 2.5))
      }

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
                          axis.title = element_text(size=16))+
                    guides(linetype = guide_legend(keywidth = 2.5))
    }

  }else{

    # fit model
    lm_sign <-c()
    lm_slope <- c()
    plotdata <- c()
    plotdata_text_part <- c()
    part_fit <- c()

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

        ran <- c()
        for(bb in rownames(coef(MODEL_G)$Gvariable)){
          ran <- rbind(ran, range(filter(subdata, Gvariable==bb)$Xvariable))
        }
        slope_intercept <- data.frame(Order_q=rep(paste("q = ",qq,sep=""),nrow(coef(MODEL_G)$Gvariable)*4),
                                      type=rep(c("Gamma","Alpha","Beta","Synchrony"), each=nrow(coef(MODEL_G)$Gvariable)),
                                      Gvariable=rep(rownames(coef(MODEL_G)$Gvariable), 4),
                                      intercept=c(coef(MODEL_G)$Gvariable[,1],
                                                  coef(MODEL_A)$Gvariable[,1],
                                                  coef(MODEL_BA)$Gvariable[,1],
                                                  coef(MODEL_S)$Gvariable[,1]),
                                      slope=c(coef(MODEL_G)$Gvariable[,2],
                                              coef(MODEL_A)$Gvariable[,2],
                                              coef(MODEL_BA)$Gvariable[,2],
                                              coef(MODEL_S)$Gvariable[,2]),
                                      rbind(ran,ran,ran,ran))
        colnames(slope_intercept)[6:7] <- c("x_min","x_max")
        part_fit <- rbind(part_fit, slope_intercept)

        plotdata_text_part <- rbind(plotdata_text_part,
                                    data.frame(slope_gamma=paste("slope = ", format(round(coef(MODEL_G)$Gvariable[,2],4), nsmall=4), sep=""),
                                               slope_alpha=paste("slope = ", format(round(coef(MODEL_A)$Gvariable[,2],4), nsmall=4), sep=""),
                                               slope_beta=paste("slope = ", format(round(coef(MODEL_BA)$Gvariable[,2],4), nsmall=4), sep=""),
                                               slope_synchrony=paste("slope = ", format(round(coef(MODEL_S)$Gvariable[,2],4), nsmall=4), sep=""),
                                               Order_q=rep(paste("q = ",qq,sep=""), nrow(coef(MODEL_G)$Gvariable)),
                                               Gvariable=rownames(coef(MODEL_G)$Gvariable)))
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

    plotdata_Stab <- data.frame(Site = rep(plotdata$Site,4),
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


    slope_text_Stab <- data.frame(slope = paste("slope = ",format(round(as.vector(as.matrix(lm_slope)),4), nsmall=4),sep=""),
                                  Order_q = rep(rownames(lm_slope),4),
                                  type = rep(c("Gamma","Alpha","Beta","Synchrony"), each = nrow(lm_slope)))
    slope_text_Stab$type <- factor(slope_text_Stab$type, levels = c("Gamma","Alpha","Beta","Synchrony"))
    slope_text_Stab$Order_q <- as.factor(slope_text_Stab$Order_q)

    if(model=="LMM"){
      plotdata_text_part <- data.frame(slope = c(plotdata_text_part$slope_gamma,
                                                 plotdata_text_part$slope_alpha,
                                                 plotdata_text_part$slope_beta,
                                                 plotdata_text_part$slope_synchrony),
                                       Order_q = rep(plotdata_text_part$Order_q, 4),
                                       type = rep(c("Gamma","Alpha","Beta","Synchrony"), each=nrow(plotdata_text_part)),
                                       Gvariable = rep(plotdata_text_part$Gvariable, 4))
      plotdata_text_part$type <- as.factor(plotdata_text_part$type)
      plotdata_text_part$Order_q <- as.factor(plotdata_text_part$Order_q)
      plotdata_text_part$Gvariable <- as.factor(plotdata_text_part$Gvariable)
      part_fit$Gvariable <- as.factor(part_fit$Gvariable)
      part_fit$Order_q <- as.factor(part_fit$Order_q)
      part_fit$type <- as.factor(part_fit$type)
    }

    # slope_text_Syn <- data.frame(slope = paste("slope = ",round(lm_slope[,5],4),sep=""),
    #                              Order_q = rownames(lm_slope))
    # slope_text_Syn$Order_q <- as.factor(slope_text_Syn$Order_q)

    tyG <- min(filter(plotdata_Stab, type=="Gamma")$Stability, na.rm = TRUE)
    tyA <- min(filter(plotdata_Stab, type=="Alpha")$Stability, na.rm = TRUE)
    # tyBM <- max(filter(plotdata_Stab, type=="Beta (multiplicative)")$Stability)
    tyBA <- max(filter(plotdata_Stab, type=="Beta")$Stability, na.rm = TRUE)
    tyS <- min(filter(plotdata_Stab, type=="Synchrony")$Stability, na.rm = TRUE)

    if(is.null(by_group)==FALSE){

      # Check if the number of unique 'Assemblage' is 4 or less
      if (length(unique(plotdata_Stab$Gvariable)) <= 4){
        cbPalette <- c("#EA0000","#64A600","#0066CC","#FFA042")
      }else{
        # If there are more than 4 assemblages, start with the same predefined color palette
        # Then extend the palette by generating additional colors using the 'ggplotColors' function
        cbPalette <- c(c("#EA0000","#64A600","#0066CC","#FFA042"),
                       ggplotColors(length(unique(plotdata_Stab$Gvariable))-4))
      }

      if(length(unique(plotdata$Gvariable)) == 4){
        gnum <- 4
        hhjust1 <- 1
        vvjust1 <- rep(c(-3.5,-3.5,4.1,-3.5),each=2)
        hhjust2 <- rep(c(2.1,1),16)
        vvjust2 <- c(rep(c(-1.7,-1.7,0.1,0.1),4), rep(c(1.1,1.1,2.6,2.6),2), rep(c(-1.7,-1.7,0.1,0.1),2))
      }else{
        gnum <- length(unique(plotdata$Gvariable))
        if(gnum%%2==0){
          create <- gnum/2
        }else{
          create <- (gnum + 1)/2
        }
        hhjust1 <- 1
        vcreate_1 <- 0.1+(-1.8)*c(create:0)
        vcreate_2 <- 1.1+1.5*c(create:0)
        vvjust1 <- rep(c(vcreate_1[1], vcreate_1[1], vcreate_2[1], vcreate_1[1]), each=2)
        if(gnum%%2==1){
          hhjust2 <- rep(c(2.1,1),create)[-length(rep(c(2.1,1),create))]
          hhjust2 <- rep(hhjust2,8)
          vvjust2_1 <- rep(vcreate_1[-1], each=2)[-length(rep(vcreate_1[-1], each=2))]
          vvjust2_2 <- rev(rep(vcreate_2[-1], each=2))[-length(rev(rep(vcreate_2[-1], each=2)))]
          vvjust2 <- c(rep(vvjust2_1, 4), rep(vvjust2_2, 2), rep(vvjust2_1, 2))
        }else{
          hhjust2 <- rep(rep(c(2.1,1),create),8)
          vvjust2 <- c(rep(rep(vcreate_1[-1], each=2),4), rep(rep(rev(vcreate_2[-1]), each=2),2),
                       rep(rep(vcreate_1[-1], each=2),2))
        }
      }


      if(model=="LMM"){
        plotout1 <- ggplot()+
          geom_point(data = plotdata_Stab, aes(x=Xvariable, y=Stability, color=Gvariable), size=2.7)+
          geom_segment(data=part_fit,
                       aes(x=x_min, xend=x_max, y=intercept+slope*x_min, yend=intercept+slope*x_max, color=Gvariable))+
          geom_line(data = plotdata_Stab, aes(x=Xvariable, y=pred, linetype = sign), linewidth=1.2, color="black")+
          scale_linetype_manual(values=c("solid","dashed"), drop = FALSE)+
          scale_color_manual(values = cbPalette)+
          facet_grid(type~Order_q, scales = "free_y")+
          labs(linetype="", color=by_group)+
          xlab(label=x_variable)+
          ylab(label="Stability and Synchrony")+ theme_bw()+
          geom_text(data=slope_text_Stab, aes(x = -Inf, y = -Inf, label = slope),
                    x=max(plotdata_Stab$Xvariable, na.rm = TRUE),
                    y=rep(c(tyG, tyA, tyBA, tyS),each=2), size=3.5,
                    hjust=hhjust1, vjust=vvjust1)+
          geom_text(data=plotdata_text_part,
                    aes(x = -Inf, y = -Inf, label=slope, color=Gvariable),
                    x=rep(max(plotdata_Stab$Xvariable, na.rm = TRUE),gnum*8),
                    y=rep(c(tyG, tyA, tyBA, tyS),each=gnum*2), size=3.5,
                    hjust=hhjust2, vjust=vvjust2, key_glyph=draw_key_path)+
          theme(strip.text = element_text(size=13), legend.title = element_text(size=13),
                legend.text = element_text(size=12), legend.key.size = unit(0.8, 'cm'),
                axis.title = element_text(size=16))+
          guides(linetype = guide_legend(keywidth = 2.5))

      }else{
        plotout1 <- ggplot()+
          geom_point(data = plotdata_Stab, aes(x=Xvariable, y=Stability, color=Gvariable), size=2.7)+
          geom_line(data = plotdata_Stab, aes(x=Xvariable, y=pred, linetype = sign), linewidth=1.2, color="black")+
          scale_linetype_manual(values=c("solid","dashed"), drop = FALSE)+
          scale_color_manual(values = cbPalette)+
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
                axis.title = element_text(size=16))+
          guides(linetype = guide_legend(keywidth = 2.5))
      }



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
                                hjust=1, vjust=rep(c(0.1,0.1,1.1,0.1),each=2))+
                      theme(strip.text = element_text(size=13), legend.title = element_text(size=13),
                            legend.text = element_text(size=12), legend.key.size = unit(0.8, 'cm'),
                            axis.title = element_text(size=16))+
                      guides(linetype = guide_legend(keywidth = 2.5))

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



# Generate Color Palette for ggplot2
#
# This function creates a color palette suitable for ggplot2 visualizations by evenly spacing colors in the HCL color space. The function ensures that the colors are well-distributed and visually distinct, making it ideal for categorical data where each category needs to be represented by a different color.
#
# @param g An integer indicating the number of distinct colors to generate. This value should be a positive integer, with higher values resulting in a broader range of colors.
# @return A vector of color codes in hexadecimal format, suitable for use in ggplot2 charts and plots. The length of the vector will match the input parameter `g`.
# @examples
# # Generate a palette of 5 distinct colors
# ggplotColors(5)
#
# # Use the generated colors in a ggplot2 chart
# library(ggplot2)
# df <- data.frame(x = 1:5, y = rnorm(5), group = factor(1:5))
# ggplot(df, aes(x, y, color = group)) +
#   geom_point() +
#   scale_color_manual(values = ggplotColors(5))
#
ggplotColors <- function(g){
  d <- 360/g # Calculate the distance between colors in HCL color space
  h <- cumsum(c(15, rep(d,g - 1))) # Create cumulative sums to define hue values
  hcl(h = h, c = 100, l = 65) # Convert HCL values to hexadecimal color codes
}





## ========== no visible global function definition for R CMD check ========== ##
utils::globalVariables(c("COLOR", "Order_q", "Site", "Stability", "Synchrony"
))




