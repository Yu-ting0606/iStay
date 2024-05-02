#' Calculate stability of the time series data for single assemblage.
#'
#' \code{Stab_Single} is a function that calculate stability of the time series data (like biomass, productivity, etc.) for single assemblage.
#'
#' @param data can be input as a \code{vector} of time series data, or \code{data.frame/matrix} (assemblages by times). The rownames can be added with "_" and property, that the property will be used as a factor to plot different colors in \code{Stab_Syn_ggplot} function.
#' @param order.q a numerical vector specifying the orders of stability. Default is c(1,2).
#' @param Alltime \code{TRUE} or \code{FALSE}, to decide whether to use all the times in the data.
#' @param start_T (argument only for \code{Alltime = FALSE}) a positive integer specifying the start column of time in the data.
#' @param end_T (argument only for \code{Alltime = FALSE}) a positive integer specifying the end column of time in the data.
#'
#'
#'
#' @return a dataframe with columns "Assemblage", "Order_q" and "Stability".
#'
#' @examples
#' data("Jena_experiment_plant_data")
#' output <- Stab_Single(data=Jena_experiment_plant_data[[1]],
#'                       order.q=seq(0.2,2,0.2), Alltime=TRUE)
#' output
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
  if(sum(which(apply(data,1,sum)==0))){
    stop('There is at least one assemblage whose time series data are all zeros.')
  }

  if(Alltime!=TRUE){
    if(is.null(start_T) | is.null(end_T)){
       stop('Need to set the length of time series for calculating.')
    }
    if((start_T>=end_T) | length(start_T)!=1 | length(end_T)!=1){
       stop('Starting and ending time need to be a number, and ending time needs larger than starting time.')
    }
  }

  stability <- function(vector,q){
    K <- length(vector)
    if(q==1){
      vector <- vector[vector!=0]
      H <- sum((vector/sum(vector))*log(vector/sum(vector)))*(-1)
      out <- H/log(K)
    }else{
      up <- 1-sum((vector/sum(vector))^q)
      out <- up/(1-K^(1-q))
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
  return(result)
}





#' Calculate stability and synchrony of the time series data for multiple assemblages.
#'
#' \code{Stab_Syn_Multiple} is a function that calculate (Gamma, Alpha and Beta) stability and synchrony of the time series data (like biomass, productivity, etc.) for multiple assemblages.
#'
#' @param data can be input as a \code{data.frame/matrix} (assemblages by times), or a \code{list} of \code{data.frames/matrices} with each dataframe representing a assemblages-by-times data. The names of dataframes can be added with "_" and property, that the property will be used as a factor to plot different colors in \code{Stab_Syn_ggplot} function.
#' @param order.q a numerical vector specifying the orders of stability and synchrony. Default is c(1,2).
#' @param Alltime \code{TRUE} or \code{FALSE}, to decide whether to use all the times in (every) dataframe.
#' @param start_T (argument only for \code{Alltime = FALSE}) a positive integer specifying the start column of time in (every) dataframe.
#' @param end_T (argument only for \code{Alltime = FALSE}) a positive integer specifying the end column of time in (every) dataframe.
#'
#'
#' @return a dataframe with columns "Assemblage", "Order_q", "Gamma Stability", "Alpha Stability", "Beta Stability (use multiple decomposition)", "Beta Stability (use additional decomposition)" and "Synchrony".
#'
#' @examples
#' data("Jena_experiment_plant_data")
#' output <- Stab_Syn_Multiple(Jena_experiment_plant_data,
#'                             order.q=seq(0.2,2,0.2), Alltime=TRUE)
#' output
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

  if(is.data.frame(data) | is.matrix(data)){
    if(nrow(data)<=1){
      stop('There is at least one list that has only a single assemblage.')
    }
    if(sum(which(apply(data,1,sum)==0))!=0){
      stop('There is at least one assemblage whose time series data are all zeros.')
    }
  }else{
    num_row <- do.call(rbind, lapply(data, function(ZZ) nrow(ZZ)))
    if(sum(which(as.vector(num_row)<=1))!=0){
      stop('There is at least one list that has only a single assemblage.')
    }
    zero <- do.call(rbind, lapply(data, function(ZZ) sum(which(apply(ZZ,1,sum)==0))))
    if(sum(which(as.vector(zero)!=0))!=0){
      stop('There is at least one assemblage whose time series data are all zeros.')
    }
  }

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
    ZZ[is.na(ZZ)] <- 0

    z_iplus <- apply(ZZ,1,sum)
    if(length(which(z_iplus==0))!=0){
      ZZ <- ZZ[-which(z_iplus==0),]
      z_iplus <- z_iplus[z_iplus!=0]
    }

    z_plusk <- apply(ZZ,2,sum)
    z_plusplus <- sum(ZZ)
    p_i <- apply(ZZ,2,function(w)w/z_iplus)
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
    return(c(gamma, alpha, gamma/alpha, gamma-alpha))
  }

  #type="C"or"U"
  Synchrony <- function(ZZ,q,type="C"){
    M <- nrow(ZZ)
    K <- ncol(ZZ)
    ZZ[is.na(ZZ)] <- 0
    z_iplus <- apply(ZZ,1,sum)
    if(length(which(z_iplus==0))!=0){
      ZZ <- ZZ[-which(z_iplus==0),]
      z_iplus <- z_iplus[z_iplus!=0]
    }

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
    colnames(result)[3:7] <- c("Stab_Gamma", "Stab_Alpha", "Stab_Beta_multiple", "Stab_Beta_additional", "Synchrony")

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
                            colnames(result)[2:6] <- c("Stab_Gamma", "Stab_Alpha", "Stab_Beta_multiple", "Stab_Beta_additional", "Synchrony")
                            return(result)
                          })
                  cal2 <- do.call(rbind, cal)
                  calcal <- data.frame(Place=names(data), cal2)
                  return(calcal)
           })
    result <- do.call(rbind, out)
  }
  return(result)

}





#' Calculate stability and synchrony of the time series data for multiple assemblages.
#'
#' \code{Stab_Syn_ggplot} is a function that calculate (Gamma, Alpha and Beta) stability and synchrony of the time series data (like biomass, productivity, etc.) for multiple assemblages.
#'
#' @param output the output obtained from \code{Stab_Single} or \code{Stab_Syn_Multiple}.
#'
#'
#' @import ggplot2
#' @import stringr
#'
#'
#' @return For an \code{Stab_Single} object, this function return a figure of order of q vs. stability.
#' For an \code{Stab_Syn_Multiple} object, this function return two figures that are order of q vs. (Gamma, Alpha, Beta) stability and order of q vs. synchrony.
#'
#' @examples
#' data("Jena_experiment_plant_data")
#'
#' output_single <- Stab_Single(data=Jena_experiment_plant_data[[1]],
#'                              order.q=seq(0.2,2,0.2), Alltime=TRUE)
#' Stab_Syn_ggplot(output_single)
#'
#' output_multi <- Stab_Syn_Multiple(Jena_experiment_plant_data,
#'                                   order.q=seq(0.2,2,0.2), Alltime=TRUE)
#' Stab_Syn_ggplot(output_multi)
#'
#' @export

Stab_Syn_ggplot <- function(output){

  if((ncol(output)!=3) & (ncol(output)!=7)){
    stop('Please put the complete output of "Stab_Single" or "Stab_Syn_Multiple" function.')
  }

  if(length(which(colnames(output)=="Stability"))!=0){

    if(sum(colnames(output)==c("Assemblage","Order_q","Stability"))!=3){
      stop('Colnames of output need to be "Assemblage","Order_q" and "Stability".')
    }

    split <- str_split(output$Assemblage, pattern = "_", n = Inf, simplify = TRUE)
    coloruse <- split[,ncol(split)]

    output <- data.frame(output, COLOR=coloruse)
    output$COLOR <- as.factor(output$COLOR)
    plotout <- ggplot(output, aes(x=Order_q, y=Stability, color=COLOR))+
                  geom_point(size=2)+
                  geom_line(linewidth=1.2)+
                  theme_bw()+
                  scale_colour_hue(l = 80, c = 170)+
                  ggtitle("Order.q vs. Stability")+
                  theme(plot.title = element_text(hjust = 0.5, size=15), legend.title=element_blank())

  }else{

    if(sum(colnames(output)==c("Place","Order_q","Stab_Gamma", "Stab_Alpha", "Stab_Beta_multiple", "Stab_Beta_additional", "Synchrony"))!=7){
      stop('Colnames of output need to be "Place","Order_q","Stab_Gamma", "Stab_Alpha", "Stab_Beta_multiple", "Stab_Beta_additional" and "Synchrony".')
    }

    split <- str_split(output$Place, pattern = "_", n = Inf, simplify = TRUE)
    coloruse <- split[,ncol(split)]

    output <- data.frame(output, COLOR=coloruse)
    output$COLOR <- factor(output$COLOR, levels=unique(coloruse))

    plotdata <- data.frame(rbind(output[,1:2],output[,1:2],output[,1:2],output[,1:2]),
                           Stability=c(output[,3],output[,4],output[,5],output[,6]),
                           Stab_type=rep(c("Gamma","Alpha","Beta(multiple)","Beta(additional)"), each=nrow(output)),
                           COLOR=rep(output$COLOR, 4))
    plotdata$Stab_type <- factor(plotdata$Stab_type, levels=c("Gamma","Alpha","Beta(multiple)","Beta(additional)"))

    plotout1 <- ggplot(plotdata, aes(x=Order_q, y=Stability, color=COLOR, group=Place))+
                  geom_point(size=2)+
                  geom_line(linewidth=1.2)+
                  theme_bw()+
                  facet_wrap(.~Stab_type, ncol=2, scales = "free")+
                  scale_colour_hue(l = 80, c = 170)+
                  ggtitle("Order.q vs. (Gamma, Alpha and Beta) Stability")+
                  theme(plot.title = element_text(hjust = 0.5, size=15), legend.title=element_blank())

    plotout2 <- ggplot(output, aes(x=Order_q, y=Synchrony, color=COLOR, group=Place))+
                  geom_point(size=2)+
                  geom_line(linewidth=1.2)+
                  theme_bw()+
                  scale_colour_hue(l = 80, c = 170)+
                  ggtitle("Order.q vs. Synchrony")+
                  theme(plot.title = element_text(hjust = 0.5, size=15), legend.title=element_blank())

    plotout <- list()
    plotout[[1]] <- plotout1
    plotout[[2]] <- plotout2
  }

  return(plotout)
}






## ========== no visible global function definition for R CMD check ========== ##
utils::globalVariables(c("COLOR", "Order_q", "Place", "Stability", "Synchrony"
))




