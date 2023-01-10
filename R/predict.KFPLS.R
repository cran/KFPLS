#' @title Prediction by KFPLS
#'
#' @description Prediction of the scalar response by KFPLS.
#'
#' @param object A KFPLS object obtained from \code{\link[KFPLS]{KFPLS}}.
#' @param newdata An \code{array} with three indices denoting the new observations of the functional predictors. The (i, j, k)-th element of it corresponds to the measurment of the i-th subject for the k-th functional predictor at j-th observation grid.
#' @param ... Not used.
#'
#' @return A \code{vector} denoting the prediction of the scalar response.
#' @export
predict.KFPLS <- function(object, newdata, ...){

  n <- object$n
  p <- object$p
  T <- object$T
  U <- object$U
  K <- object$K
  K_c <- object$K_c
  Xfd_list <- object$Xfd_list
  XX_list <- object$XX_list
  Y_c <- object$Y_c
  meanY <- object$meanY
  obser_time <- object$obser_time
  basis <- object$basis
  sigm <- object$sigm

  n_test <- dim(newdata)[1]

  X_test_fd_list <- list()
  XXt_list <- list()
  XtXt_list <- list()
  for(k in 1:p){
    X_test_fd_list[[k]] <- fda::smooth.basis(obser_time, t(newdata[, , k]), basis)$fd
    XXt_list[[k]] <- fda::inprod(X_test_fd_list[[k]], Xfd_list[[k]])
    XtXt_list[[k]] <- fda::inprod(X_test_fd_list[[k]], X_test_fd_list[[k]])
  }

  # RBF
  K_t <- matrix(0, nrow = n_test, ncol = n)
  for(i in 1:n_test){
    for(j in 1:n){
      XX_total <- 0
      for(k in 1:p){
        XX_total <- XX_total + XtXt_list[[k]][i, i] + XX_list[[k]][j, j] - 2 * XXt_list[[k]][i, j]
      }
      K_t[i, j] <- exp(-sigm * XX_total)
    }
  }
  K_t_c <- (K_t - rep(1, n_test) %*% t(rep(1, n)) %*% K /n) %*% (diag(n) - rep(1, n) %*% t(rep(1, n)) /n)
  Y_pre_test <- K_t_c %*% U %*% solve(t(T) %*% K_c %*% U) %*% t(T) %*% Y_c + meanY

  return(Y_pre_test)

}
