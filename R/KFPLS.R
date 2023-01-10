#' @title Kernel functional partial least squares method
#'
#' @description Kernel functional partial least squares (KFPLS) method for functional nonlinear models with scalar response and functional predictors. The Gaussian kernel is used.
#'
#' @param X An \code{array} with three indices. The (i, j, k)-th element of it corresponds to the measurment of the i-th subject for the k-th functional predictor at j-th observation grid.
#' @param Y A \code{vector} with length n, where n is the sample size. The i-th element of it corresponds to the measurement of the scalar response for the i-th subject.
#' @param obser_time A \code{vector} denoting the observation times of the functional predictors.
#' @param nfold An \code{integer} denoting the number of folds for the selection of the tuning parameters by cross-validation.
#' @param n_comp A \code{vector} denoting the candidates of the number of components.
#' @param sigm_list A \code{vector} denoting the candidates of the tuning parameter for the Gaussian kernel.
#' @param basis A basis object denoting the basis that used for the smoothing of the functional predictors. It is created by functions in \code{fda} package, such as \code{\link[fda]{create.bspline.basis}}.
#'
#' @return A \code{list} containing the following components:
#' \item{n}{A \code{scalar} denoting the sample size.}
#' \item{p}{A \code{scalar} denoting the number of functional predictors.}
#' \item{nk}{A \code{scalar} denoting the selected number of components.}
#' \item{T}{A \code{matrix} denoting the value of T at convergence.}
#' \item{U}{A \code{matrix} denoting the value of U at convergence.}
#' \item{K}{A \code{matrix} denoting the Gram matrix.}
#' \item{K_c}{A \code{matrix} denoting the centralized Gram matrix.}
#' \item{Xfd_list}{A \code{list} of length \code{p}. The k-th entry corresponds to the functional data object of the k-th functional predictor.}
#' \item{XX_list}{A \code{list} of length \code{p}. The k-th entry corresponds to the matrix that denotes the inner product of the k-th functional predictor for all subjects.}
#' \item{Y_c}{A \code{vector} denoting the centralized scalar response.}
#' \item{meanY}{A \code{scalar} denoting the sample mean of the scalar response.}
#' \item{Y_hat}{A \code{vector} denoting the prediction of the scalar response.}
#' \item{obser_time}{A \code{vector} denoting the observation times of the functional predictors.}
#' \item{basis}{A basis object denoting the basis that used for the smoothing of the functional predictors.}
#' \item{sigm}{A \code{scalar} denoting the selected tuning parameter for the Gaussian kernel.}
#' \item{CVscore}{A \code{matrix} denoting the CV scores.}
#' \item{time}{A \code{scalar} denoting the computation time.}
#' @export
#'
#' @examples
#' # Generate data
#' n <- 200
#' t_range <- c(0, 1)
#' obser_time <- seq(0, 1, length.out = 51)
#' beta_fun <- function(t){2 * sin(2 * pi * t)}
#' basis <- fda::create.bspline.basis(t_range, nbasis = 13, norder = 4,
#' breaks = seq(0, 1, length.out = 11))
#' beta_fd <- fda::smooth.basis(obser_time, beta_fun(obser_time), basis)$fd
#' X_basis <- fda::create.bspline.basis(t_range, nbasis = 23, norder = 4,
#' breaks = seq(0, 1, length.out = 21))
#' Bbeta <- fda::inprod(X_basis, beta_fd)
#' Xi_B <- splines::bs(obser_time, knots = seq(0, 1, length.out = 21)[-c(1, 21)],
#' degree = 3, intercept = TRUE)
#' a <- array(0, dim = c(n, 23, 1))
#' X <- array(0, dim = c(n, 51, 1))
#' Y <- NULL
#' for(i in 1:n){
#' a[i, , 1] <- stats::rnorm(23)
#' X[i, , 1] <- Xi_B %*% a[i, , 1]
#' aBbeta <- as.numeric(t(a[i, , 1]) %*% Bbeta)
#' Y[i] <- aBbeta + stats::rnorm(1, mean = 0, sd = 0.05)
#' }
#' # KFPLS
#' KFPLS_list <- KFPLS(X, Y, obser_time, nfold = 5, n_comp = 5, sigm_list = 0.005, basis)
#' plot(KFPLS_list$Y_hat, Y)
#' lines(Y, Y)
KFPLS <- function(X, Y, obser_time, nfold, n_comp, sigm_list, basis){

  start_time <- Sys.time()

  n <- length(Y)
  p <- dim(X)[3]

  sample_id <- rep(nfold, n)
  nid_cand <- 1:n
  n_cv <- round(n/nfold)
  for(i in 1:(nfold - 1)){

    id_select <- sample(nid_cand, size = n_cv, replace = F)
    sample_id[id_select] <- i

    nid_cand <- nid_cand[-which(nid_cand %in% id_select)]

  }

  #comoute CV score
  CV_score <- matrix(NA, nrow = length(n_comp), ncol = length(sigm_list))
  for(ncomp_id in 1:(length(n_comp))){

    nk <- n_comp[ncomp_id]

    for(sigm_id in 1:(length(sigm_list))){

      sigm <- sigm_list[sigm_id]

      CV_score_ind <- NULL

      for(id in 1:nfold){

        ##########dataset#############
        if(p > 1){

          X_train <- X[which(sample_id != id), ,]
          Y_train <- Y[which(sample_id != id)]
          n_train <- length(which(sample_id != id))

          X_test <- X[which(sample_id == id), ,]
          Y_test <- Y[which(sample_id == id)]
          n_test <- length(which(sample_id == id))

        }else{

          n_train <- length(which(sample_id != id))
          X_train <- array(X[which(sample_id != id), ,], c(n_train, length(obser_time), p))
          Y_train <- Y[which(sample_id != id)]

          n_test <- length(which(sample_id == id))
          X_test <- array(X[which(sample_id == id), ,], c(n_test, length(obser_time), p))
          Y_test <- Y[which(sample_id == id)]

        }

        ############Computation###########
        Xfd_list <- list()
        XX_list <- list()
        for(j in 1:p){
          Xfd_list[[j]] <- fda::smooth.basis(obser_time, t(X_train[, , j]), basis)$fd
          XX_list[[j]] <- fda::inprod(Xfd_list[[j]], Xfd_list[[j]])
        }

        # RBF
        K <- matrix(0, ncol = n_train, nrow = n_train)
        for(i in 1:n_train){
          for(j in 1:n_train){
            XX_total <- 0
            for(k in 1:p){
              XX_total <- XX_total + XX_list[[k]][i, i] + XX_list[[k]][j, j] - 2 * XX_list[[k]][i, j]
            }
            kij <- exp(-sigm * XX_total)
            K[i, j] <- kij
            K[j, i] <- kij
          }
        }

        K_c <- (diag(n_train) - rep(1, n_train) %*%
                  t(rep(1, n_train))/n_train) %*% K %*%
          (diag(n_train) - rep(1, n_train) %*% t(rep(1, n_train))/n_train)
        meanY <- mean(Y_train)
        Y_c <- Y_train - meanY

        K_run <- K_c
        Y_run <- Y_c
        T <- matrix(0, ncol = nk, nrow = n_train)
        U <- matrix(0, ncol = nk, nrow = n_train)
        for(comp_id in 1:nk){

          u <- stats::rnorm(n_train)
          t <- stats::rnorm(n_train)
          t <- t/sqrt(sum(t^2))
          diff_t <- 1
          m <- 1
          while(diff_t > 10^(-5)){

            t_old <- t

            t <- K_run %*% u
            t <- t/sqrt(sum(t^2))
            c <- t(Y_run) %*% t
            u <- Y_run %*% c
            u <- u/sqrt(sum(u^2))

            diff_t <- norm(t - t_old, type = "F")

            m <- m + 1

          }

          T[, comp_id] <- t
          U[, comp_id] <- u

          #deflation
          tt <- t %*% t(t)
          K_run <- K_run - tt %*% K_run - K_run %*% tt + tt %*% K_run %*% tt
          Y_run <- Y_run - tt %*% Y_run

        }

        Y_test_c <- Y_test - meanY
        X_test_fd_list <- list()
        XXt_list <- list()
        XtXt_list <- list()
        for(k in 1:p){
          X_test_fd_list[[k]] <- fda::smooth.basis(obser_time, t(X_test[, , k]), basis)$fd
          XXt_list[[k]] <- fda::inprod(X_test_fd_list[[k]], Xfd_list[[k]])
          XtXt_list[[k]] <- fda::inprod(X_test_fd_list[[k]], X_test_fd_list[[k]])
        }

        # RBF
        K_t <- matrix(0, nrow = n_test, ncol = n_train)
        for(i in 1:n_test){
          for(j in 1:n_train){
            XX_total <- 0
            for(k in 1:p){
              XX_total <- XX_total + XtXt_list[[k]][i, i] + XX_list[[k]][j, j] - 2 * XXt_list[[k]][i, j]
            }
            K_t[i, j] <- exp(-sigm * XX_total)
          }
        }
        K_t_c <- (K_t - rep(1, n_test) %*% t(rep(1, n_train)) %*% K /n_train) %*%
          (diag(n_train) - rep(1, n_train) %*% t(rep(1, n_train)) /n_train)
        Y_pre_test <- K_t_c %*% U %*% solve(t(T) %*% K_c %*% U) %*% t(T) %*% Y_c

        CV_score_ind[id] <- mean((Y_pre_test - Y_test_c)^2)

      }
      CV_score[ncomp_id, sigm_id] <- mean(CV_score_ind)


    }


  }

  min_id <- which(CV_score == min(CV_score))
  sigm_id <- ceiling(min_id/length(n_comp))
  nk_id <- min_id - (sigm_id - 1) * length(n_comp)
  nk <- n_comp[nk_id]
  sigm <- sigm_list[sigm_id]


  #############Computation for full data############
  Xfd_list <- list()
  XX_list <- list()
  for(j in 1:p){
    Xfd_list[[j]] <- fda::smooth.basis(obser_time, t(X[, , j]), basis)$fd
    XX_list[[j]] <- fda::inprod(Xfd_list[[j]], Xfd_list[[j]])
  }

  # RBF
  K <- matrix(0, ncol = n, nrow = n)
  for(i in 1:n){
    for(j in 1:n){
      XX_total <- 0
      for(k in 1:p){
        XX_total <- XX_total + XX_list[[k]][i, i] + XX_list[[k]][j, j] - 2 * XX_list[[k]][i, j]
      }
      kij <- exp(-sigm * XX_total)
      K[i, j] <- kij
      K[j, i] <- kij
    }
  }


  K_c <- (diag(n) - rep(1, n) %*% t(rep(1, n))/n) %*% K %*% (diag(n) - rep(1, n) %*% t(rep(1, n))/n)
  meanY <- mean(Y)
  Y_c <- Y - meanY

  K_run <- K_c
  Y_run <- Y_c
  T <- matrix(0, ncol = nk, nrow = n)
  U <- matrix(0, ncol = nk, nrow = n)
  for(comp_id in 1:nk){

    u <- stats::rnorm(n)
    t <- stats::rnorm(n)
    t <- t/sqrt(sum(t^2))
    diff_t <- 1
    m <- 1
    while(diff_t > 10^(-5)){

      t_old <- t

      t <- K_run %*% u
      t <- t/sqrt(sum(t^2))
      c <- t(Y_run) %*% t
      u <- Y_run %*% c
      u <- u/sqrt(sum(u^2))

      diff_t <- norm(t - t_old, type = "F")

      m <- m + 1

    }

    T[, comp_id] <- t
    U[, comp_id] <- u

    #deflation
    tt <- t %*% t(t)
    K_run <- K_run - tt %*% K_run - K_run %*% tt + tt %*% K_run %*% tt
    Y_run <- Y_run - tt %*% Y_run

  }

  Y_hat <- K_c %*% U %*% solve(t(T) %*% K_c %*% U) %*% t(T) %*% Y_c + meanY

  end_time <- Sys.time()

  result <- list()
  result$n <- n
  result$p <- p
  result$nk <- nk
  result$T <- T
  result$U <- U
  result$K <- K
  result$K_c <- K_c
  result$Xfd_list <- Xfd_list
  result$XX_list <- XX_list
  result$Y_c <- Y_c
  result$meanY <- meanY
  result$Y_hat <- Y_hat
  result$obser_time <- obser_time
  result$basis <- basis
  result$sigm <- sigm
  result$CVscore <- CV_score
  result$time <- end_time - start_time

  class(result) <- 'KFPLS'

  return(result)


}
