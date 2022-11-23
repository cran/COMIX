#' This function creates an object that summarizes the Heidelberg-Welch 
#' convergence diagnostic.
#'
#' @param res An object of class \code{COMIX} or \code{tidyChainCOMIX}.
#' @param params A character vector naming the parameters to compute the 
#' Heidelberg-Welch diagnostic for.
#' @param eps Target value for ratio of halfwidth to sample mean.
#' @param pvalue Significance level to use.
#' @return An \code{heidelParamsCOMIX} object which is a named list,
#' with a named element for each requested parameter. Each element is a data
#' frame that includes the Heidelberg-Welch diagnostic and results of a
#' stationarity test for the parameter.
#' @examples
#' library(COMIX)
#' # Number of observations for each sample (row) and cluster (column):
#' njk <- 
#'   matrix(
#'     c(
#'       150, 300,
#'       250, 200
#'     ),
#'     nrow = 2,
#'     byrow = TRUE
#'   )
#' 
#' # Dimension of data:
#' p <- 3
#' 
#' # Scale and skew parameters for first cluster:
#' Sigma1 <- matrix(0.5, nrow = p, ncol = p) + diag(0.5, nrow = p)
#' alpha1 <- rep(0, p)
#' alpha1[1] <- -5
#' # location parameter for first cluster in first sample:
#' xi11 <- rep(0, p)
#' # location parameter for first cluster in second sample (aligned with first):
#' xi21 <- rep(0, p)
#' 
#' # Scale and skew parameters for second cluster:
#' Sigma2 <- matrix(-1/3, nrow = p, ncol = p) + diag(1 + 1/3, nrow = p)
#' alpha2 <- rep(0, p)
#' alpha2[2] <- 5
#' # location parameter for second cluster in first sample:
#' xi12 <- rep(3, p)
#' # location parameter for second cluster in second sample (misaligned with first):
#' xi22 <- rep(4, p)
#' 
#' # Sample data:
#' set.seed(1)
#' Y <- 
#'   rbind(
#'     sn::rmsn(njk[1, 1], xi = xi11, Omega = Sigma1, alpha = alpha1),
#'     sn::rmsn(njk[1, 2], xi = xi12, Omega = Sigma2, alpha = alpha2),
#'     sn::rmsn(njk[2, 1], xi = xi21, Omega = Sigma1, alpha = alpha1),
#'     sn::rmsn(njk[2, 2], xi = xi22, Omega = Sigma2, alpha = alpha2)
#'   )
#' 
#' C <- c(rep(1, rowSums(njk)[1]), rep(2, rowSums(njk)[2]))
#' 
#' prior <- list(zeta = 1, K = 10)
#' pmc <- list(naprt = 5, nburn = 200, nsave = 200) # Reasonable usage
#' pmc <- list(naprt = 5, nburn = 2, nsave = 5) # Minimal usage for documentation
#' # Fit the model:
#' res <- comix(Y, C, pmc = pmc, prior = prior)
#' 
#' # Relabel to resolve potential label switching issues:
#' res_relab <- relabelChain(res)
#' effssz <- effectiveSampleSize(res_relab, "w")
#' # Or:
#' tidy_chain <- tidyChain(res_relab, "w")
#' hd <- heidelParams(tidy_chain, "w")
#' # (see vignette for a more detailed example)
#' @export
heidelParams <- function(res, params = c("w", "xi", "xi0", "psi", "G", "E", "eta"), 
                         eps = 0.1, pvalue = 0.05) {
  
  stopifnot(is(res, "COMIX") | is(res, "tidyChainCOMIX"))
  
  if (is(res, "COMIX")) {
    tidy_chain <- tidyChain(res, params)
  } else {
    tidy_chain <- res
  }
  
  n <- attributes(tidy_chain)$n
  P <- attributes(tidy_chain)$p
  nsave <- attributes(tidy_chain)$nsave
  K <- attributes(tidy_chain)$K
  J <- attributes(tidy_chain)$J
  
  non_trivial_k <- attributes(tidy_chain)$non_trivial_k
  non_triv_j_k <- attributes(tidy_chain)$non_triv_j_k

  glob_freq_t <- attributes(tidy_chain)$glob_freq_t
  local_freq_t <- attributes(tidy_chain)$local_freq_t
  
  heidelParams <- list()
  class(heidelParams) <- "heidelParamsCOMIX"
  attributes(heidelParams)$n <- n
  attributes(heidelParams)$p <- P
  attributes(heidelParams)$nsave <- nsave
  attributes(heidelParams)$K <- K
  attributes(heidelParams)$J <- J
  attributes(heidelParams)$non_trivial_k <- non_trivial_k
  attributes(heidelParams)$non_triv_j_k <- non_triv_j_k
  attributes(heidelParams)$eps <- eps
  attributes(heidelParams)$pvalue <- pvalue
  attributes(heidelParams)$glob_freq_t <- glob_freq_t
  
  # w -----
  if ("w" %in% params) {
    tc <- tidy_chain$w
    tc$triv <- TRUE
    for (j in 1:J) {
      tc$triv[tc$j == j & tc$k %in% non_triv_j_k[[as.character(j)]]] <- FALSE
    }
    tc <- tc %>% filter(!.data$triv) %>% select(-.data$triv)

    a <- 
      tc %>% 
      pivot_wider(names_from = c(.data$k, .data$j), values_from = c(.data$W)) %>%
      select(-.data$iter)
    aa <- mcmc(data = a, start = 1)
    hd <- heidel.diag(x = aa, eps = eps, pvalue = pvalue)
    class(hd) <- "matrix"
    kj <- apply(str_split(rownames(hd), "_", simplify = TRUE), 2, as.character)
    dkj <- tc %>% select(.data$k, .data$j) %>% distinct()
    stopifnot(all.equal(dkj %>% as.matrix() %>% unname(), kj))
    
    heidelParams$w <- 
      tibble(dkj, as_tibble(hd)) %>%
      mutate(kj = rownames(hd)) %>%
      mutate(
        stest = factor(.data$stest, levels = c(1, 0), labels = c("passed", "failed")),
        htest = factor(.data$htest, levels = c(1, 0), labels = c("passed", "failed"))
        )
    
    dplyr.summarise.inform <- options()$dplyr.summarise.inform
    options(dplyr.summarise.inform = FALSE)
    attributes(heidelParams$w)$meanW <- 
      tidy_chain$w %>% 
      group_by(.data$j, .data$k) %>% 
      summarize(meanW = mean(.data$W))
    options(dplyr.summarise.inform = dplyr.summarise.inform)
    
  }
  
  # xi0 -----
  if ("xi0" %in% params) {
    tc <- tidy_chain$xi0 %>% filter(.data$k %in% non_trivial_k)

    a <-
      tc %>%
      pivot_wider(names_from = c(.data$k, .data$p), values_from = c(.data$xi0)) %>%
      select(-.data$iter)
    
    aa <- mcmc(data = a, start = 1)
    hd <- heidel.diag(x = aa, eps = eps, pvalue = pvalue)
    class(hd) <- "matrix"
    kp <- apply(str_split(rownames(hd), "_", simplify = TRUE), 2, as.character)
    dkp <- tc %>% select(.data$k, .data$p) %>% distinct()
    stopifnot(all.equal(dkp %>% as.matrix() %>% unname(), kp))
    
    heidelParams$xi0 <- 
      tibble(dkp, as_tibble(hd)) %>%
      mutate(kp = rownames(hd)) %>%
      mutate(
        stest = factor(.data$stest, levels = c(1, 0), labels = c("passed", "failed")),
        htest = factor(.data$htest, levels = c(1, 0), labels = c("passed", "failed"))
      ) %>%
      left_join(glob_freq_t %>% select(.data$k, .data$frq_t), by = "k")
    
    dplyr.summarise.inform <- options()$dplyr.summarise.inform
    options(dplyr.summarise.inform = FALSE)
    attributes(heidelParams$xi0)$meanXi0 <- 
      tidy_chain$xi0 %>% 
      group_by(.data$k, .data$p) %>% 
      summarize(meanXi0 = mean(.data$xi0))
    options(dplyr.summarise.inform = dplyr.summarise.inform)
  }
  
  
  # xi -----
  if ("xi" %in% params) {
    tc <- tidy_chain$xi
    tc$triv <- TRUE
    for (j in 1:J) {
      tc$triv[tc$j == j & tc$k %in% non_triv_j_k[[as.character(j)]]] <- FALSE
    }
    tc <- tc %>% filter(!.data$triv) %>% select(-.data$triv)

    a <-
      tc %>%
      pivot_wider(names_from = c(.data$k, .data$p, .data$j), values_from = c(.data$xi)) %>%
      select(-.data$iter)
    
    aa <- mcmc(data = a, start = 1)
    hd <- heidel.diag(x = aa, eps = eps, pvalue = pvalue)
    class(hd) <- "matrix"
    kpj <- apply(str_split(rownames(hd), "_", simplify = TRUE), 2, as.character)
    dkpj <- tc %>% select(.data$k, .data$p, .data$j) %>% distinct()
    stopifnot(all.equal(dkpj %>% as.matrix() %>% unname(), kpj))
    
    local_freq_t <- 
      local_freq_t %>% 
      ungroup() %>% 
      mutate(j = factor(.data$j), k = factor(.data$k), .data$frq_t) %>%
      select(.data$j, .data$k, .data$frq_t)
    
    heidelParams$xi <- 
      tibble(dkpj, as_tibble(hd)) %>%
      mutate(kpj = rownames(hd)) %>%
      mutate(
        stest = factor(.data$stest, levels = c(1, 0), labels = c("passed", "failed")),
        htest = factor(.data$htest, levels = c(1, 0), labels = c("passed", "failed"))
      ) %>%
      left_join(local_freq_t, by = c("j", "k"))
    
    dplyr.summarise.inform <- options()$dplyr.summarise.inform
    options(dplyr.summarise.inform = FALSE)
    attributes(heidelParams$xi)$meanXi <- 
      tidy_chain$xi %>% 
      group_by(.data$k, .data$j, .data$p) %>% 
      summarize(meanXi = mean(.data$xi))
    options(dplyr.summarise.inform = dplyr.summarise.inform)
  }
  
  # psi -----
  if ("psi" %in% params) {
    tc <- tidy_chain$psi %>% filter(.data$k %in% non_trivial_k)

    a <-
      tc %>%
      pivot_wider(names_from = c(.data$k, .data$p), values_from = c(.data$psi)) %>%
      select(-.data$iter)
    
    aa <- mcmc(data = a, start = 1)
    hd <- heidel.diag(x = aa, eps = eps, pvalue = pvalue)
    class(hd) <- "matrix"
    kp <- apply(str_split(rownames(hd), "_", simplify = TRUE), 2, as.character)
    dkp <- tc %>% select(.data$k, .data$p) %>% distinct()
    stopifnot(all.equal(dkp %>% as.matrix() %>% unname(), kp))
    
    heidelParams$psi <- 
      tibble(dkp, as_tibble(hd)) %>%
      mutate(kp = rownames(hd)) %>%
      mutate(
        stest = factor(.data$stest, levels = c(1, 0), labels = c("passed", "failed")),
        htest = factor(.data$htest, levels = c(1, 0), labels = c("passed", "failed"))
      ) %>%
      left_join(glob_freq_t %>% select(.data$k, .data$frq_t), by = "k")
    
    dplyr.summarise.inform <- options()$dplyr.summarise.inform
    options(dplyr.summarise.inform = FALSE)
    attributes(heidelParams$psi)$meanPsi <- 
      tidy_chain$psi %>% 
      group_by(.data$k, .data$p) %>% 
      summarize(meanPsi = mean(.data$psi))
    options(dplyr.summarise.inform = dplyr.summarise.inform)
  }
  
  # G -----
  if ("G" %in% params) {
    tc <- tidy_chain$G %>% filter(.data$k %in% non_trivial_k)

    a <-
      tc %>%
      pivot_wider(names_from = c(.data$k, .data$p1, .data$p2), values_from = c(.data$G)) %>%
      select(-.data$iter)
    
    aa <- mcmc(data = a, start = 1)
    hd <- heidel.diag(x = aa, eps = eps, pvalue = pvalue)
    class(hd) <- "matrix"
    kp1p2 <- apply(str_split(rownames(hd), "_", simplify = TRUE), 2, as.character)
    dkp1p2 <- tc %>% select(.data$k, .data$p1, .data$p2) %>% distinct()
    stopifnot(all.equal(dkp1p2 %>% as.matrix() %>% unname(), kp1p2))
    
    heidelParams$G <-
      tibble(dkp1p2, as_tibble(hd)) %>%
      mutate(p1p2 = rownames(hd)) %>%
      mutate(
        stest = factor(.data$stest, levels = c(1, 0), labels = c("passed", "failed")),
        htest = factor(.data$htest, levels = c(1, 0), labels = c("passed", "failed"))
      ) %>%
      left_join(glob_freq_t %>% select(.data$k, .data$frq_t), by = "k")
    
    dplyr.summarise.inform <- options()$dplyr.summarise.inform
    options(dplyr.summarise.inform = FALSE)
    attributes(heidelParams$G)$meanG <- 
      tidy_chain$G %>% 
      group_by(.data$k, .data$p1, .data$p2) %>% 
      summarize(meanG = mean(.data$G))
    options(dplyr.summarise.inform = dplyr.summarise.inform)
  }
  
  # E -----
  if ("E" %in% params) {
    tc <- tidy_chain$E %>% filter(.data$k %in% non_trivial_k)

    a <-
      tc %>%
      pivot_wider(names_from = c(.data$k, .data$p1, .data$p2), values_from = c(.data$E)) %>%
      select(-.data$iter)
    
    aa <- mcmc(data = a, start = 1)
    hd <- heidel.diag(x = aa, eps = eps, pvalue = pvalue)
    class(hd) <- "matrix"
    kp1p2 <- apply(str_split(rownames(hd), "_", simplify = TRUE), 2, as.character)
    dkp1p2 <- tc %>% select(.data$k, .data$p1, .data$p2) %>% distinct()
    stopifnot(all.equal(dkp1p2 %>% as.matrix() %>% unname(), kp1p2))
    
    heidelParams$E <-
      tibble(dkp1p2, as_tibble(hd)) %>%
      mutate(p1p2 = rownames(hd)) %>%
      mutate(
        stest = factor(.data$stest, levels = c(1, 0), labels = c("passed", "failed")),
        htest = factor(.data$htest, levels = c(1, 0), labels = c("passed", "failed"))
      ) %>%
      left_join(glob_freq_t %>% select(.data$k, .data$frq_t), by = "k")
    
    dplyr.summarise.inform <- options()$dplyr.summarise.inform
    options(dplyr.summarise.inform = FALSE)
    attributes(heidelParams$E)$meanE <- 
      tidy_chain$E %>% 
      group_by(.data$k, .data$p1, .data$p2) %>% 
      summarize(meanE = mean(.data$E))
    options(dplyr.summarise.inform = dplyr.summarise.inform)
  }
  
  # eta -----
  if ("eta" %in% params) {
    eta <- mcmc(data = tidy_chain$eta$eta, start = 1)
    hd <- heidel.diag(x = eta, eps = eps, pvalue = pvalue)
    class(hd) <- "matrix"
    
    heidelParams$eta <-
      as_tibble(hd) %>%
      mutate(
        stest = factor(.data$stest, levels = c(1, 0), labels = c("passed", "failed")),
        htest = factor(.data$htest, levels = c(1, 0), labels = c("passed", "failed"))
      )
  }
    
  return(heidelParams)
}

#' This function creates plots for the Heidelberg-Welch diagnostic and 
#' results of test of stationarity for the parameters of the model. 
#'
#' @param hd An object of class \code{heidelParamsCOMIX} as created
#' by the function \code{heidelParams}.
#' @param param Character, naming the parameter to create a plot of the 
#' Heidelberg-Welch diagnostic for.
#' @return A \code{ggplot2} plot containing the Heidelberg-Welch diagnostic plot.
#' @examples
#' library(COMIX)
#' # Number of observations for each sample (row) and cluster (column):
#' njk <- 
#'   matrix(
#'     c(
#'       150, 300,
#'       250, 200
#'     ),
#'     nrow = 2,
#'     byrow = TRUE
#'   )
#' 
#' # Dimension of data:
#' p <- 3
#' 
#' # Scale and skew parameters for first cluster:
#' Sigma1 <- matrix(0.5, nrow = p, ncol = p) + diag(0.5, nrow = p)
#' alpha1 <- rep(0, p)
#' alpha1[1] <- -5
#' # location parameter for first cluster in first sample:
#' xi11 <- rep(0, p)
#' # location parameter for first cluster in second sample (aligned with first):
#' xi21 <- rep(0, p)
#' 
#' # Scale and skew parameters for second cluster:
#' Sigma2 <- matrix(-1/3, nrow = p, ncol = p) + diag(1 + 1/3, nrow = p)
#' alpha2 <- rep(0, p)
#' alpha2[2] <- 5
#' # location parameter for second cluster in first sample:
#' xi12 <- rep(3, p)
#' # location parameter for second cluster in second sample (misaligned with first):
#' xi22 <- rep(4, p)
#' 
#' # Sample data:
#' set.seed(1)
#' Y <- 
#'   rbind(
#'     sn::rmsn(njk[1, 1], xi = xi11, Omega = Sigma1, alpha = alpha1),
#'     sn::rmsn(njk[1, 2], xi = xi12, Omega = Sigma2, alpha = alpha2),
#'     sn::rmsn(njk[2, 1], xi = xi21, Omega = Sigma1, alpha = alpha1),
#'     sn::rmsn(njk[2, 2], xi = xi22, Omega = Sigma2, alpha = alpha2)
#'   )
#' 
#' C <- c(rep(1, rowSums(njk)[1]), rep(2, rowSums(njk)[2]))
#' 
#' prior <- list(zeta = 1, K = 10)
#' pmc <- list(naprt = 5, nburn = 200, nsave = 200) # Reasonable usage
#' pmc <- list(naprt = 5, nburn = 2, nsave = 5) # Minimal usage for documentation
#' # Fit the model:
#' res <- comix(Y, C, pmc = pmc, prior = prior)
#' 
#' # Relabel to resolve potential label switching issues:
#' res_relab <- relabelChain(res)
#' effssz <- effectiveSampleSize(res_relab, "w")
#' # Or:
#' tidy_chain <- tidyChain(res_relab, "w")
#' hd <- heidelParams(tidy_chain, "w")
#' plotHeidelParams(hd, "w")
#' # (see vignette for a more detailed example)
#' @export
plotHeidelParams <- function(hd, param) {
  stopifnot(is(hd, "heidelParamsCOMIX"))
  stopifnot(length(param) == 1)
  stopifnot(param %in% c("w", "xi0", "xi", "psi", "G", "E"))
  
  J <- attributes(hd)$J
  
  j_names <- paste0("Sample ", 1:J)
  names(j_names) <- 1:J

  non_trivial_k <- attributes(hd)$non_trivial_k
  glob_freq_t <- attributes(hd)$glob_freq_t
  frq_t <- glob_freq_t$frq_t[glob_freq_t$k %in% non_trivial_k]
  k_names_frq <- paste0("Cluster ", non_trivial_k, "\n(Est. Freq. = ", round(frq_t, 2), ")")
  names(k_names_frq) <- non_trivial_k
  
  scm <- 
    scale_color_manual(
      name = "Heidelberg-Welch\nStationarity test", 
      labels = c("Passed", "Failed"), 
      values = c("passed" = "#00ba38", "failed" = "#f8766d")
      ) 
    
  # w -----
  if (param == "w") {
    g <-
      hd$w %>% 
      left_join(attributes(hd$w)$meanW, by = c("k", "j")) %>%
      mutate(mean_na_replace = ifelse(!is.na(.data$mean), .data$mean, .data$meanW)) %>% 
      mutate(start_na_replace = ifelse(!is.na(.data$start), as.character(.data$start), "")) %>% 
      ggplot(aes(x = .data$k, y = .data$mean_na_replace, color = .data$stest, label = .data$start_na_replace)) + 
      geom_point()  +
      geom_segment(aes(xend = .data$k, y = 0, yend = .data$mean_na_replace)) + 
      geom_text(
        aes(y = .data$meanW + 0.13 * (max(.data$meanW) - min(.data$meanW)) * sign(.data$meanW)),
        size = 2.5, 
        color = "black"
        ) +
      scm +
      ylab("Estimated weight\n(start chain from)") + 
      xlab("Cluster Number") +
      facet_wrap(~ .data$j, labeller = labeller(j = j_names))
      
    return(g)
  }
  
  # xi0 -----
  if (param == "xi0") {
    g <-
      hd$xi0 %>% 
      left_join(attributes(hd$xi0)$meanXi0, by = c("k", "p")) %>%
      mutate(mean_na_replace = ifelse(!is.na(.data$mean), .data$mean, .data$meanXi0)) %>% 
      mutate(start_na_replace = ifelse(!is.na(.data$start), as.character(.data$start), "")) %>% 
      ggplot(aes(x = .data$p, y = .data$mean_na_replace, color = .data$stest, label = .data$start_na_replace)) + 
      geom_segment(aes(xend = .data$p, y = 0, yend = .data$mean_na_replace)) + 
      geom_point() +
      geom_text(
        aes(y = .data$meanXi0 + 0.13 * (max(.data$meanXi0) - min(.data$meanXi0)) * sign(.data$meanXi0)),
        size = 2.5, 
        color = "black"
      ) +
      scm + 
      ylab("Estimated grand location\n(start chain from)") + 
      xlab("Margin") +
      facet_wrap(~ .data$k, labeller = labeller(k = k_names_frq))
  
    return(g)
  }

  
  # xi -----
  if (param == "xi") {
    g <-
      hd$xi %>% 
      left_join(attributes(hd$xi)$meanXi, by = c("k", "j", "p")) %>%
      mutate(mean_na_replace = ifelse(!is.na(.data$mean), .data$mean, .data$meanXi)) %>% 
      mutate(start_na_replace = ifelse(!is.na(.data$start), as.character(.data$start), "")) %>% 
      ggplot(aes(x = .data$p, y = .data$mean_na_replace, color = .data$stest, label = .data$start_na_replace)) + 
      geom_segment(aes(xend = .data$p, y = 0, yend = .data$mean_na_replace)) + 
      geom_point() + 
      geom_text(
        aes(y = .data$meanXi + 0.13 * (max(.data$meanXi) - min(.data$meanXi)) * sign(.data$meanXi)),
        size = 2.5, 
        color = "black"
      ) +
      scm + 
      ylab("Estimated cluster-specific location\n(start chain from)") +
      xlab("Margin") +
      facet_grid(.data$j ~ .data$k, labeller = labeller(j = j_names, k = k_names_frq))
      
    return(g)
  }
  
  
  # psi -----
  if (param == "psi") {
    g <-
      hd$psi %>% 
      left_join(attributes(hd$psi)$meanPsi, by = c("k", "p")) %>%
      mutate(mean_na_replace = ifelse(!is.na(.data$mean), .data$mean, .data$meanPsi)) %>% 
      mutate(start_na_replace = ifelse(!is.na(.data$start), as.character(.data$start), "")) %>% 
      ggplot(aes(x = .data$p, y = .data$mean_na_replace, color = .data$stest, label = .data$start_na_replace)) + 
      geom_segment(aes(xend = .data$p, y = 0, yend = .data$mean_na_replace)) + 
      geom_point() +
      geom_text(
        aes(y = .data$meanPsi + 0.13 * (max(.data$meanPsi) - min(.data$meanPsi)) * sign(.data$meanPsi)),
        size = 2.5, 
        color = "black"
      ) +
      scm + 
      ylab("Estimated \U03C8\n(start chain from)") + 
      xlab("Margin") +
      facet_wrap(~ .data$k, labeller = labeller(k = k_names_frq))
    
    return(g)
  }
  
  
  # G -----
  if (param == "G") {
    g <-
      hd$G %>% 
      mutate(
        Diagonal = 
          factor(
            .data$p1 == .data$p2,
            levels = c(TRUE, FALSE),
            labels = c("Diagonal", "Off-Diagonal")
          )
      ) %>%
      mutate(p1p2 = paste0(.data$p1, ", ", .data$p2)) %>%
      left_join(attributes(hd$G)$meanG, by = c("k", "p1", "p2")) %>%
      mutate(mean_na_replace = ifelse(!is.na(.data$mean), .data$mean, .data$meanG)) %>% 
      mutate(start_na_replace = ifelse(!is.na(.data$start), as.character(.data$start), "")) %>% 
      ggplot(aes(x = .data$p1p2, y = .data$mean_na_replace, color = .data$stest, label = .data$start_na_replace)) + 
      geom_segment(aes(xend = .data$p1p2, y = 0, yend = .data$mean_na_replace)) + 
      geom_point(aes(shape = .data$Diagonal)) +
      geom_text(
        aes(y = .data$meanG + 0.13 * (max(.data$meanG) - min(.data$meanG)) * sign(.data$meanG)),
        size = 2.5, 
        color = "black"
      ) +
      scm + 
      ylab("Estimated G\n(start chain from)") + 
      xlab("Margin") +
      facet_wrap(~ .data$k, labeller = labeller(k = k_names_frq)) +
      theme(axis.text.x = element_text(angle = -90, vjust = .5, hjust = 1))
    
    return(g)
  }
  
  # E -----
  if (param == "E") {
    g <-
      hd$E %>% 
      mutate(
        Diagonal = 
          factor(
            .data$p1 == .data$p2,
            levels = c(TRUE, FALSE),
            labels = c("Diagonal", "Off-Diagonal")
          )
      ) %>%
      mutate(p1p2 = paste0(.data$p1, ", ", .data$p2)) %>%
      left_join(attributes(hd$E)$meanE, by = c("k", "p1", "p2")) %>%
      mutate(mean_na_replace = ifelse(!is.na(.data$mean), .data$mean, .data$meanE)) %>% 
      mutate(start_na_replace = ifelse(!is.na(.data$start), as.character(.data$start), "")) %>% 
      ggplot(aes(x = .data$p1p2, y = .data$mean_na_replace, color = .data$stest, label = .data$start_na_replace)) + 
      geom_segment(aes(xend = .data$p1p2, y = 0, yend = .data$mean_na_replace)) + 
      geom_point(aes(shape = .data$Diagonal)) +
      geom_text(
        aes(y = .data$meanE + 0.13 * (max(.data$meanE) - min(.data$meanE)) * sign(.data$meanE)),
        size = 2.5,
        color = "black"
      ) +
      scm + 
      ylab("Estimated G\n(start chain from)") + 
      xlab("Margin") +
      facet_wrap(~ .data$k, labeller = labeller(k = k_names_frq)) +
      theme(axis.text.x = element_text(angle = -90, vjust = .5, hjust = 1))
    
    return(g)
  }

  return(NULL)  
}
