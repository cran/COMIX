#' This function creates an object that summarizes the Geweke convergence 
#' diagnostic.
#'
#' @param res An object of class \code{COMIX} or \code{tidyChainCOMIX}.
#' @param params A character vector naming the parameters to compute the 
#' Geweke diagnostic for. 
#' @param frac1 Double, fraction to use from beginning of chain.
#' @param frac2 Double, fraction to use from end of chain.
#' @param probs A vector of 2 doubles, probabilities denoting the limits
#' of a confidence interval for the convergence diagnostic.
#' @return An \code{gewekeParamsCOMIX} object which is a named list,
#' with a named element for each requested parameter. Each element is a data
#' frame that includes the Geweke diagnostic and result of a stationarity test 
#' for the parameter.
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
#' gwk <- gewekeParams(tidy_chain, "w")
#' # (see vignette for a more detailed example)
#' @export
gewekeParams <- function(res, params = c("w", "xi", "xi0", "psi", "G", "E", "eta"), 
                         frac1 = 0.1, frac2 = 0.5, probs = c(0.025, 0.975)) {
    
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
  
  gewekeParams <- list()
  class(gewekeParams) <- "gewekeParamsCOMIX"
  attributes(gewekeParams)$n <- n
  attributes(gewekeParams)$p <- P
  attributes(gewekeParams)$nsave <- nsave
  attributes(gewekeParams)$K <- K
  attributes(gewekeParams)$J <- J
  attributes(gewekeParams)$non_trivial_k <- non_trivial_k
  attributes(gewekeParams)$non_triv_j_k <- non_triv_j_k
  attributes(gewekeParams)$frac <- c(frac1, frac2)
  attributes(gewekeParams)$glob_freq_t <- glob_freq_t
  
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
    gd <- geweke.diag(x = aa, frac1 = frac1, frac2 = frac2)
    kj <- apply(str_split(names(gd$z), "_", simplify = TRUE), 2, as.character)
    dkj <- tc %>% select(.data$k, .data$j) %>% distinct()
    stopifnot(all.equal(dkj %>% as.matrix() %>% unname(), kj))
    gewekeParams$w <- tibble(dkj, geweke = gd$z)

    dplyr.summarise.inform <- options()$dplyr.summarise.inform
    options(dplyr.summarise.inform = FALSE)
    meanW <- tc %>% group_by(.data$j, .data$k) %>% summarize(meanW = mean(.data$W))
    stopifnot(all(meanW %>% select(.data$k, .data$j) %>% as.matrix() == kj))
    options(dplyr.summarise.inform = dplyr.summarise.inform)
    
    gewekeParams$w <-
      gewekeParams$w %>% 
      mutate(stationary = .data$geweke > qnorm(probs[1]) & .data$geweke < qnorm(probs[2])) %>%
      left_join(meanW, by = c("j", "k"))
  }
  
  # xi0 -----
  if ("xi0" %in% params) {
    tc <- tidy_chain$xi0 %>% filter(.data$k %in% non_trivial_k)

    a <-
      tc %>%
      pivot_wider(names_from = c(.data$k, .data$p), values_from = c(.data$xi0)) %>%
      select(-.data$iter)
    
    aa <- mcmc(data = a, start = 1)
    gd <- geweke.diag(x = aa, frac1 = frac1, frac2 = frac2)
    kp <- apply(str_split(names(gd$z), "_", simplify = TRUE), 2, as.character)
    dpk <- tc %>% select(.data$k, .data$p) %>% distinct()
    stopifnot(all.equal(dpk %>% as.matrix() %>% unname(), kp))
    gewekeParams$xi0 <- tibble(dpk, geweke = gd$z)

    dplyr.summarise.inform <- options()$dplyr.summarise.inform
    options(dplyr.summarise.inform = FALSE)
    meanXi0 <- tc %>% group_by(.data$p, .data$k) %>% summarize(meanXi0 = mean(.data$xi0))
    options(dplyr.summarise.inform = dplyr.summarise.inform)
    
    gewekeParams$xi0 <-
      gewekeParams$xi0 %>% mutate(stationary = .data$geweke > qnorm(probs[1]) & .data$geweke < qnorm(probs[2])) %>%
      left_join(glob_freq_t %>% select(.data$k, .data$frq_t), by = "k") %>%
      left_join(meanXi0, by = c("p", "k"))
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
    gd <- geweke.diag(x = aa, frac1 = frac1, frac2 = frac2)
    kpj <- apply(str_split(names(gd$z), "_", simplify = TRUE), 2, as.character)
    dkpj <- tc %>% select(.data$k, .data$p, .data$j) %>% distinct()
    stopifnot(all.equal(dkpj %>% as.matrix() %>% unname(), kpj))
    gewekeParams$xi <- tibble(dkpj, geweke = gd$z)
    
    local_freq_t <- 
      local_freq_t %>% 
      ungroup() %>% 
      mutate(j = factor(.data$j), k = factor(.data$k), .data$frq_t) %>%
      select(.data$j, .data$k, .data$frq_t)
    
    dplyr.summarise.inform <- options()$dplyr.summarise.inform
    options(dplyr.summarise.inform = FALSE)
    meanXi <- tc %>% group_by(.data$j, .data$p, .data$k) %>% summarize(meanXi = mean(.data$xi))
    options(dplyr.summarise.inform = dplyr.summarise.inform)
    
    gewekeParams$xi <-
      gewekeParams$xi %>% 
      mutate(stationary = .data$geweke > qnorm(probs[1]) & .data$geweke < qnorm(probs[2])) %>%
      left_join(local_freq_t, by = c("j", "k")) %>%
      left_join(meanXi, by = c("j", "p", "k"))
    }
  
  # psi -----
  if ("psi" %in% params) {
    tc <- tidy_chain$psi %>% filter(.data$k %in% non_trivial_k)

    a <-
      tc %>%
      pivot_wider(names_from = c(.data$k, .data$p), values_from = c(.data$psi)) %>%
      select(-.data$iter)
    
    aa <- mcmc(data = a, start = 1)
    gd <- geweke.diag(x = aa, frac1 = frac1, frac2 = frac2)
    kp <- apply(str_split(names(gd$z), "_", simplify = TRUE), 2, as.character)
    dpk <- tc %>% select(.data$k, .data$p) %>% distinct()
    stopifnot(all.equal(dpk %>% as.matrix() %>% unname(), kp))
    gewekeParams$psi <- tibble(dpk, geweke = gd$z)
    
    dplyr.summarise.inform <- options()$dplyr.summarise.inform
    options(dplyr.summarise.inform = FALSE)
    meanPsi <- tc %>% group_by(.data$p, .data$k) %>% summarize(meanPsi = mean(.data$psi))
    options(dplyr.summarise.inform = dplyr.summarise.inform)
    
    gewekeParams$psi <-
      gewekeParams$psi %>% 
      mutate(stationary = .data$geweke > qnorm(probs[1]) & .data$geweke < qnorm(probs[2])) %>%
      left_join(glob_freq_t %>% select(.data$k, .data$frq_t), by = "k") %>%
      left_join(meanPsi, by = c("p", "k"))
  }
  
  # G -----
  if ("G" %in% params) {
    tc <- tidy_chain$G %>% filter(.data$k %in% non_trivial_k)

    a <-
      tc %>%
      pivot_wider(names_from = c(.data$k, .data$p1, .data$p2), values_from = c(.data$G)) %>%
      select(-.data$iter)
    
    aa <- mcmc(data = a, start = 1)
    gd <- geweke.diag(x = aa, frac1 = frac1, frac2 = frac2)
    kp1p2 <- apply(str_split(names(gd$z), "_", simplify = TRUE), 2, as.character)
    dkp1p2 <- tc %>% select(.data$k, .data$p1, .data$p2) %>% distinct()
    stopifnot(all.equal(dkp1p2 %>% as.matrix() %>% unname(), kp1p2))
    gewekeParams$G <- tibble(dkp1p2, geweke = gd$z)
    
    dplyr.summarise.inform <- options()$dplyr.summarise.inform
    options(dplyr.summarise.inform = FALSE)
    meanG <- tc %>% group_by(.data$k, .data$p1, .data$p2) %>% summarize(meanG = mean(.data$G))
    options(dplyr.summarise.inform = dplyr.summarise.inform)
    
    gewekeParams$G <-
      gewekeParams$G %>% 
      mutate(stationary = .data$geweke > qnorm(probs[1]) & .data$geweke < qnorm(probs[2])) %>%
      left_join(glob_freq_t %>% select(.data$k, .data$frq_t), by = "k") %>%
      left_join(meanG, by = c("k", "p1", "p2")) 
  }
  
  # E -----
  if ("E" %in% params) {
    tc <- tidy_chain$E %>% filter(.data$k %in% non_trivial_k)

    a <-
      tc %>%
      pivot_wider(names_from = c(.data$k, .data$p1, .data$p2), values_from = c(.data$E)) %>%
      select(-.data$iter)
    
    aa <- mcmc(data = a, start = 1)
    gd <- geweke.diag(x = aa, frac1 = frac1, frac2 = frac2)
    kp1p2 <- apply(str_split(names(gd$z), "_", simplify = TRUE), 2, as.character)
    dkp1p2 <- tc %>% select(.data$k, .data$p1, .data$p2) %>% distinct()
    stopifnot(all.equal(dkp1p2 %>% as.matrix() %>% unname(), kp1p2))
    gewekeParams$E <- tibble(dkp1p2, geweke = gd$z)
    
    dplyr.summarise.inform <- options()$dplyr.summarise.inform
    options(dplyr.summarise.inform = FALSE)
    meanE <- tc %>% group_by(.data$k, .data$p1, .data$p2) %>% summarize(meanE = mean(.data$E))
    options(dplyr.summarise.inform = dplyr.summarise.inform)
    
    gewekeParams$E <-
      gewekeParams$E %>% 
      mutate(stationary = .data$geweke > qnorm(probs[1]) & .data$geweke < qnorm(probs[2])) %>%
      left_join(glob_freq_t %>% select(.data$k, .data$frq_t), by = "k") %>%
      left_join(meanE, by = c("k", "p1", "p2")) 
  }
  
  # eta -----
  if ("eta" %in% params) {
    eta <- mcmc(data = tidy_chain$eta$eta, start = 1)
    gd <- geweke.diag(x = eta, frac1 = frac1, frac2 = frac2)
    gewekeParams$eta <- tibble(geweke = unname(gd$z))
    gewekeParams$eta <- 
      gewekeParams$eta %>% mutate(stationary = .data$geweke > qnorm(probs[1]) & .data$geweke < qnorm(probs[2]))
  }
    
  return(gewekeParams)
}

#' This function creates plots for the Geweke diagnostic and results of test of
#' stationarity for the parameters of the model. 
#'
#' @param gwk An object of class \code{gewekeParamsCOMIX} as created
#' by the function \code{gewekeParams}.
#' @param param Character, naming the parameter to create a plot of the Geweke
#' diagnostic for.
#' @return A \code{ggplot2} plot containing the Geweke diagnostic plot.
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
#' gwk <- gewekeParams(tidy_chain, "w")
#' plotGewekeParams(gwk, "w")
#' # (see vignette for a more detailed example)
#' @export
plotGewekeParams <- function(gwk, param) {
  stopifnot(class(gwk) == "gewekeParamsCOMIX")
  stopifnot(length(param) == 1)
  stopifnot(param %in% c("w", "xi0", "xi", "psi", "G", "E"))
  
  J <- attributes(gwk)$J
  
  j_names <- paste0("Sample ", 1:J)
  names(j_names) <- 1:J
  
  non_trivial_k <- attributes(gwk)$non_trivial_k
  glob_freq_t <- attributes(gwk)$glob_freq_t
  frq_t <- glob_freq_t$frq_t[glob_freq_t$k %in% non_trivial_k]
  k_names_frq <- paste0("Cluster ", non_trivial_k, "\n(Est. Freq. = ", round(frq_t, 2), ")")
  names(k_names_frq) <- non_trivial_k
  
  scm <- 
    scale_color_manual(
      name = "Geweke\nStationarity test", 
      labels = c("Passed", "Failed"), 
      values = c("TRUE" = "#00ba38", "FALSE" = "#f8766d")
    ) 
  
  # w -----
  if (param == "w") {
    g <-
      gwk$w %>% 
      ggplot(aes(x = .data$k, y = .data$meanW, color = .data$stationary)) + 
      geom_point()  +
      geom_segment(aes(xend = .data$k, y = 0, yend = .data$meanW)) + 
      scm +
      ylab("Estimated weight") + 
      xlab("Cluster Number") +
      facet_wrap(~ .data$j, labeller = labeller(j = j_names))
    
    return(g)
  }
  
  # xi0 -----
  if (param == "xi0") {
    g <-
      gwk$xi0 %>% 
      ggplot(aes(x = .data$p, y = .data$meanXi0, color = .data$stationary)) + 
      geom_segment(aes(xend = .data$p, y = 0, yend = .data$meanXi0)) + 
      geom_point() +
      scm + 
      ylab(expression(xi[0])) + 
      xlab("Margin") +
      facet_wrap(~ .data$k, labeller = labeller(k = k_names_frq))
    
    return(g)
  }
  
  
  # xi -----
  if (param == "xi") {
    g <-
      gwk$xi %>% 
      ggplot(aes(x = .data$p, y = .data$meanXi, color = .data$stationary)) + 
      geom_segment(aes(xend = .data$p, y = 0, yend = .data$meanXi)) + 
      geom_point() + 
      scm + 
      ylab(expression(xi)) +
      xlab("Margin") +
      facet_grid(.data$j ~ .data$k, labeller = labeller(j = j_names, k = k_names_frq))
    
    return(g)
  }
  
  
  # psi -----
  if (param == "psi") {
    g <-
      gwk$psi %>% 
      ggplot(aes(x = .data$p, y = .data$meanPsi, color = .data$stationary)) + 
      geom_segment(aes(xend = .data$p, y = 0, yend = .data$meanPsi)) + 
      geom_point() +
      scm + 
      ylab(expression(psi)) + 
      xlab("Margin") +
      facet_wrap(~ .data$k, labeller = labeller(k = k_names_frq))
    
    return(g)
  }
  
  
  # G -----
  if (param == "G") {
    g <-
      gwk$G %>% 
      mutate(
        Diagonal = 
          factor(
            .data$p1 == .data$p2,
            levels = c(TRUE, FALSE),
            labels = c("Diagonal", "Off-Diagonal")
          )
      ) %>%
      mutate(p1p2 = paste0(.data$p1, ", ", .data$p2)) %>%
      ggplot(aes(x = .data$p1p2, y = .data$meanG, color = .data$stationary)) + 
      geom_segment(aes(xend = .data$p1p2, y = 0, yend = .data$meanG)) + 
      geom_point(aes(shape = .data$Diagonal)) +
      scm + 
      ylab("G") + 
      xlab("Margin") +
      facet_wrap(~ .data$k, labeller = labeller(k = k_names_frq)) +
      theme(axis.text.x = element_text(angle = -90, vjust = .5, hjust = 1))
    
    return(g)
  }
  
  # E -----
  if (param == "E") {
    g <-
      gwk$E %>% 
      mutate(
        Diagonal = 
          factor(
            .data$p1 == .data$p2,
            levels = c(TRUE, FALSE),
            labels = c("Diagonal", "Off-Diagonal")
          )
      ) %>%
      mutate(p1p2 = paste0(.data$p1, ", ", .data$p2)) %>%
      ggplot(aes(x = .data$p1p2, y = .data$meanE, color = .data$stationary)) + 
      geom_segment(aes(xend = .data$p1p2, y = 0, yend = .data$meanE)) + 
      geom_point(aes(shape = .data$Diagonal)) +
      scm + 
      ylab("E") + 
      xlab("Margin") +
      facet_wrap(~ .data$k, labeller = labeller(k = k_names_frq)) +
      theme(axis.text.x = element_text(angle = -90, vjust = .5, hjust = 1))
    
    return(g)
  }
  
  return(NULL)
}