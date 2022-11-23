#' This function creates an object that summarizes the effective sample size
#' for the parameters of the model. 
#'
#' @param res An object of class \code{COMIX} or \code{tidyChainCOMIX}.
#' @param params A character vector naming the parameters to compute the 
#' effective sample size for. 
#' @return An \code{effectiveSampleSizeCOMIX} object which is a named list,
#' with a named element for each requested parameter. Each element is a data
#' frame that includes the effective sample size for the parameter.
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
#' effssz <- effectiveSampleSize(tidy_chain, "w")
#' # (see vignette for a more detailed example)
#' @export
effectiveSampleSize <- function(res, params = c("w", "xi", "xi0", "psi", "G", "E", "eta")) {
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
  
  effssz <- list()
  class(effssz) <- "effectiveSampleSizeCOMIX"
  attributes(effssz)$n <- n
  attributes(effssz)$p <- P
  attributes(effssz)$nsave <- nsave
  attributes(effssz)$K <- K
  attributes(effssz)$J <- J
  attributes(effssz)$non_trivial_k <- non_trivial_k
  attributes(effssz)$non_triv_j_k <- non_triv_j_k
  attributes(effssz)$glob_freq_t <- attributes(tidy_chain)$glob_freq_t
  attributes(effssz)$local_freq_t <- attributes(tidy_chain)$local_freq_t
  
  # w -----
  if ("w" %in% params) {
    effssz_w <- 
      tibble(
        k = rep(1:K, times = J),
        j = rep(1:J, each = K),
        effssz = 0,
        triv = TRUE,
        W = 0
      )
    
    dplyr.summarise.inform <- options()$dplyr.summarise.inform
    options(dplyr.summarise.inform = FALSE)
    effssz_w$W <- 
      tidy_chain$w %>% 
      group_by(.data$k, .data$j) %>% 
      summarize(W = mean(.data$W)) %>% 
      arrange(.data$j, .data$k) %>% 
      pull(.data$W)
    options(dplyr.summarise.inform = dplyr.summarise.inform)
    
    a <- 
      tidy_chain$w %>% 
      pivot_wider(names_from = c(.data$k, .data$j), values_from = c(.data$W)) %>%
      select(-.data$iter)
    
    aa <- mcmc(data = a, start = 1)
    eff_raw <- effectiveSize(aa)
    kj <- apply(str_split(names(eff_raw), "_", simplify = TRUE), 2, as.integer)
    stopifnot(all.equal(effssz_w %>% select(.data$k, .data$j) %>% as.matrix() %>% unname(), kj))
    effssz_w$effssz <- eff_raw
    for (j in 1:J) {
      effssz_w$triv[effssz_w$j == j][non_triv_j_k[[as.character(j)]]] <- FALSE
    }
    
    effssz$w <- effssz_w
  }
  
  
  # xi0 -----
  if ("xi0" %in% params) {
    effssz_xi0 <- 
      tibble(
        k = rep(1:K, times = P),
        p = rep(1:P, each = K),
        effssz = 0,
        triv = TRUE
      )
  
    a <- 
      tidy_chain$xi0 %>% 
      pivot_wider(names_from = c(.data$k, .data$p), values_from = c(.data$xi0)) %>%
      select(-.data$iter)
    
    aa <- mcmc(data = a, start = 1)
    eff_raw <- effectiveSize(aa)
    kp <- apply(str_split(names(eff_raw), "_", simplify = TRUE), 2, as.integer)
    stopifnot(all.equal(effssz_xi0 %>% select(.data$k, .data$p) %>% as.matrix() %>% unname(), kp))
    effssz_xi0$effssz <- eff_raw
    effssz_xi0$triv[effssz_xi0$k %in% non_trivial_k] <- FALSE

    effssz$xi0 <- effssz_xi0
  }
  
  # xi -----
  if ("xi" %in% params) {
    effssz_xi <- 
      tibble(
        expand_grid(k = 1:K, p = 1:P, j = 1:J),
        effssz = 0,
        triv = TRUE
      )
    
    a <- 
      tidy_chain$xi %>% 
      pivot_wider(names_from = c(.data$k, .data$p, .data$j), values_from = c(.data$xi)) %>%
      select(-.data$iter)
    
    aa <- mcmc(data = a, start = 1)
    eff_raw <- effectiveSize(aa)
    kpj <- apply(str_split(names(eff_raw), "_", simplify = TRUE), 2, as.integer)
    stopifnot(all.equal(effssz_xi %>% select(.data$k, .data$p, .data$j) %>% as.matrix() %>% unname(), kpj))
    effssz_xi$effssz <- eff_raw
    for (j in 1:J) {
      effssz_xi$triv[effssz_xi$j == j & effssz_xi$k %in% non_triv_j_k[[as.character(j)]]] <- FALSE
    }

    effssz$xi <- effssz_xi
  }

  # psi -----
  if ("psi" %in% params) {
    effssz_psi <- 
      tibble(
        k = rep(1:K, times = P),
        p = rep(1:P, each = K),
        effssz = 0,
        triv = TRUE
      )
    
    a <- 
      tidy_chain$psi %>% 
      pivot_wider(names_from = c(.data$k, .data$p), values_from = c(.data$psi)) %>%
      select(-.data$iter)
    
    aa <- mcmc(data = a, start = 1)
    eff_raw <- effectiveSize(aa)
    kp <- apply(str_split(names(eff_raw), "_", simplify = TRUE), 2, as.integer)
    stopifnot(all.equal(effssz_psi %>% select(.data$k, .data$p) %>% as.matrix() %>% unname(), kp))
    effssz_psi$effssz <- eff_raw
    effssz_psi$triv[effssz_psi$k %in% non_trivial_k] <- FALSE
    
    effssz$psi <- effssz_psi
  }

  # G -----
  if ("G" %in% params) {
    a <- 
      tidy_chain$G %>% 
      pivot_wider(names_from = c(.data$k, .data$p1, .data$p2), values_from = c(.data$G)) %>%
      select(-.data$iter)
    
    aa <- mcmc(data = a, start = 1)
    eff_raw <- effectiveSize(aa)
    kpp <- apply(str_split(names(eff_raw), "_", simplify = TRUE), 2, as.integer)
    colnames(kpp) <- c("k", "p1", "p2")
    effssz_G <- 
      tibble(
        as_tibble(kpp),
        effssz = eff_raw,
        triv = TRUE
      )
    effssz_G$triv[effssz_G$k %in% non_trivial_k] <- FALSE
    
    effssz$G <- effssz_G
  }

  
  # E -----
  if ("E" %in% params) {
    a <- 
      tidy_chain$E %>% 
      pivot_wider(names_from = c(.data$k, .data$p1, .data$p2), values_from = c(.data$E)) %>%
      select(-.data$iter)
    
    aa <- mcmc(data = a, start = 1)
    eff_raw <- effectiveSize(aa)
    kpp <- apply(str_split(names(eff_raw), "_", simplify = TRUE), 2, as.integer)
    colnames(kpp) <- c("k", "p1", "p2")
    effssz_E <- 
      tibble(
        as_tibble(kpp),
        effssz = eff_raw,
        triv = TRUE
      )
    effssz_E$triv[effssz_E$k %in% non_trivial_k] <- FALSE
    
    effssz$E <- effssz_E
  }

  # eta -----
  if ("eta" %in% params) {
    eta <- pull(tidy_chain$eta)
    effssz$eta <- unname(effectiveSize(mcmc(eta)))
  }
  
  return(effssz)
}

#' This function creates plots for the effective sample size 
#' for the parameters of the model. 
#'
#' @param effssz An object of class \code{effectiveSampleSizeCOMIX} as created
#' by the function \code{effectiveSampleSize}.
#' @param param Character, naming the parameter to create a plot of effective
#' sample sizes.
#' @return A \code{ggplot2} plot containing the effective sample size plot.
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
#' effssz <- effectiveSampleSize(tidy_chain, "w")
#' plotEffectiveSampleSize(effssz, "w")
#' # (see vignette for a more detailed example)
#' @export
plotEffectiveSampleSize <- function(effssz, param) {
  stopifnot(is(effssz, "effectiveSampleSizeCOMIX"))
  stopifnot(length(param) == 1)
  stopifnot(param %in% c("w", "xi0", "xi", "psi", "G", "E"))
  
  J <- attributes(effssz)$J
  P <- attributes(effssz)$p
  nsave <- attributes(effssz)$nsave

  local_freq_t <- attributes(effssz)$local_freq_t
  glob_freq_t <- attributes(effssz)$glob_freq_t
  non_trivial_k <- attributes(effssz)$non_trivial_k
  non_triv_j_k <- attributes(effssz)$non_triv_j_k
  
  k_names <- paste0("Cluster ", non_trivial_k)
  names(k_names) <- non_trivial_k

  j_names <- paste0("Sample ", names(non_triv_j_k))
  names(j_names) <- names(non_triv_j_k)
  
  # w -----
  if (param == "w") {
    ggeffssz_w <- 
      effssz$w %>% 
      filter(!.data$triv) %>% 
      mutate(k = factor(.data$k), j = factor(.data$j)) %>%
      mutate(effssz = ifelse(.data$effssz <= nsave, .data$effssz, nsave))

    g <-    
      ggplot(ggeffssz_w) +
        geom_bar(
          mapping = 
            aes(
              x = .data$k, 
              y = .data$effssz,
            ),
          fill = rgb(1, 0, 0, 0.6),
          stat = "identity"
        ) +
      geom_bar(
        mapping = 
          aes(
            x = .data$k, 
            y = .data$W * nsave,
          ),
        fill = rgb(0, 0, 1, 0.6),
        stat = "identity"
      ) +
      scale_y_continuous(limits = c(0, nsave), sec.axis = ~ . / nsave) +
      facet_wrap(~ .data$j, labeller = labeller(j = j_names)) +
      ylab("Effective Sample Size (Red, Left Axis)\nEstimated Weight of Cluster (Blue, Right Axis)") +
      xlab("(Non-empty) Cluster number") +
      theme(legend.position = "none")
    return(g)
  }

  # xi0 -----  
  if ("xi0" == param) {
    ggeffssz_xi0 <- 
      effssz$xi0 %>% 
      filter(!.data$triv) %>%
      mutate(k = factor(.data$k)) %>%
      mutate(effssz = ifelse(.data$effssz <= nsave, .data$effssz, nsave))

    g <-
      ggplot(ggeffssz_xi0) +
        geom_bar(
          mapping =
            aes(
              x = .data$p,
              y = .data$effssz
            ),
          fill = rgb(1, 0, 0, 0.6),
          stat = "identity"
        ) +
        geom_bar(
          data = glob_freq_t,
          mapping = aes(
            x = .data$x,
            y = .data$frq_t * nsave
          ),
          stat = "identity",
          width = P - 0.1,
          fill = rgb(0, 0, 1, 0.6)
        ) +
        scale_y_continuous(limits = c(0, nsave), sec.axis = ~ . / nsave) +
        facet_wrap(~ .data$k, labeller = labeller(k = k_names)) +
        ylab("Effective Sample Size (Red, Left Axis)\nEstimated Relative Frequency\nof Cluster (Blue, Right Axis)") +
        xlab("Margin") +
        theme(legend.position = "none", axis.title = element_text(size = 10))
    
        return(g)
  }

  # xi -----
  if ("xi" == param) {
    ggeffssz_xi <- 
      effssz$xi %>% 
      filter(!.data$triv) %>%
      mutate(k = factor(.data$k)) %>%
      mutate(effssz = ifelse(.data$effssz <= nsave, .data$effssz, nsave))
    
    g <-
      ggplot(ggeffssz_xi) +
      geom_bar(
        mapping =
          aes(
            x = .data$p,
            y = .data$effssz
          ),
        fill = rgb(1, 0, 0, 0.6),
        stat = "identity"
      ) +
      geom_bar(
        data = local_freq_t,
        mapping = aes(
          x = .data$x,
          y = .data$frq_t * nsave
        ),
        stat = "identity",
        width = P - 0.1,
        fill = rgb(0, 0, 1, 0.6)
      ) +
      scale_y_continuous(limits = c(0, nsave), sec.axis = ~ . / nsave) +
      facet_grid(.data$j ~ .data$k, labeller = labeller(j = j_names, k = k_names)) +
      ylab("Effective Sample Size (Red, Left Axis)\nEstimated Relative Frequency\nof Cluster (Blue, Right Axis)") +
      xlab("Margin") +
      theme(legend.position = "none", axis.title = element_text(size = 10))
    
    return(g)
  }
  
  # psi -----
  if ("psi" == param) {
    ggeffssz_psi <- 
      effssz$psi %>% 
      filter(!.data$triv) %>%
      mutate(k = factor(.data$k)) %>%
      mutate(effssz = ifelse(.data$effssz <= nsave, .data$effssz, nsave))
    
    g <-
      ggplot(ggeffssz_psi) +
      geom_bar(
        mapping =
          aes(
            x = .data$p,
            y = .data$effssz
          ),
        fill = rgb(1, 0, 0, 0.6),
        stat = "identity"
      ) +
      geom_bar(
        data = glob_freq_t,
        mapping = aes(
          x = .data$x,
          y = .data$frq_t * nsave
        ),
        stat = "identity",
        width = P - 0.1,
        fill = rgb(0, 0, 1, 0.6)
      ) +
      scale_y_continuous(limits = c(0, nsave), sec.axis = ~ . / nsave) +
      facet_wrap(~ .data$k, labeller = labeller(k = k_names)) +
      ylab("Effective Sample Size (Red, Left Axis)\nEstimated Relative Frequency\nof Cluster (Blue, Right Axis)") +
      xlab("Margin") +
      theme(legend.position = "none", axis.title = element_text(size = 10))
    
    return(g)
  }

  # G -----  
  if ("G" == param) {
    glob_freq_t_matrix <- glob_freq_t
    glob_freq_t_matrix$x <- (((P + 1) * P / 2) + 1) / 2

    ggeffssz_G <- 
      effssz$G %>% 
      filter(!.data$triv) %>%
      mutate(k = factor(.data$k)) %>%
      mutate(pp = paste0(.data$p1, ",", .data$p2)) %>%
      mutate(
        Diagonal = 
          factor(
            .data$p1 == .data$p2,
            levels = c(TRUE, FALSE),
            labels = c("Diagonal", "Off-Diagonal")
            )
        ) %>%
      mutate(effssz = ifelse(.data$effssz <= nsave, .data$effssz, nsave))
    
    g <-
      ggplot(ggeffssz_G) +
      geom_bar(
        mapping =
          aes(
            x = reorder(.data$pp, as.integer(.data$Diagonal)),
            y = .data$effssz,
            fill = .data$Diagonal
          ),
        stat = "identity"
      ) +
      scale_fill_manual(values = c(rgb(1, 0, 0, 0.6), rgb(0.7, 0, 0.7, 0.6))) +
      geom_bar(
        data = glob_freq_t_matrix,
        mapping = aes(
          x = .data$x,
          y = .data$frq_t * nsave
        ),
        stat = "identity",
        width = (P + 1) * P / 2 - 0.1,
        fill = rgb(0, 0, 1, 0.6)
      ) +
      scale_y_continuous(limits = c(0, nsave), sec.axis = ~ . / nsave) +
      facet_wrap(~ .data$k, labeller = labeller(k = k_names)) +
      ylab("Effective Sample Size (Left Axis)\nEstimated Relative Frequency\nof Cluster (Blue, Right Axis)") +
      xlab("Margins") +
      theme(axis.text.x = element_text(angle = -90, vjust = .5, hjust = 1)) +
      theme(
        legend.title = element_blank(),
        legend.text = element_text(size = 6)
        )    
    return(g)
  }

  # E -----  
  if ("E" == param) {
    glob_freq_t_matrix <- glob_freq_t
    glob_freq_t_matrix$x <- (((P + 1) * P / 2) + 1) / 2

    ggeffssz_E <- 
      effssz$E %>% 
      filter(!.data$triv) %>%
      mutate(k = factor(.data$k)) %>%
      mutate(pp = paste0(.data$p1, ",", .data$p2)) %>%
      mutate(
        Diagonal = 
          factor(
            .data$p1 == .data$p2,
            levels = c(TRUE, FALSE),
            labels = c("Diagonal", "Off-Diagonal")
          )
      ) %>%
      mutate(effssz = ifelse(.data$effssz <= nsave, .data$effssz, nsave))
    
    g <-
      ggplot(ggeffssz_E) +
      geom_bar(
        mapping =
          aes(
            x = reorder(.data$pp, as.integer(.data$Diagonal)),
            y = .data$effssz,
            fill = .data$Diagonal
          ),
        stat = "identity"
      ) +
      scale_fill_manual(values = c(rgb(1, 0, 0, 0.6), rgb(0.7, 0, 0.7, 0.6))) +
      geom_bar(
        data = glob_freq_t_matrix,
        mapping = aes(
          x = .data$x,
          y = .data$frq_t * nsave
        ),
        stat = "identity",
        width = (P + 1) * P / 2 - 0.1,
        fill = rgb(0, 0, 1, 0.6)
      ) +
      scale_y_continuous(limits = c(0, nsave), sec.axis = ~ . / nsave) +
      facet_wrap(~ .data$k, labeller = labeller(k = k_names)) +
      ylab("Effective Sample Size (Left Axis)\nEstimated Relative Frequency\nof Cluster (Blue, Right Axis)") +
      xlab("Margins") +
      theme(axis.text.x = element_text(angle = -90, vjust = .5, hjust = 1)) +
      theme(
        legend.title = element_blank(),
        legend.text = element_text(size = 6)
      )    
    return(g)
  }
  
  return(NULL)
}
