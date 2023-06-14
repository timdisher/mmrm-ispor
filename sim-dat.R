# Next steps:
#  - Allow for a user specified correlation matrix for within study
#  - Plot of observed vs predicted
library(MASS)
library(R2WinBUGS)
library(dplyr)
simulate_meta_data <- function(num_trials,
                               mus = c(1, 1.5, 2, 3, 4, 4.5, 5, 7),
                               se_sd = 1) {
 num_timepoints <- length(mus)
 df <- num_timepoints # Degrees of freedom for Wishart distribution
 
 decay <- 0.5 # Decay parameter for AR(1) structure (between 0 and 1)
 
 # Create an empty 8x8 matrix
 scale_matrix <- matrix(0, nrow = num_timepoints, ncol = num_timepoints)
 
 # Fill the scale matrix with values according to the AR(1) structure
 for (i in 1:num_timepoints) {
   for (j in 1:num_timepoints) {
     scale_matrix[i, j] <- decay^abs(i - j)
   }
 }
  
  # Set up storage for output
  sigma_sample <- rWishart(1, df, scale_matrix)
  cor <- cov2cor(sigma_sample[,,1])
  ses <- mus/1.7
  

 sim.trials <- purrr::map(1:num_trials, ~ {
    
    ses <- abs(MASS::mvrnorm(1, mu = ses, Sigma = diag(rep(se_sd, length(mus)))))
    cov <- diag(ses) %*% cor %*% diag(ses)
    ys <- mvrnorm(1, mu = mus, Sigma = cov)
  
    list(y = ys, se = ses)  
  })
  
 y <- purrr::map(sim.trials, "y") %>%
   do.call(rbind,.)
  
 se <- purrr::map(sim.trials, "se") %>%
   do.call(rbind,.)
 # Generate data for each trial

  
  # Return list of data for all trials
  return(list(y = y, se = se, cor = cor))
}


set.seed(950)
true_mean <- c(1, 1.5, 2, 3, 4, 4.5, 5, 7)
data <- simulate_meta_data(5, true_mean)
y.long <- as.data.frame(data$y) %>%
  dplyr::mutate(trial = 1:dplyr::n()) %>%
  tidyr::pivot_longer(cols = -trial, values_to ="y")
se.long <- as.data.frame(data$se) %>%
  dplyr::mutate(trial = 1:dplyr::n()) %>%
  tidyr::pivot_longer(cols = -trial, values_to ="se")


pass <- FALSE
while(pass == FALSE){
  
(data.in <- dplyr::full_join(y.long, se.long) %>%
  dplyr::group_by(trial) %>%
  dplyr::mutate(dur = sample(2:8,1),
                name = 1:dplyr::n(),
                y = ifelse(name > dur, NA, y),
                mean.se.trial = mean(se, na.rm = TRUE),
                se = ifelse(name > dur, mean.se.trial, se)
                #se = ifelse(name > dur, NA, se)
               ) %>%
  dplyr::ungroup() %>%
  dplyr::select(-c(dur, mean.se.trial)) %>%
  tidyr::pivot_wider(names_from = name, values_from = c(y, se)) %>%
  dplyr::select(-trial) %>%
  as.matrix())

if(sum(!is.na(data.in[,8])) > 1){
  pass <- TRUE
} 
  }
data.in
num_trials <- nrow(data.in)
num_timepoints <- ncol(data.in)/2

model_data <- list(num_trials = num_trials,
                   num_timepoints = num_timepoints,
                   y = data.in[,1:num_timepoints],
                   sd = data.in[,(num_timepoints+1):(num_timepoints*2)],
                   R = data$cor
                   )

# MCMC settings
n_chains <- 3
n_burnin <- 20000
n_iter <- n_burnin + 20000

n_thin <- 1
parameters <- c("mu", "y")
# Run WinBUGS using R2WinBUGS
bugs_fit <- R2WinBUGS::bugs(model.file = "model.txt", data = model_data, 
                 inits = NULL, 
                 n.chains = n_chains,
                 n.iter = n_iter, 
                 n.burnin = n_burnin, 
                 n.thin = n_thin,
                 DIC = FALSE, 
                 parameters = parameters, 
                 debug = FALSE,
                 bugs.directory = Sys.getenv("bugs"))


Rnocor <- diag(1, num_timepoints)

model_data$R <- Rnocor
bugs_fit2 <- R2WinBUGS::bugs(model.file = "model.txt", 
                            data = model_data, 
                            inits = NULL, 
                            n.chains = n_chains,
                            n.iter = n_iter, 
                            n.burnin = n_burnin, 
                            n.thin = n_thin,
                            DIC = FALSE, 
                            parameters = parameters, 
                            debug = FALSE,
                            bugs.directory = Sys.getenv("bugs"))


# Summarize posterior results
imp.means <- bugs_fit$summary %>%
  tibble::as_tibble(rownames = "param") %>%
  dplyr::select(
    param,
    sd,
    lwr = `2.5%`,
    point = `50%`,
    upr = `97.5%`
  ) %>%
  dplyr::mutate(
    trial = as.numeric(dplyr::case_when(grepl("mu", param) ~ "10",
                             TRUE ~ stringr::str_extract(param, "\\d"))),
    timepoint = stringr::str_extract(param, "\\d+(?=[^\\d]*$)"),
    param = stringr::str_extract(param, "[[:alpha:]]+"),
    type = ifelse(grepl("y", param), "Imputed", "Mean")) 
  

imp.means.nocor <- bugs_fit2$summary %>%
  tibble::as_tibble(rownames = "param") %>%
  dplyr::select(
    param,
    sd,
    lwr = `2.5%`,
    point = `50%`,
    upr = `97.5%`
  ) %>%
  dplyr::mutate(
    trial = as.numeric(dplyr::case_when(grepl("mu", param) ~ "10",
                                        TRUE ~ stringr::str_extract(param, "\\d"))),
    timepoint = stringr::str_extract(param, "\\d+(?=[^\\d]*$)"),
    param = stringr::str_extract(param, "[[:alpha:]]+"),
    type = ifelse(grepl("y", param), "Imputed", "Mean")) 

true.dat <- data$y %>% as.data.frame() %>%
  purrr::set_names(as.character(1:8)) %>%
  dplyr::mutate(trial = 1:dplyr::n()) %>%
  tidyr::pivot_longer(cols = -trial, values_to = "point") %>%
  dplyr::rename(timepoint = name) %>%
  dplyr::mutate(param = "y")

estools::es_pal(sel = 3)

pdat <- 
  data.in %>%
  as.data.frame() %>%
  dplyr::mutate(trial = 1:dplyr::n()) %>%
  tidyr::pivot_longer(cols = -trial, values_to = "point") %>%
  tidyr::separate(name, into = c("param", "timepoint")) %>%
  dplyr::filter(param == "y") %>%
  dplyr::mutate(type = "Observed") %>%
  dplyr::bind_rows(imp.means %>% dplyr::filter(type != "Mean")) %>%
  dplyr::mutate(timepoint = as.numeric(timepoint)) %>%
  dplyr::left_join(., imp.means.nocor %>% dplyr::filter(type == "Mean") %>%
  dplyr::select(c(timepoint, nocor = point)) %>%
    dplyr::mutate(timepoint = as.numeric(timepoint))
) %>%
    dplyr::group_by(trial) %>%
    dplyr::mutate(na.num = cumsum(!is.na(point)),
                  nocor = dplyr::case_when(is.na(point) ~ nocor,
                                           dplyr::lead(na.num,1) == na.num ~ point,
                                           TRUE ~ NA))
pdat %>% 
  tidyr::drop_na(point) %>%
  ggplot2::ggplot(ggplot2::aes(x = timepoint, y = point)) +
  ggplot2::facet_grid(~trial,
                      labeller = ggplot2::labeller(trial= c(
                        "1" = "Trial 1",
                        "2" = "Trial 2",
                        "3" = "Trial 3",
                        "4" = "Trial 4",
                        "5" = "Trial 5"
                      ))) +
  ggplot2::geom_line(ggplot2::aes(group = trial, colour = type))+
  ggplot2::geom_point(ggplot2::aes(colour = type)) +
  ggplot2::geom_line(ggplot2::aes(x = timepoint, y = nocor, linetype = "Imputed (no correlation)"),
                     colour = "#007398",
                     data = pdat %>% tidyr::drop_na(nocor)) +
  ggplot2::geom_line(ggplot2::aes(x = as.numeric(timepoint), y = point, group = trial,
                                  linetype = "True"), data = true.dat, inherit.aes = FALSE) +
  ggplot2::geom_ribbon(ggplot2::aes(x = as.numeric(timepoint), 
                                    ymin = lwr, ymax = upr,group = trial),
                       data = imp.means %>%
                         dplyr::filter(param == "y", type == "Imputed"),
                       alpha = 0.2, fill = "#ED8B00") +
  estools::scale_colour_es(sel = c(2,1)) +
  estools::scale_fill_es(sel = 2) +
  ggplot2::labs(x = "Time", y = "Y",
                linetype = "",
                colour = "") +
  ggplot2::theme_minimal(base_size = 14) +
  ggplot2::theme(legend.position = "bottom")



sum.tab <- imp.means.nocor %>%
  dplyr::filter(param == "mu") %>%
  dplyr::select(
    param,
    timepoint,
    sd.nocor = sd,
    point.nocor = point
    
  ) %>% 
  dplyr::left_join(imp.means %>%
                      dplyr::select(param,timepoint, sd, point)) %>%
  dplyr::mutate_if(is.numeric, round, 2) %>%
  dplyr::mutate(true_mean = true_mean,
                mmrm = glue::glue("{point} ({sd})"),
                nocor = glue::glue("{point.nocor} ({sd.nocor})")
               )

write.table(sum.tab$nocor, file = "clipboard",sep = "\t", row.names = FALSE)
