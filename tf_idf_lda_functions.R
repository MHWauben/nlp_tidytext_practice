
### Traditional tf-idf on pre-existing categories ----
tf_idf_out <- function(data, grouping_var){
  grouping_var <- enquo(grouping_var)
  tf_idf_out <- data %>%
    count(!!grouping_var, word, sort = TRUE) %>%
    group_by(!!grouping_var) %>%
    mutate(total = sum(n)) %>%
    arrange(-n) %>%
    mutate(rank = row_number(),
           `term frequency` = n/total) %>%
    bind_tf_idf(word, !!grouping_var, n) %>%
    arrange(desc(tf_idf)) %>%
    ungroup() %>%
    mutate(word = factor(word, levels = rev(unique(word))))
  return(tf_idf_out)
}

vis_tf_idf <- function(data, grouping_var){
  grouping_var <- enquo(grouping_var)
  data %>%
    group_by(!!grouping_var) %>%
    top_n(10, tf_idf) %>%
    ungroup() %>%
    mutate(word = reorder(word, tf_idf)) %>%
    ggplot(aes(word, tf_idf, fill = !!grouping_var)) +
    geom_col(show.legend = FALSE) +
    labs(title = paste0("TF-IDF by ", as_label(grouping_var)),
         x = NULL, y = "tf-idf") +
    facet_wrap(vars(!!grouping_var), ncol = 2, scales = "free") +
    coord_flip()
}

### LDA functions ----

lda_dtm <- function(data, grouping_var, k = 2){
  grouping_var <- enquo(grouping_var)
  prep_data <- tf_idf_out(data, !!grouping_var)
  dtm_inc <- prep_data %>%
    cast_dtm(document = !!grouping_var, term = word, value = n)
  
  lda_out <- topicmodels::LDA(dtm_inc, k = k, control = list(seed = 1234))
  betas   <- tidy(lda_out, matrix = "beta")
  gammas  <- tidy(lda_out, matrix = "gamma")
  return(list(lda_out = lda_out,
              beta = betas,
              gamma = gammas))
}

lda_perinc <- function(data, grouping_var){
  grouping_var <- enquo(grouping_var)
  groups <- as_label(substitute(grouping_var))
  data_spread <- data %>%
    mutate(topic = paste0(groups, "_", topic)) %>%
    tidyr::spread(key = "topic", value = "gamma")
}

vis_lda <- function(data, n = 10){
  data %>%
    group_by(topic) %>%
    top_n(n, beta) %>%
    ungroup() %>%
    arrange(topic, -beta) %>%
    mutate(term = reorder(term, beta)) %>%
    ggplot(aes(term, beta, fill = factor(topic))) +
    geom_col(show.legend = FALSE) +
    facet_wrap(~ topic, scales = "free") +
    coord_flip()
}

comp_topics <- function(data, k = 4){
  comparison <- lapply(seq_len(ncol(combn(1:k, 2))), function(i) combn(1:k, 2)[,i])
  data_out <- data.frame()
  for(i in 1:length(comparison)) {
    first_comp <- data %>%
      dplyr::filter(topic %in% comparison[[i]]) %>%
      dplyr::mutate(topic = ifelse(topic == max(topic), 2, 1)) %>%
      dplyr::mutate(topic = paste0("topic", topic)) %>%
      tidyr::spread(topic, beta) %>%
      dplyr::filter(topic1 > .001 | topic2 > .001) %>%
      dplyr::mutate(log_ratio = log2(topic2 / topic1),
                    comparison = paste(comparison[[i]], collapse = " vs "))
    data_out <- rbind(data_out, first_comp)
  }
  return(data_out)
}

vis_betas <- function(data, n = 20){
  all_comps <- unique(data$comparison)
  if (length(all_comps) > 8) {
    all_comps <- all_comps[1:8]
  } else {
    all_comps <- all_comps
  }
  n <- ifelse(n > 20, 20, n)
  data %>%
    filter(comparison %in% all_comps) %>%
    group_by(comparison) %>%
    top_n(n, abs(log_ratio)) %>%
    ungroup() %>%
    mutate(term = factor(term, levels = unique(term[order(log_ratio)]))) %>%
    arrange(desc(log_ratio)) %>%
    ggplot(aes(y = log_ratio, x = term, fill = log_ratio > 0))+
    geom_bar(stat = "identity")+
    coord_flip()+
    facet_wrap(~comparison, scales = "free")+
    theme(legend.position = "none")
}

vis_gammas <- function(data, n){
  n <- ifelse(n > 20, 20, n)
  data %>%
    filter(document %in% sample(levels(as.factor(document)), n, replace = TRUE)) %>%
    ungroup() %>%
    ggplot(aes(x = document, y = gamma, fill = as.factor(topic)))+
    geom_bar(stat = "identity")+
    labs(title = "Topic make-up",
         x = " ", fill = "Topic")+
    coord_flip()
}

full_LDA <- function(data, grouping_var, k, n, comparison){
  grouping_var <- enquo(grouping_var)
  
  lda_out <- lda_dtm(data, !!grouping_var, k)
  comp_outputs <- comp_topics(lda_out$beta, k)
  perinc_out <- lda_perinc(lda_out$gamma, !!grouping_var)
  
  vis_lda <- vis_lda(lda_out$beta, n)
  beta_vis <- vis_betas(comp_outputs, n)
  gamma_vis <- vis_gammas(lda_out$gamma, n)
  
  return(list(lda = lda_out,
              comps = comp_outputs,
              per_inc = perinc_out,
              vis = list(full = vis_lda,
                         beta = beta_vis,
                         gamma = gamma_vis)))
}