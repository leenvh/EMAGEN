
##################################

###  co-occurrence analysis ######

##################################


library(readxl)
library(dplyr)
library(tidyr)
library(purrr)
library(writexl)

#### crt vs k13 

df <- read_excel("data/Pfcrt_vs_Pfk13.xlsx", sheet = "Pfcrt_vs_Pfk13")

mut_cols <- setdiff(names(df), 
                    c("District", "Region", "Sample_ID", 
                      "API)", 
                      "API_strata", "Geography_cluster"))

# Function: co-occurrence analysis for a pair
cooccurrence_with_table <- function(data, var1, var2) {
  tab <- table(data[[var1]], data[[var2]])
  N <- sum(tab)
  
  chisq_val <- NA; df_val <- NA; chisq_p <- NA; fisher_p <- NA
  
  if (all(dim(tab) == c(2, 2))) {
    chisq_test <- suppressWarnings(chisq.test(tab, correct = FALSE))
    fisher_test <- fisher.test(tab)
    
    chisq_val <- unname(chisq_test$statistic)
    df_val <- chisq_test$parameter
    chisq_p <- chisq_test$p.value
    fisher_p <- fisher_test$p.value
  }
  
  tibble(
    Mutation1 = var1,
    Mutation2 = var2,
    ChiSq = round(chisq_val, 3),
    df = df_val,
    ChiSq_p = round(chisq_p, 4),
    Fisher_p = round(fisher_p, 4),
    Significant = ifelse(!is.na(chisq_p) & chisq_p < 0.05, "Yes", "No"),
    N = N,
    ContingencyTable = list(as.data.frame.matrix(tab))
  )
}

# Wrapper to run by group
analyze_group <- function(data, group_var) {
  data %>%
    group_by(.data[[group_var]]) %>%
    group_modify(~ {
      combn(mut_cols, 2, function(x) {
        cooccurrence_with_table(.x, x[1], x[2])
      }, simplify = FALSE) %>% bind_rows()
    }, .keep = TRUE)
}


results_api <- analyze_group(df, "API_strata")
results_geo <- analyze_group(df, "Geography_cluster")

#  contingency tables in tidy form
tables_api <- results_api %>%
  select(API_strata, Mutation1, Mutation2, ContingencyTable) %>%
  unnest(ContingencyTable)

tables_geo <- results_geo %>%
  select(Geography_cluster, Mutation1, Mutation2, ContingencyTable) %>%
  unnest(ContingencyTable)

# Save 
write_xlsx(
  list(
    Results_API = results_api %>% select(-ContingencyTable),
    Tables_API = tables_api,
    Results_Geography = results_geo %>% select(-ContingencyTable),
    Tables_Geography = tables_geo
  ),
  "Pfcrt_Pfk13_Cooccurrence.xlsx"
)


###### for crt vs mdr1

df <- read_excel("data/Pfcrt_vs_Pfmdr1.xlsx", sheet = "Pfcrt_vs_Pfmdr1")


mut_cols <- setdiff(names(df), 
                    c("District", "Region", "Sample_ID", 
                      "API)", 
                      "API_strata", "Geography_cluster"))

# Function: co-occurrence analysis for a pair
cooccurrence_with_table <- function(data, var1, var2) {
  tab <- table(data[[var1]], data[[var2]])
  N <- sum(tab)
  
  chisq_val <- NA; df_val <- NA; chisq_p <- NA; fisher_p <- NA
  
  if (all(dim(tab) == c(2, 2))) {
    chisq_test <- suppressWarnings(chisq.test(tab, correct = FALSE))
    fisher_test <- fisher.test(tab)
    
    chisq_val <- unname(chisq_test$statistic)
    df_val <- chisq_test$parameter
    chisq_p <- chisq_test$p.value
    fisher_p <- fisher_test$p.value
  }
  
  tibble(
    Mutation1 = var1,
    Mutation2 = var2,
    ChiSq = round(chisq_val, 3),
    df = df_val,
    ChiSq_p = round(chisq_p, 4),
    Fisher_p = round(fisher_p, 4),
    Significant = ifelse(!is.na(chisq_p) & chisq_p < 0.05, "Yes", "No"),
    N = N,
    ContingencyTable = list(as.data.frame.matrix(tab))
  )
}

# Wrapper to run by group
analyze_group <- function(data, group_var) {
  data %>%
    group_by(.data[[group_var]]) %>%
    group_modify(~ {
      combn(mut_cols, 2, function(x) {
        cooccurrence_with_table(.x, x[1], x[2])
      }, simplify = FALSE) %>% bind_rows()
    }, .keep = TRUE)
}

# 
results_api <- analyze_group(df, "API_strata")
results_geo <- analyze_group(df, "Geography_cluster")

# contingency tables in tidy form
tables_api <- results_api %>%
  select(API_strata, Mutation1, Mutation2, ContingencyTable) %>%
  unnest(ContingencyTable)

tables_geo <- results_geo %>%
  select(Geography_cluster, Mutation1, Mutation2, ContingencyTable) %>%
  unnest(ContingencyTable)

# Save 
write_xlsx(
  list(
    Results_API = results_api %>% select(-ContingencyTable),
    Tables_API = tables_api,
    Results_Geography = results_geo %>% select(-ContingencyTable),
    Tables_Geography = tables_geo
  ),
  "Pfcrt_Pfmdr1_Cooccurrence.xlsx"
)



###### for crt vs mdr1


df <- read_excel("data/Pfmdr1_vs_Pfk13.xlsx", sheet = "Pfmdr1_vs_Pfk13")


mut_cols <- setdiff(names(df), 
                    c("District", "Region", "Sample_ID", 
                      "API)", 
                      "API_strata", "Geography_cluster"))

# Function: co-occurrence analysis for a pair
cooccurrence_with_table <- function(data, var1, var2) {
  tab <- table(data[[var1]], data[[var2]])
  N <- sum(tab)
  
  chisq_val <- NA; df_val <- NA; chisq_p <- NA; fisher_p <- NA
  
  if (all(dim(tab) == c(2, 2))) {
    chisq_test <- suppressWarnings(chisq.test(tab, correct = FALSE))
    fisher_test <- fisher.test(tab)
    
    chisq_val <- unname(chisq_test$statistic)
    df_val <- chisq_test$parameter
    chisq_p <- chisq_test$p.value
    fisher_p <- fisher_test$p.value
  }
  
  tibble(
    Mutation1 = var1,
    Mutation2 = var2,
    ChiSq = round(chisq_val, 3),
    df = df_val,
    ChiSq_p = round(chisq_p, 4),
    Fisher_p = round(fisher_p, 4),
    Significant = ifelse(!is.na(chisq_p) & chisq_p < 0.05, "Yes", "No"),
    N = N,
    ContingencyTable = list(as.data.frame.matrix(tab))
  )
}

# Wrapper to run by group
analyze_group <- function(data, group_var) {
  data %>%
    group_by(.data[[group_var]]) %>%
    group_modify(~ {
      combn(mut_cols, 2, function(x) {
        cooccurrence_with_table(.x, x[1], x[2])
      }, simplify = FALSE) %>% bind_rows()
    }, .keep = TRUE)
}


results_api <- analyze_group(df, "API_strata")
results_geo <- analyze_group(df, "Geography_cluster")

#  contingency tables in tidy form
tables_api <- results_api %>%
  select(API_strata, Mutation1, Mutation2, ContingencyTable) %>%
  unnest(ContingencyTable)

tables_geo <- results_geo %>%
  select(Geography_cluster, Mutation1, Mutation2, ContingencyTable) %>%
  unnest(ContingencyTable)

# Save 
write_xlsx(
  list(
    Results_API = results_api %>% select(-ContingencyTable),
    Tables_API = tables_api,
    Results_Geography = results_geo %>% select(-ContingencyTable),
    Tables_Geography = tables_geo
  ),
  "Pfmdr1_Pfk13_Cooccurrence.xlsx"
)



