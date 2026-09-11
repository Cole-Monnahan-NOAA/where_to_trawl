#process_simulation_results
# very early attempt at processing simulation results

library(dplyr)
library(ggplot2)
library(stringr)

rmspe <- function(actual, predicted, na.rm = TRUE) {
  sqrt(mean(((actual - predicted) / actual)^2, na.rm = na.rm)) * 100
}

# let's summarize
df<-result %>%
  filter(str_starts(category, "pollock")) %>%
  # Group by replicate, method, and class to keep faceting structure intact
  group_by(method, category) %>%
  summarize(A_survey=sum(A_survey), A_true=sum(A_true))%>%
  group_by(method)%>% summarise(rmse =  rmspe(sum(A_true), sum(A_survey)))


result%>%ggplot(aes(x = method, y = rmspe(A_true,A_survey), fill = method)) +
  geom_violin(trim = FALSE, alpha = 0.7) +
  geom_boxplot(width = 0.1, color = "black", outlier.shape = NA) +
  facet_wrap(~ category) +
  labs(
    title = "Summed Percent Error for Pollock Categories by Replicate",
    x = "Method",
    y = "Summed Percent Error",
    fill = "Method"
  ) +
  theme_minimal()



# summarize the results
# make a violin plot
df<-result %>%
  filter(str_starts(category, "pollock")) %>%
  # Group by replicate, method, and class to keep faceting structure intact
  group_by(method) %>%
  summarise(discrep =  (abs(sum(A_survey)-sum(A_true))/sum(A_true)) )


## START with this #########################################

### trying a different way
df_by_size<-result %>%
  filter(str_starts(category, "pollock")) %>%
  group_by(method, iter,category) %>%  # add category here if desired
  summarize(A_survey=sum(A_survey), A_true=sum(A_true), 
            rmspe=rmspe(A_true,A_survey))%>%
  group_by(method,category)%>%
  summarise(mean_rmspe = mean(rmspe),
    sd_rmspe   = sd(rmspe),
    se_rmspe   = sd(rmspe) / sqrt(n()),
    ci_lower   = mean_rmspe - 1.96 * se_rmspe,
    ci_upper   = mean_rmspe + 1.96 * se_rmspe
  )%>%
  arrange(category, method)

df_by_all<-result %>%
  filter(str_starts(category, "pollock")) %>%
  group_by(method, iter) %>%  # add category here if desired
  summarize(A_survey=sum(A_survey), A_true=sum(A_true), 
            rmspe=rmspe(A_true,A_survey))%>%
  group_by(method)%>%
  summarise(mean_rmspe = mean(rmspe),
            sd_rmspe   = sd(rmspe),
            se_rmspe   = sd(rmspe) / sqrt(n()),
            ci_lower   = mean_rmspe - 1.96 * se_rmspe,
            ci_upper   = mean_rmspe + 1.96 * se_rmspe
  ) 