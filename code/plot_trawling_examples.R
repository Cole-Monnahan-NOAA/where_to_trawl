#plot_trawling_examples
# let's plot trawling examples
# shows backscatter along a survey with trawl locations for various stratedgies
# idea is to visualize the basic results

library(ggplot2)
library(here)

# read in an sa vector from a simulated survey

# read in a 
df<-read_csv("data/test_survey_500_EDSUs.csv")
sa=df$sa #sa Values

n_hauls =5

# get trawl locations
ind_max=get_trawl_locs_max(sa, n_hauls, thresh_dist =  round(length(sa)/(n_hauls*4)))
ind_random=get_trawl_locs_random(sa, n_hauls, thresh_dist =  round(length(sa)/(n_hauls*4)))
ind_cumsum=get_trawl_locs_cum(sa, n_hauls)
ind_systematic=get_trawl_locs_systematic(sa, n_hauls)


# Basic plot overlaying both datasets
ggplot() +
  geom_point(data = survey, aes(x = x, y = sa, color = "Survey")) +
  geom_point(data = ind_max, aes(x = index, y = value, color = "max"),size=4,shape=0, stroke=2) +
  geom_point(data = ind_random, aes(x = index-.3, y = value, color = "random"),size=4,shape=1,stroke=2) +
  geom_point(data = ind_cumsum, aes(x = index-.2, y = value, color = "cumsum"),size=4,shape=5,stroke=2) +
  geom_point(data = ind_systematic, aes(x = index+.3, y = value, color = "systematic"),size=4,shape=6,stroke=2) +
  scale_color_manual(
    name = "Trawl placement",
    values = c("max" = "red","random"="blue", "cumsum"="green", "systematic"="orange" )) +
  labs(x = "EDSU", y = "Backscatter (sa, m^2/m^2)") +
  theme_minimal()
 

ggsave(here("results/plot_trawling_example.png"), width = 10, height = 6, dpi = 150)
