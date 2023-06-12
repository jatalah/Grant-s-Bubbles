library(tidyverse)
library(ggpubr)
library(janitor)
library(sjPlot)
library(glmmTMB)
library(emmeans)
library(broom)

d <- 
  read_csv('data/Waikawa dataset.csv') %>% 
  clean_names() %>% 
  mutate(time_month =factor(time_months)) %>% 
  pivot_longer(cols = c(bare_light_biofilm:macrofouling)) %>% 
  mutate(name = fct_relevel(name, "bare_light_biofilm", "thick_biofilm_algal_fuzz"))


labs <- as_labeller(c("bare_light_biofilm" = "Bare space",
                      "thick_biofilm_algal_fuzz" = "Biofilm", 
                      "macrofouling"= "Macrofouling"))

# boxplot---------
ggplot(d, aes(time_month, value / 100, fill = treatment)) +
  geom_boxplot(alpha = .5) +
  theme_javier() +
  scale_y_continuous(labels = scales::percent_format()) +
  ggplot2::labs(x = "Months", fill = "Treatment", y = NULL) +
  facet_wrap(~ name, labeller = labs) +
  stat_compare_means(label = "p.signif",
                     size = 3,
                     method = "wilcox.test")

ggsave(
  last_plot(),
  filename = 'figures/bubbles_boxplot.png',
  width = 9,
  height = 3,
  dpi = 300
)
  
# convert percentage data for beta glms----
betareg_scaler <-
  function(y) {
    n <- length(!is.na(y))
    (y/100 * (n-1)  + 0.5) / n
  }

# run glms----
glms <-
  d %>%
  group_by(name) %>%
  mutate(value = betareg_scaler(value)) %>%
  nest() %>%
  mutate(
    glms = map(
      data,
      ~ glmmTMB(
        value ~ time_month * treatment + (1 | block),
        data = .x,
        family = beta_family()
      )
    ),
    pairwise = map(glms, ~ emmeans(.x, ~ time_month * treatment) %>% pairs() %>% tidy())
  )

# model summary---  
tab_model(
  glms$glms,
  dv.labels = c("Bare space", "Macrofouling", "Biofilm"),
  file = "tables/glmm_summaries.doc"
)


# pairwise test summary----
glms %>% 
  select(pairwise) %>% 
  unnest(cols = c(pairwise)) %>% 
  select(-term, -df, -null.value) %>% 
  dplyr::filter(
    contrast %in% c(
      "time_month0 Bubbled - time_month0 Control" ,
      "time_month1 Bubbled - time_month1 Control",
      "time_month2 Bubbled - time_month2 Control",
      "time_month3 Bubbled - time_month3 Control",
      "time_month4 Bubbled - time_month4 Control"
    )
  ) %>% 
  separate_wider_delim(contrast, names = c("A", "B"), " - ") %>% 
  write_csv('tables/pairwise_posthoc_tests.csv')
