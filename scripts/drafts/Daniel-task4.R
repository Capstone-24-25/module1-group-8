library(tidyverse)
library(infer)
library(randomForest)
library(tidymodels)
library(modelr)
library(yardstick)
library(glmnet)
library(ggplot2)
load('data/biomarker-clean.RData')

## Lasso Regression
###################

# Create class category
biomarker <- biomarker_clean %>% 
  select(-ados) %>% 
  mutate(class = as.numeric(group == 'ASD'))

# Partition
set.seed(102824)
partitions <- biomarker %>% 
  initial_split(prop=0.8)

x_train <- training(partitions) %>% 
  select(-group, -class) %>% 
  as.matrix()

y_train <- training(partitions) %>% 
  pull(class)

# Lambda Selection
cv_out <- cv.glmnet(x_train,
                    y_train,
                    family = 'binomial',
                    alpha = 1)
cv_out_df <- tidy(cv_out)

# Fit best lambda
best_lambda <- cv_out_df %>% 
  arrange(lambda) %>% 
  filter(nzero<7) %>% 
  slice(1) %>% 
  pull(lambda)

final_lasso <- glmnet(x_train,
                      y_train,
                      family = 'binomial',
                      alpha=1,
                      lambda = best_lambda)
lasso_coefs <- coef(final_lasso)
lasso_coefs_df <- as.data.frame(as.matrix(lasso_coefs)) %>% 
  filter(s0!=0)

proteins_lasso <- lasso_coefs_df %>% 
  mutate(proteins = rownames(lasso_coefs_df)) %>% 
  filter(rownames(lasso_coefs_df) != "intercept") %>%
  slice(-1) %>% 
  pull(proteins)

## Logistic Regression
#######################

# Pull proteins
fit_train <- training(partitions) %>% 
  select(group, any_of(proteins_lasso)) %>% 
  mutate(class = (group == 'ASD')) %>% 
  select(-group)
  
fit_test <- testing(partitions) %>% 
  select(group, any_of(proteins_lasso)) %>% 
  mutate(class = (group == 'ASD')) %>% 
  select(-group)

# Fit model with training data
fit <- glm(class ~ .,
           data = fit_train,
           family = 'binomial')

# Evaluate errors on test set
class_metrics <- metric_set(sensitivity, 
                            specificity, 
                            accuracy,
                            roc_auc)

fit_test %>%
  add_predictions(fit, type = 'response') %>%
  mutate(est = as.factor(pred > 0.5), tr_c = as.factor(class)) %>%
  class_metrics(estimate = est,
                truth = tr_c, pred,
                event_level = 'second')
