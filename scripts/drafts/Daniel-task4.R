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

# LASSO estimates
fit <- glmnet(x_train, y_train, family = 'binomial')
fit_df <- tidy(fit)




# Fit best lambda
best_lambdas_chart <- cv_out_df %>% 
  arrange(lambda) %>% 
  filter(nzero<10)
best_lambda <- min(cv_out_df)
final_lasso <- glmnet(x_train,
                      y_train,
                      family = 'binomial',
                      alpha=1,
                      lambda = best_lambda)
lasso_coefs <- coef(final_lasso)
lasso_coefs_df <- as.data.frame(as.matrix(lasso_coefs))






