# Function for cross validation of the GAM models. Based on the gamclass::CVgam but adapted as that function
# doesn't allow specification of different response distributions

gam.crossvalidation <- function(mod, method = "GCV.Cp", gamma = 1.2, nfold = 10, seed = 1503){
  
  data <- mod$model %>%
    mutate(population = as.factor(population))
  formula_mod <- formula(mod)
  family_mod <- mod$family$family
  
  set.seed(seed)
  
  cvparts <- sample(1:nfold, nrow(data), replace = TRUE)
  folds <- unique(cvparts)
  
  pred_df <- matrix(0, nrow = nrow(data), ncol = nfold)
  predicted <- numeric(nrow(data))

  rmse_df <- data.frame(
    nfold = c(1:nfold),
    RMSE_train = rep(0,nfold),
    RMSE_test = rep(0,nfold)
  )
  
  for(i in folds) {
    trainrows <- cvparts != i
    testrows <- cvparts == i
    
    # Extract the original data used to fit the model
    data_train <- data[trainrows, ]
    
    possibleError <- tryCatch(
    # Fit the updated model with original family, method, and gamma
    traingam <- gam(formula = formula_mod,
                    data = data[trainrows,],
                    family = family_mod,
                    method = method,
                    gamma = gamma),
    error=function(e) e)
    
    if(inherits(possibleError, "error")) next
    
    # Predict on full data (so I can retrieve RMSE from both train and test dataset)
    pred_df[,i] <- predict(object = traingam, newdata = data, select = T, type = "response")
    
    # Get RMSE for train and test data
    rmse_df$RMSE_train[i] <- caret::RMSE(pred_df[trainrows,i], data[trainrows,1])
    rmse_df$RMSE_test[i] <- caret::RMSE(pred_df[testrows,i], data[testrows,1])
    
    
    # Save predicted test data to seperate dataframe for overal RMSE
    predicted[testrows] <- pred_df[testrows,i]
    
  }
  
  
  
  df <- data.frame(
    RMSE = caret::RMSE(predicted, data[,1]),
    RMSE_train = mean(rmse_df$RMSE_train),
    RMSE_tests = mean(rmse_df$RMSE_test),
    train = I(list(rmse_df$RMSE_train)),
    test = I(list(rmse_df$RMSE_test))
    )
  
  return(df)
  
  
}
