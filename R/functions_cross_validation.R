# Function for cross validation of the GAM models. Based on the gamclass::CVgam but adapted as that function
# doesn't allow specification of different response distributions

gam.crossvalidation <- function(mod, method = "GCV.Cp", gamma = 1.2, nfold = 10, seed = 1503){
  
  data <- mod$model
  formula_mod <- formula(mod)
  family_mod <- mod$family$family
  
  set.seed(seed)
  
  cvparts <- sample(1:nfold, nrow(data), replace = TRUE)
  
  folds <- unique(cvparts)
  predicted <- numeric(nrow(data))
  
  for(i in folds) {
    trainrows <- cvparts != i
    testrows <- cvparts == i
    
    # Extract the original data used to fit the model
    data_train <- data[trainrows, ]
    
    # Fit the updated model with original family, method, and gamma
    traingam <- gam(formula = formula_mod,
                    data = data[trainrows,],
                    family = family_mod,
                    method = method,
                    gamma = gamma)
    
    # Predict on the test set
    predicted[testrows] <- predict(object = traingam, newdata = data[testrows,], select = T, type = "response")
  }
  
  df <- data.frame(
    R2 = caret::R2(predicted, data[,1]),
    RMSE = caret::RMSE(predicted, data[,1]),
    MAE = caret::MAE(predicted, data[,1]))
  
  return(df)
  
  
}
