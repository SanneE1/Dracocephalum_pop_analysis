# Function for cross validation of the GAM models. Based on the gamclass::CVgam but adapted as that function
# doesn't allow specification of different response distributions

gam.selection.criteria <- function(mod, formula = NA, response_column, 
                                  do.cross.validation = F, subset_length = NA,
                                method = "GCV.Cp", gamma = 1.4, seed = 1503){
  if (do.cross.validation) {
  data <- mod$model %>%
    mutate(population = as.factor(population),
           year_t0 = as.integer(as.character(year_t0))) %>%
    as.data.frame
  if(is.na(formula)) {
    formula_mod <- formula(mod)
  } else {
    formula_mod <- as.formula(formula)
  }
  
  family_mod <- mod$family$family

  set.seed(seed)

  years <- unique(data$year_t0)

  if(!is.na(subset_length)) {
    years <- sample(years, subset_length)
  }

  pred_df <- matrix(0, nrow = nrow(data), ncol = length(years))
  predicted <- numeric(nrow(data))

  rmse_df <- data.frame(
    year_left_out = c(1:length(years)),
    RMSE_train = rep(0,length(years)),
    RMSE_test = rep(0,length(years))
  )

  for(i in c(1:length(years))) {
    trainrows <- data$year_t0 != years[i]
    testrows <- data$year_t0 == years[i]

    # Create training data
    data_train <- data[trainrows, ]

    # possibleError <- tryCatch(
    # Fit the updated model with original family, method, and gamma
    traingam <- gam(formula = formula_mod,
                    data = data_train,
                    family = family_mod,
                    method = method,
                    gamma = gamma)

    # Predict on full data (so I can retrieve RMSE from both train and test dataset)
    pred_df[,i] <- predict(object = traingam, newdata = data, type = "response")

    # Save predicted test data to separate dataframe for overall RMSE
    predicted[testrows] <- pred_df[testrows,i]

    not_na <- which(!is.na(pred_df[,i]))

    # Get RMSE for train and test data
    rmse_df$RMSE_train[i] <- caret::RMSE(pred_df[which(trainrows & not_na),i], data[which(trainrows & not_na), response_column])
    rmse_df$RMSE_test[i] <- caret::RMSE(pred_df[which(testrows & not_na),i], data[which(testrows & not_na), response_column])

    rmse_df$GCV[i] <- traingam$gcv.ubre
  }
  
  df <- data.frame(
    RMSE = caret::RMSE(predicted[!is.na(predicted)], data[!is.na(predicted),response_column]),
    RMSE_train = mean(rmse_df$RMSE_train),
    RMSE_tests = mean(rmse_df$RMSE_test, na.rm = T),
    GCV = mod$gcv.ubre,
    AICc = AICc(mod)
    )
  } else {
    df <- data.frame(
      GCV = mod$gcv.ubre,
      AICc = AICc(mod)
    )
}
  return(df)
  
  
}
