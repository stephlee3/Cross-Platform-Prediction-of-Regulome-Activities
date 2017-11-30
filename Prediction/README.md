# Prediction Result and Model Evaluation
We use C++ to implement this part ,since it is more computational efficienct than R code. 

- pred.cpp

This file is to calculate the prediction result.  
Input: selected predcitors, coefficient of the regression model  
Output: predicted DH value  

- corr.cpp

This file is to evaluate the prediction performance.  
Input:　predicted DH value, true DH value  
Output: cross-cell-type correlation, prediction squared error  
