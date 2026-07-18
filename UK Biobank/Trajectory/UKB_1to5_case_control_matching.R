
library(data.table)
library(dplyr)
library(MatchIt)
setwd("D:/ZYD/R/meta/biobank/time")
data <- fread("Covariates2.csv")
data$Sex <- as.factor(data$Sex)
data$Smoke <- as.factor(data$Smoke)
cols_to_match <- c("target_y", "Age", "Sex", "BMI", "TDI", "Smoke")
data <- na.omit(data, cols = cols_to_match)

match_control_bladder <- function(data) {
  set.seed(123) 
  m.out <- matchit(target_y ~ Age + BMI + TDI + Smoke, 
                   data = data,
                   method = "nearest", 
                   distance = "glm", 
                   ratio = 5,            
                   link = "logit",
                   exact = ~ Sex,         
                   caliper = c(.2, Age = 3, BMI = 5),
                   std.caliper = c(TRUE, FALSE, FALSE))
  
  m.data <- match.data(m.out)
  return(m.data)
}
matched_result <- match_control_bladder(data)
print(table(matched_result$target_y))
fwrite(matched_result, "Bladder_Cancer_Matched_1to5.csv", row.names = F)
