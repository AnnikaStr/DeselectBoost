# DeselectBoost
Function to enhance the sparsity of statistical boosting models by deselecting the base-learners with a minor impact on the total risk reduction. The deselection is done by considering the attributable (default) or cumulative risk reduction. 

# Example
```r
data("bodyfat", package = "TH.data")

### linear model 
bodyfat.glm <- glmboost(DEXfat ~ ., data = bodyfat)
cvr <- cvrisk(bodyfat.glm, grid = 1:500)
mstop(bodyfat.glm) <- mstop(cvr)

# deselection via attributable risk reduction (1% of the total risk reduction)
bodyfat.db <- DeselectBoost(bodyfat.gb, fam = Gaussian())  
coef(bodyfat.db)

# threshold value of 0.1 (10% of the total risk reduction)
bodyfat.db <- DeselectBoost(bodyfat.glm, fam = Gaussian(), tau = 0.1)    
coef(bodyfat.db)

# deselection via cumulative risk reduction
DeselectBoost(bodyfat.gb, fam = Gaussian(), method = 'cumulative')

### additive model
bodyfat.gam <- gamboost(DEXfat ~ ., data = bodyfat)
coef(bodyfat.gam)

bodyfat.db1 <- DeselectBoost(bodyfat.gam, fam = Gaussian(), data = bodyfat)


```
