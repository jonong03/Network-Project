# Network-Project
Motivated by analyzing whether correlation or coviarate effect exist in panel voting behaviour, we proposed a logistic quadratic exponential model to solve this multivariate correlated binary data problem. 
The solution is inspired and based on Ising model in the undirected graphical modelling context.

There are three functions in the R code: <br>
**rising**: random ising data generator <br>
**dising**: density function for ising data <br>
**MLE_ising**: estimate parameteres from data, methods include: <br>
1) "*Newton*": Maximizing joint probability function directly from fisher score information. 
2) "*Logit*": Logistic function for estimating parameters, with variance estimated from variance-covariance matrix of the data.
3) "Logit0": Logistic estimates from glm
4) "Logit1": Logistic estimates from glm, with variance estimated from fisher score information analytically.
5) "Poisson": Poisson regression estimates from glm
6) "gd": Gradient Descent
7) "gee": Poisson estimates from glm
8) "ridge": Ridge logistic regression, with variance estimated from vairance-covariance matrix of the data.

The recommended methods are "Newton" and "Logit" due to the proven unbiasedness, robustness and efficiencies. 
