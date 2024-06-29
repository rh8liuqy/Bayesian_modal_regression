set.seed(100)
require(graphics)

# Random Number Generation
X <- rDTP(n = 1e5,theta = 5,sigma1 = 7,sigma2 = 3,delta1 = 5,delta2 = 6)

# Plot the histogram
hist(X, breaks = 100, freq = FALSE)

# The red dashed line should match the underlining histogram
points(x = seq(-100,40,length.out = 1000),
       y = dDTP(x = seq(-100,40,length.out = 1000),
                theta = 5,sigma1 = 7,sigma2 = 3,delta1 = 5,delta2 = 6),
       type = "l",
       col = "red",
       lwd = 3,
       lty = 2)
