set.seed(100)
require(graphics)

# Random Number Generation
X <- rTPSC(n = 1e5,w = 0.7,theta = -1,sigma = 3,delta = 5)

# Plot the histogram
hist(X, breaks = 100, freq = FALSE)

# The red dashed line should match the underlining histogram
points(x = seq(-70,50,length.out = 1000),
       y = dTPSC(x = seq(-70,50,length.out = 1000),
                 w = 0.7,theta = -1,sigma = 3,delta = 5),
       type = "l",
       col = "red",
       lwd = 3,
       lty = 2)
