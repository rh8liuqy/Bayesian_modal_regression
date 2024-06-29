rm(list = ls())
library(tidyverse)

## Wind speed NO2 PM10 http://lib.stat.cmu.edu/datasets/
## The data are a subsample of 500 observations from a data set that originate in a study where air pollution at a road is related to traffic volume and meteorological variables, collected by the Norwegian Public Roads Administration. The response variable (column 1) consist of hourly values of the logarithm of the concentration of NO2 (particles), measured at Alnabru in Oslo, Norway, between October 2001 and August 2003. The predictor variables (columns 2 to 8) are the logarithm of the number of cars per hour, temperature $2$ meter above ground (degree C), wind speed (meters/second), the temperature difference between $25$ and $2$ meters above ground (degree C), wind direction (degrees between 0 and 360), hour of day and day number from October 1. 2001. Submitted by Magne Aldrin (magne.aldrin@nr.no). [28/Jul/04] (19kbytes)
rm(list = ls())
df1 <- read.table("PM10.dat")
df1 <- df1[,c(1,4)]
colnames(df1) <- c("PM10","wind_speed")
df1$PM10 <- exp(df1$PM10)
df1 %>%
  ggplot(aes(x = wind_speed, y = PM10))  +
  geom_point()
write.csv(df1,"pm10.csv",row.names = FALSE)
