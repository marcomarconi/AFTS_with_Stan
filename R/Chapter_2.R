source('common.R')

## @knitr init

ggplot()


## @knitr show_AR_model

model_file <- '../models/AR.stan'
cat(paste(readLines(model_file)), sep = '\n')

## @knitr plot_gnp

gnp <- scan(file="../data/dgnp82.txt")*100
gnp1 <- ts(gnp,frequency=4,start=c(1947,2))
plot(gnp1)
points(gnp1)
