options(shiny.port=7702)
library(clustRviz)
data("presidential_speech")
Xdat <- presidential_speech$X
carp.fit <- CARP(Xdat,obs.labels = presidential_speech$labels)
plot(carp.fit,type='interactive')
