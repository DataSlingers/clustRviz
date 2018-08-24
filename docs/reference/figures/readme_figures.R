library(clustRviz)
carp.fit <- CARP(
  presidential_speech
)
saveviz(
  carp.fit,
  file.name = 'path_dyn.gif',
  plot.type = 'path',
  image.type = 'dynamic'
)

saveviz(
  carp.fit,
  file.name = 'dend_dyn.gif',
  plot.type = 'dendrogram',
  image.type = 'dynamic',
  dend.labels.cex = 1.5,
  dend.branch.width = 2
)

saveviz.CARP
