library(clustRviz)
carp.fit <- CARP(presidential_speech)
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
  image.type = 'dynamic'
)
