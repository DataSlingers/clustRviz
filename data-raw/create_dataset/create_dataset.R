library(tm)
library(stringr)
library(SnowballC)
library(parallel)
library(Matrix)
library(tidyverse)
(f <- content_transformer(function(x, pattern) gsub(pattern, "", x)))
(f1 <- content_transformer(function(x, pattern) gsub(pattern, " ", x)))
inaug <- VCorpus(DirSource("../scrap_results/inaug"))
inaug <- tm_map(inaug,f,"<[^>]*>")
inaug <- tm_map(inaug,f1,"[^a-zA-Z]")
inaug <- tm_map(inaug,stripWhitespace)
inaug <- tm_map(inaug,content_transformer(tolower))
inaug <- tm_map(inaug,removeWords,stopwords("english"))
inaug <- tm_map(inaug,stemDocument)
dtm.inaug <- DocumentTermMatrix(inaug)

sou <- VCorpus(DirSource("../scrap_results/sou"))
sou <- tm_map(sou,f,"<[^>]*>")
sou <- tm_map(sou,f1,"[^a-zA-Z]")
sou <- tm_map(sou,stripWhitespace)
sou <- tm_map(sou,content_transformer(tolower))
sou <- tm_map(sou,removeWords,stopwords("english"))
sou <- tm_map(sou,stemDocument)
dtm.sou <- DocumentTermMatrix(sou)

conv <- VCorpus(DirSource("../scrap_results/conv"))
conv <- tm_map(conv,f,"<[^>]*>")
conv <- tm_map(conv,f1,"[^a-zA-Z]")
conv <- tm_map(conv,stripWhitespace)
conv <- tm_map(conv,content_transformer(tolower))
conv <- tm_map(conv,removeWords,stopwords("english"))
conv <- tm_map(conv,stemDocument)
dtm.conv <- DocumentTermMatrix(conv)

dtm <- c(dtm.sou,dtm.inaug,dtm.conv)
### person level
as.data.frame(as.matrix(dtm)) %>%
  tbl_df() %>%
  mutate(
    speech = rownames(dtm),
    name = str_replace(speech,":.*$","")
  ) %>%
  select(-speech) %>%
  group_by(name) %>%
  summarise_all(funs(sum)) -> dtm.df
dtm.mat <- dtm.df %>% select(-name) %>% as.matrix()
dtm.mat.log <- log(dtm.mat + 1)
wrd.var <- apply(dtm.mat.log,2,var)
top.wrd.var <- names(sort(wrd.var,decreasing = TRUE)[1:75])
dtm.mat.log <- dtm.mat.log[,colnames(dtm.mat.log) %in% top.wrd.var]
saveRDS(dtm.mat.log,"presidential_speech.rds")