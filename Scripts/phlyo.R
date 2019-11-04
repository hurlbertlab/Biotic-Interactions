library(ape)
library(gdata)
# http://birdtree.org/subsets/
# ericson all species tree
yaml <- read_yaml("data/config.yaml", fileEncoding = "UTF-8")
trees <- read.nexus("data/output.nex")
envloc <- read.csv("data/envloc.csv", header =TRUE)
AOU2 = AOU[,c("AOU_OUT", "CRC_SCI_NAME")] %>%
  distinct()
envloc$FocalAOU[envloc$FocalAOU == 4812] = 4810

envloc2 <- envloc %>%
  left_join(AOU2, by = c("FocalAOU" = "AOU_OUT")) %>%
  filter(!FocalAOU %in% c(7221, 7350, 7360, 7380)) %>%
  mutate(match_names = gsub(" ", "_", CRC_SCI_NAME)) %>%
  select(match_names, COMPSC) %>%
  distinct()
envloc2$match_names <- as.character(envloc2$match_names)
# envloc2$match_names[envloc2$FocalAOU == 4810] = "Aphelocoma_californica"
# envloc2$match_names[envloc2$FocalAOU == 7350] = "Poecile_atricapillus"
# envloc2$match_names[envloc2$FocalAOU == 7360] = "Poecile_carolinensis"
# envloc2$match_names[envloc2$FocalAOU == 7380] = "Poecile_gambeli"     
# envloc2$V <- "Y"
# no_match_tree$V1 <- as.character(no_match_tree$V1)
# no_match_tree$V1 <- gsub("_", " ", no_match_tree$V1)
# no_match_tree$V2 <- "Y"
# test <- full_join(envloc2[,c("CRC_SCI_NAME", "V")], no_match_tree, by = c("CRC_SCI_NAME" = "V1"))
# # test <- filter(no_match_tree$V2 %in% c(envloc2$CRC_SCI_NAME))
# test <- envloc[is.na(match(envloc2$CRC_SCI_NAME,no_match_tree$V1)),]
envloc2$match_names <- reorder.factor(envloc2$match_names, new.order=tree1$tip.label)

envloc3 <- envloc2 %>%
  arrange(match_names) %>%
  remove_rownames %>% 
  column_to_rownames(var="match_names") 
envloc3 <- as.matrix(envloc3)

# prunedphy <- prune.sample(envloc3, trees)
# kTest(envloc2, trees, 100)
tree1 <- trees[[1]]


newdf <- c()
for(i in 1:100){
  tsub = trees[[i]]
  n = phylosignal(envloc3, tsub)
newdf = rbind(newdf, n)
}
newdf = data.frame(newdf)
mean(newdf$K)
