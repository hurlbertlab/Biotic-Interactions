library(ape)
# http://birdtree.org/subsets/
# ericson all species tree
yaml <- read_yaml("data/config.yaml", fileEncoding = "UTF-8")
trees <- read.nexus("data/output.nex")

AOU2 = AOU[,c("AOU_OUT", "CRC_SCI_NAME")] %>%
  distinct()
envloc$FocalAOU[envloc$FocalAOU == 4812] = 4810

envloc2 <- envloc %>%
  left_join(AOU2, by = c("FocalAOU" = "AOU_OUT")) %>%
  select(COMPSC, CRC_SCI_NAME) %>%
  distinct()
envloc2$CRC_SCI_NAME <- as.character(envloc2$CRC_SCI_NAME)
envloc2$V <- "Y"
no_match_tree$V1 <- as.character(no_match_tree$V1)
no_match_tree$V1 <- gsub("_", " ", no_match_tree$V1)
no_match_tree$V2 <- "Y"

test <- full_join(envloc2[,c("CRC_SCI_NAME", "V")], no_match_tree, by = c("CRC_SCI_NAME" = "V1"))
# test <- filter(no_match_tree$V1 %in% c(envloc2$CRC_SCI_NAME))
test <- envloc[is.na(match(envloc2$CRC_SCI_NAME,no_match_tree$V1)),]

prunedphy <- prune.sample(envloc2, trees)
kTest(envloc2, trees, 100)


phylosignal(envloc2, prunedphy[[1]], checkdata=TRUE)
