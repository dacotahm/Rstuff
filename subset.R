library(sqldf)

db <- dbConnect(SQLite(), './cuffData.db')
dbListTables(db)
expData <- (dbReadTable(db, 'geneExpDiffData'))

geneDiff <- read.csv(file = "gene_exp.diff", header = TRUE, sep = "\t")
drop <- c("gene", "test_stat", "status", "q_value", "test_id")
geneDiff <- geneDiff[ , !(names(geneDiff) %in% drop)]
geneDiff <- geneDiff[,0] = sub(pattern = ":.*", replacement = "", geneDiff$locus)


m <- subset(expData, 
              expData$significant == 'yes' & 
              (expData$sample_1 == "OLNA" | 
              expData$sample_2 == "OLNA") &
              expData$log2_fold_change > 5 & 
              expData$log2_fold_change != "Inf" &
              expData$log2_fold_change != "-Inf")

#pull genes from specific comparision 
m <- subset(geneDiff, 
            geneDiff$significant == "yes" &
            geneDiff$sample_1 == "OLNA" &
            geneDiff$sample_2 == "OLA215" &
            abs(geneDiff$log2.fold_change.) > 3 &
            geneDiff$log2.fold_change. != "Inf" &
            geneDiff$log2.fold_change. != "-Inf")

m <- m[order(m$log2.fold_change., decreasing = TRUE), ]

plot(log(m$value_1), log(m$value_2))
write.csv(m, file = "OLNAxOLA215.csv")


test[,0] = sub(pattern = ":.*", replacement = "", test$locus)

n <- subset(m, m$sample_1 == "OLNA" & m$sample_2 == "OLA215")
