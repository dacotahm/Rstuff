library(sqldf)

#This opens a connection to SQL database.  Make sure the file pathway is correct with getwd()
db <- dbConnect(SQLite(), './cuffData.db')
dbListTables(db)
expData <- (dbReadTable(db, 'geneExpDiffData'))

#This pulls data from the table that meets the following conditions.
m <- subset(expData, 
              expData$significant == 'yes' & 
              (expData$sample_1 == "OLNA" | 
              expData$sample_2 == "OLNA") &
              expData$log2_fold_change > 5 & 
              expData$log2_fold_change != "Inf")

#Display it in order of fold-change value			  
m[order(m$log2_fold_change, decreasing = TRUE), ]

plot(log(m$value_1), log(m$value_2))



geneExpDiff <- read.csv(file = './gene_exp.diff', header = TRUE, sep = "")
geneExpDiff[head(as.numeric(rownames(m))),]

lookUp <- as.numeric(rownames(m))
screamers <- geneExpDiff[lookUp,]

library(stringr)
s1 <- str_split_fixed(m$transcript, ":", 2)
m1$transcript <- s1[,1]
m <- m1;

OLNAup <- findSimilar(cuffdiff, c(1,1,1,1,1,1000,1,1,1,1,1), n=20)
expressionPlot(OLNAup, logMode = TRUE, showErrorbars = TRUE)

testing = NULL
testing = data.frame(replicate(2, sample(0:5, 20, rep=TRUE)))
testing$row = seq.int(1,20)
testing$OLNA = dtf$OLNA
testing = testing[,3:4]


o = expData[is.element(expData$gene_id, lookUp),]


o = subset(o, (o$sample_1 == 'OLNA' | o$sample_2 == 'OLNA') & o$significant == 'yes')
lookUp = unique(o[,1])
refGenes = geneExpDiff[is.element(geneExpDiff$test_id, lookUp),]
splitLocus = str_split_fixed(refGenes$locus, ":", 2)
refGenes$locus = splitLocus[,1]
write.table(as.factor(refGenes$locus), file = './OLNAUpRegulatedTranscriptList.txt', row.names = FALSE)

#Find genes based on expression profile, sort and pull the most extreme
OLNAdown <- findSimilar(cuffdiff, c(1000, 1000, 1000, 1000, 1000, 1, 1000, 1000, 1000, 1000, 1000), n=20)
OLNAdown@fpkm[order(OLNAdown@fpkm$sample_name),]
i <- getGene(cuffdiff, "XLOC_024117")

