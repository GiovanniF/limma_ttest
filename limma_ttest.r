
###STARTTEXTSCRIPT

library(RCurl)
library(genefilter)
library(limma)

setwd(tempdir())


##https://drive.google.com/file/d/0B__nP63GoFhMd042V09GcHFvVFE/edit?usp=sharing
x = getBinaryURL("https://docs.google.com/uc?export=download&id=0B__nP63GoFhMd042V09GcHFvVFE", followlocation = TRUE, ssl.verifypeer = FALSE)
writeBin(x, "std_var.txt", useBytes = TRUE)
std_var =  as.matrix(read.table("std_var.txt",
header = FALSE, sep = "\t", quote="\"", comment.char = ""
))

##https://drive.google.com/file/d/0B__nP63GoFhMc2tFT3hXVW52Q0E/edit?usp=sharing
x = getBinaryURL("https://docs.google.com/uc?export=download&id=0B__nP63GoFhMc2tFT3hXVW52Q0E", followlocation = TRUE, ssl.verifypeer = FALSE)
writeBin(x, "mean_var.txt", useBytes = TRUE)
mean_var = as.matrix(read.table("mean_var.txt",
header = FALSE, sep = "\t", quote="\"", comment.char = ""
))


NO_OF_GROUPS = 2

NO_OF_SAMP_PER_GROUP_GROUND_TRUTH = 100

NO_OF_GENES = NROW(mean_var)

gxexprs_GroundTruth = matrix(0, NO_OF_GENES, NO_OF_GROUPS*NO_OF_SAMP_PER_GROUP_GROUND_TRUTH)


std_var_list = list(std_var, std_var,seq(.01,.2, length.out= NO_OF_GENES))

simulated_mean = matrix(c(8,9),NO_OF_GENES,NO_OF_GROUPS, byrow=TRUE)

mean_var_list = list(mean_var, simulated_mean, simulated_mean)

main_list = list("micorarray mean and std", "means 8 and 9 and micorarray std",  "means 8 and 9 and uniform std from 0.01 to 0.2")

for(count_i in 1:3)
{

mean_var = mean_var_list[[count_i]]
std_var = std_var_list[[count_i]]


for(i in 1:NO_OF_GENES)
{
	if ((i%%1000)==1)
	{
		cat(i, "\n")
	}
	for(j in 1:NO_OF_GROUPS)
	{
		gxexprs_GroundTruth[i, NO_OF_SAMP_PER_GROUP_GROUND_TRUTH*(j-1) + 1:NO_OF_SAMP_PER_GROUP_GROUND_TRUTH]=
		rnorm(NO_OF_SAMP_PER_GROUP_GROUND_TRUTH, mean_var[i,j], std_var[i])
	}
}

final_choice = c(1,2)

Group = factor(sprintf("S%02d",
rep(final_choice,each=NO_OF_SAMP_PER_GROUP_GROUND_TRUTH)
))

t_test_pvals = rowttests(gxexprs_GroundTruth, fac=Group)$p.value



design <- model.matrix(~0+Group)
colnames(design) <- gsub("Group","",colnames(design))
contrastsX = sprintf("S01 - S%02d",2)

#fit <- lmFit(gxexprs_GroundTruth[, final_col_dest], design)
fit <- lmFit(gxexprs_GroundTruth, design)

contrast.matrix <- makeContrasts(contrasts=contrastsX,levels=design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)

TTable = topTable(fit2, coef=1, adjust="fdr", number = nrow(gxexprs_GroundTruth))
limma_pvals = TTable$P.Value[as.numeric(rownames(TTable))]


dev.new()
plot(-log10(limma_pvals), -log10(t_test_pvals), main=main_list[[count_i]])

}

###ENDTEXTSCRIPT



