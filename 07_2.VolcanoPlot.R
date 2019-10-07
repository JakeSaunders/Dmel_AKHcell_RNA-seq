# ---------------------------------------------
# set working dir, load packages, and load objects ----
setwd("C:/Users/saundecj/Google Drive/fly/AKHcells.rnaSeq")

## Load objects
load("outs/06.dds.rld.res.RData")
sig.gene <- read.csv(file = "outs/06.padj.less.0.01.csv")

## load packages 
#install.packages("BiocManager")
library(BiocManager)
#install(pkgs = "DESeq2")
library("DESeq2")
library(scales)


# ---------------------------------------------
# rename objects ----
diff.exp <- res
# ---------------------------------------------
# Save plot 

pdf("plots/07_2.volcano.plot.pdf", width = 10,height = 8.5)

# ---------------------------------------------
# make volcano plot ----
    with(
        diff.exp,
        plot(
            x = log2FoldChange,
            y = -log10(padj),
            pch = 20,
            main = "Volcano plot: Fed vs Starved AKH cells",
            col=alpha("black",0.25),
            xlim=c(-10,10),
            bty="n"
        )
    )

# ---------------------------------------------
# add color to sig points ----

    # light green for log2FoldChange > 2
    with(
        subset(diff.exp, abs(log2FoldChange)>2 ), 
        points(
            log2FoldChange, 
            -log10(padj), 
            pch=20, 
            col=alpha("darkolivegreen3",0.25)
        )
    )   
   
    # dark green for log2FoldChange > 5 about 50fold
    with(
        subset(diff.exp, abs(log2FoldChange)>5 ), 
        points(
            log2FoldChange, 
            -log10(padj), 
            pch=20, 
            col="darkgreen"
        )
    )

    # gold for padj < 0.1
    with(
        subset(diff.exp, padj< 0.1 ), 
        points(
            log2FoldChange, 
            -log10(padj), 
            pch=20, 
            col="gold"
        )
    )    
    
    # orange for padj < 0.05
    with(
        subset(diff.exp, padj< 0.05 ), 
        points(
            log2FoldChange, 
            -log10(padj), 
            pch=20, 
            col="darkorange"
        )
    )
    
    # red for padj < 0.001
    with(
        subset(diff.exp, padj< 0.001 ), 
        points(
            log2FoldChange, 
            -log10(padj), 
            pch=20, 
            col="red"
        )
    )
#---------------------------------------------------------
# label points ----
#    install.packages("calibrate")
    library(calibrate)

# Label high expressers points with the textxy function 
    with(
        sig.gene, 
        textxy(log2FoldChange, -log10(padj), labs=X, cex=.5)
    )

# ------------------------------------------------------------
# legend ----
legend(
    x = "topleft",
    legend = c(
        "2 fold change", 
        "50 fold change", 
        "p < 0.1",
        "p < 0.05", 
        "p < 0.001"
    ),
    pch = 20,
    col = c(
        "darkolivegreen3",
        "darkgreen",
        "gold",
        "darkorange",
        "red"
    ),bty = "n",cex=0.8)

# --------------------------------------------------
# Close pdf plot ----
dev.off()

# helpful: http://www.gettinggeneticsdone.com/2014/05/r-volcano-plots-to-visualize-rnaseq-microarray.html ----

    
    
    