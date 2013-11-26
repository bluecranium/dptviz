# ###################################################################
# Author: Jeffrey Bhasin <jeffb@case.edu>
# ###################################################################

# ###################################################################
# Main Report Functions

# -------------------------------------------------------------------
# ' Summary plots of data preprocessed in bins by iDPT
# '
# ' iDPT performs a preprocessing stage to convert a set of alignments in a BAM file to average coverage values in windows of a user-specified size genome-wide. This function saves a PDF file containing plots that summarize the distribution of these average coverage values and their relationship to CpG density on a specified reference genome.
# ' @param dptpath Path to a folder containing an iDPT run (the folder than contains the "DPT_ws" folder)
# ' @param chrs A vector of chromosome names to use (by default, the autosomes and sex chromosomes from hg19)
# ' @param bsgenome A BSgenome object for the reference genome used in the alignments contained in the input BAM files given to iDPT (i.e. for UCSC's hg19, you would load "library(BSgenome.Hsapiens.UCSC.hg19)" and give BSgenome.Hsapiens.UCSC.hg19 to this argument). This is needed to compute CpG density in the windows used by iDPT.
# ' @param file Filename to save the PDF report to (if NULL, will use the name of the folder as a prefix and save "prefix.reportPrepro.pdf" in the current working directory")
# ' @param ncore Number of concurrent threads to use
# ' @export
reportPrepro <- function(dptpath, chrs=c(sapply(seq(1,22),function(x) paste("chr",x,sep="")),"chrX","chrY"), bsgenome, file=NULL, ncore=1)
{
	# Set number of cores for %dopar%
	registerDoMC(ncore)

	# Set filename using dptpath if none given as argument
	if(is.null(file))
	{
		prefix <- str_split(dptpath, "/")[[1]]
		prefix <- prefix[length(prefix)]
		file <- paste(prefix, ".reportPrepro.pdf", sep="")
	}

	# Load data as matrix
	message("Loading data matrix")
	mat <- dptviz:::loadPrepro(dptpath, chrs)

	# Load bin locations as GRanges
	message("Extracting bin coordinates")
	bins.gr <- dptviz:::getBins(dptpath, chrs)

	# Calculate CpG density in reference genome for each bin
	message("Calculating CpG density for each bin")
	cpgctvec <- dptviz:::getCpG(bins.gr, bsgenome)

	# Make the average coverage into integer counts
	mat.r <- round(mat,0)

	# Start Plotting

	# Open PDF graphics device
	pdf(file=file, width=10.5, height=8)
	message(paste("Plotting as PDF: ",file, sep=""))

	# PLOT: Proportional Bar Graph of Zero to Nonzero
	df <- data.frame(sample=rownames(mat), zero=rowSums(mat==0), "less.than.10"=rowSums((mat<=10)&(mat!=0)), "from.11.to.100"=rowSums((mat>10)&(mat<=100)), "greater.than.100"=rowSums(mat>100))
	df <- melt(df)
	print(ggplot(df, aes(x=sample, y=value, fill=variable)) + geom_bar(stat="identity") + dptviz:::ggnice() + coord_flip() + labs(title="Proportion of Windows with Zero Reads", x="Sample", y="Number of Windows"))

	# PLOT: Kernel Density Plot of Read Counts by Windows
	dat <- melt(mat[,colSums(mat)>0])
	dat <- dat[dat$value>0,]
	print(ggplot(dat, aes(x=value, y=..count..)) + geom_density(aes(group=X1, colour=X1, fill=X1), alpha=0.3) + dptviz:::ggnice() + labs(title="Unrounded Read Distribution, Windows with > 0 Reads", x="log10(Read Count)", y="Number of Windows with Count"))

	# PLOT: Zoomed Kernel Density
	datsub <- dat[dat$value<25,]
	print(ggplot(datsub, aes(x=value, y=..count..)) + geom_density(aes(group=X1, colour=X1, fill=X1), alpha=0.3) + dptviz:::ggnice() + labs(title="Unrounded Read Distribution, Windows with > 0 Reads", x="Read Count", y="Number of Windows with Count"))

	# PLOT: Histogram of exact counts at read levels near the initfilter cutoff zone
	df <- foreach(i=1:nrow(mat), .combine="rbind") %dopar%
	{
		matsub <- mat.r[i,(mat.r[i,]>0)&(mat.r[i,]<=20)]
		h <- hist(matsub, breaks=c(1:21), include.lowest=TRUE, right=FALSE, plot=FALSE)
		bins <- h$breaks[1:(length(h$breaks)-1)]
		bins <- factor(bins, levels=bins)
		ret <- data.frame(sample=rownames(mat)[i], reads=bins, windows=h$counts)
		ret
	}
	print(ggplot(df, aes(x=reads, y=windows, fill=sample)) + geom_bar(stat="identity") + dptviz:::ggnice() + labs(title="Histogram of Rounded Read Counts", x="Read Count", y="Number of Windows at Count") + facet_grid(sample ~ .))

	# PLOT: Histogram of CpG count in windows with and without signal
	df1 <- foreach(i=1:nrow(mat), .combine="rbind") %dopar%
	{
		wins <- mat.r[i,]>0
		cg <- cpgctvec[wins]
		h <- hist(cg, breaks=0:(max(cpgctvec)+1), include.lowest=TRUE, right=FALSE, plot=FALSE)
		bins <- h$breaks[1:(length(h$breaks)-1)]
		bins <- factor(bins, levels=bins)
		ret <- data.frame(sample=rownames(mat)[i], type="greater_than_zero", ncpg=bins, nwindows=h$counts, per=(h$counts/sum(wins))*100)
		ret
	}
	df2 <- foreach(i=1:nrow(mat), .combine="rbind") %dopar%
	{
		wins <- mat.r[i,]==0
		cg <- cpgctvec[wins]
		h <- hist(cg, breaks=0:(max(cpgctvec)+1), include.lowest=TRUE, right=FALSE, plot=FALSE)
		bins <- h$breaks[1:(length(h$breaks)-1)]
		bins <- factor(bins, levels=bins)
		ret <- data.frame(sample=rownames(mat)[i], type="zero", ncpg=bins, nwindows=h$counts, per=(h$counts/sum(wins))*100)
		ret
	}
	df3 <- foreach(i=1:nrow(mat), .combine="rbind") %dopar%
	{
		wins <- mat.r[i,]>3
		cg <- cpgctvec[wins]
		h <- hist(cg, breaks=0:(max(cpgctvec)+1), include.lowest=TRUE, right=FALSE, plot=FALSE)
		bins <- h$breaks[1:(length(h$breaks)-1)]
		bins <- factor(bins, levels=bins)
		ret <- data.frame(sample=rownames(mat)[i], type="greater_than_3", ncpg=bins, nwindows=h$counts, per=(h$counts/sum(wins))*100)
		ret
	}
	df <- rbind(df1,df3,df2)
	print(ggplot(df, aes(x=ncpg, y=per, fill=sample)) + geom_bar(stat="identity") + dptviz:::ggnice() + labs(title="Percent of Zero and Nonzero Signal Windows at CpG Counts", x="Reference Genome CpG Count", y="Percent of Windows at CpG Count") + facet_grid(sample ~ type))

	# PLOT: CpG to nReads Relationship by Sample and Pooled Overall
	df <- foreach(i=1:nrow(mat), .combine="rbind") %dopar%
	{
		df <- data.frame(sample=rownames(mat)[i], ncpg=cpgctvec[cpgctvec>0], reads=mat.r[i,cpgctvec>0])
		df$logreads <- log10(df$reads)
		df <- df[df$reads>3,]
		df <- df[df$reads<100,]
		df$ncpg.fact <- factor(df$ncpg)
		df
	}
	print(ggplot(df, aes(x=ncpg.fact, y=reads, fill=sample)) + geom_boxplot() + ggnice() + facet_grid(sample ~ .) + labs(title="Read Count Distributions at Increasing CpG Counts", x="CpG Count", y="Read Count"))
	cor <- cor(df$ncpg,df$reads,use="pairwise.complete.obs")
	reg <- lm(df$reads ~ df$ncpg)
	print(smoothScatter(df$ncpg, df$reads,nrpoints=1000,main=paste("nCpG vs Reads, Pooled for All Samples", " (cor=",round(cor,digits=4),")",sep=""),xlab="CpG", ylab="Reads"))
	print(abline(reg))

	# Close Graphics Device
	dev.off()

	NULL
}
# -------------------------------------------------------------------

# -------------------------------------------------------------------
# ' Summary plots of sites detected by iDPT

# -------------------------------------------------------------------

# -------------------------------------------------------------------
# ' Annotation (via the goldmine package) of the genomic contexts of sites detected by iDPT

# -------------------------------------------------------------------

# ###################################################################

# ###################################################################
# Utility Functions

# Nice clean ggplot2 theme
ggnice <- function()
{
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.text = element_text(color="black"))
}

# Loads bcrvg.RData and converts object to a matrix (rows=samples, cols=windows genome-wide in order of chrs)
loadPrepro <- function(dptpath, chrs)
{
	# Load iDPT workspace
	bcvrg.path <- paste(dptpath, "/DPT_ws/Preprocess/bcvrg.RData", sep="")
	load(bcvrg.path)

	# Build matrix of coverage values
	samples <- names(bcvrg$data)

	mat <- foreach(sample=samples, .combine=rbind, .verbose=FALSE) %dopar%
	{
		# get genome-wide vector
		foreach(chr=chrs, .combine=c) %do%	
		{
			print(paste("Doing",sample, "chr", chr, sep=" "))
			vec <- as.vector(bcvrg$data[[sample]][[chr]])
			vec
		}
	}
	rownames(mat) <- samples
	rownames(mat) <- str_replace(rownames(mat),".bam","")
	mat
}

# Returns GenomicRanges object of each bin iDPT used
getBins <- function(dptpath, chrs)
{
	bcvrg.path <- paste(dptpath, "/DPT_ws/Preprocess/bcvrg.RData", sep="")
	load(bcvrg.path)
	bcvrg.gr2 <- GRanges(seqnames=seqnames(bcvrg.rd), ranges=IRanges(start=start(bcvrg.rd), end=end(bcvrg.rd)))
	gbins <- bcvrg.gr2[seqnames(bcvrg.gr2) %in% chrs]
	seqlevels(gbins) <- chrs
	gbins <- gbins[order(seqnames(gbins))]
}

# Returns vector of CpG frequency in each region of a GenomicRanges object
getCpG <- function(gr, bsgenome)
{	
	# Parallel by seqnames (chrs)
	out <- foreach(chr=as.character(unique(seqnames(gr))), .combine="c") %dopar%
	{
		print(paste("Doing",chr,sep=" "))
		seq <- getSeq(bsgenome, names=gr[seqnames(gr)==chr])
		difreq <- dinucleotideFrequency(seq, as.prob=FALSE)
		cpgfreq <- difreq[,colnames(difreq)=="CG"]
		cpgfreq
	}
	out
}
# ###################################################################
