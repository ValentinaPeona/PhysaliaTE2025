## Script to get basic info about a repeat library in fasta format
# Usage: Rscript --vanilla librarySummary.R input.lib
# it writes the results in input.lib.info

# 15 May 2023 - Valentina Peona

args = commandArgs(trailingOnly=TRUE)

filename = args[1]

## General info: number of sequences, max, min and mean length, nucleotide composition
## Info per category: number of sequences, max, min and mean length, nucleotide composition

### FUNCTIONS
sequenceLength = function(lens){

	output = paste("Maximum length:", max(lens), "bp." , "\nMinimum length:", min(lens), "bp.", "\nMean length:", mean(lens), "bp\n")

	return(output)

}

nucleotideComposition = function(library){

	A = lengths(regmatches(library[sequences,1], gregexpr("A|a", library[sequences,1])))
	T = lengths(regmatches(library[sequences,1], gregexpr("T|t", library[sequences,1])))
	C = lengths(regmatches(library[sequences,1], gregexpr("C|c", library[sequences,1])))
	G = lengths(regmatches(library[sequences,1], gregexpr("G|g", library[sequences,1])))
	N = lengths(regmatches(library[sequences,1], gregexpr("N|n", library[sequences,1])))

	output = paste("Number of A:", sum(A), "bp\nNumber of C:", sum(C), "bp\nNumber of G:", sum(G), "bp\nNumber of T:", sum(T), "bp\nNumber of N:", sum(N), "bp\nPercentage of A:", (sum(A)/sum(lens))*100, "%\nPercentage of C:", (sum(C)/sum(lens))*100, "%\nPercentage of G:", (sum(G)/sum(lens))*100, "%\nPercentage of T:", (sum(T)/sum(lens))*100, "%\nPercentage of N:", (sum(N)/sum(lens))*100, "%\n")

	return(output)
}


###

# read input fasta
library = read.table(file = filename, stringsAsFactors = FALSE, header = FALSE, comment.char = "", sep = "\t")

# lines of headers
headers = seq(from = 1, to = nrow(library), by = 2)

# lines of sequences
sequences = seq(from = 2, to = nrow(library), by = 2)

# write to file
outfile = paste0(filename, ".info")
sink(outfile)

cat("Library:", filename, "\n")

lens = nchar(library[sequences,])

cat("General length of consensus sequences:\nThe library comprises a total of", sum(lens), "bases in", length(lens), "sequences\n")

cat(sequenceLength(lens))

cat("General nucleotide composition:\n")

cat(nucleotideComposition(library))

cat("-------------------\n")

categories = sub(pattern = " .*", x = sub(pattern = "^>.*#", x = sapply(strsplit(x = library[headers,1], split = "/"), "[", 1), replacement = ""), replacement = "")
uniquecat = unique(categories)

cat("Info by category of repeats\n")
cat("Repeat categories:\n")
cat(uniquecat, "\n")

for(i in 1:length(uniquecat)){

	cat("\n", uniquecat[i], "\n")
	cat("Number of", uniquecat[i], ":", sum(categories == uniquecat[i]), "\n")
	lens = nchar(library[sequences,][categories == uniquecat[i]])
	cat(sequenceLength(lens))

	subLib = vector(length = 2 * length(which(categories == uniquecat[i])))
	subLib[seq(from = 1, to = length(subLib), by = 2)] = library[headers,1][which(categories == uniquecat[i])]
	subLib[seq(from = 2, to = length(subLib), by = 2)] = library[sequences,1][which(categories == uniquecat[i])]
	subLib = as.data.frame(subLib)
	cat(nucleotideComposition(subLib))
	cat("-------------------\n")

}

sink()
