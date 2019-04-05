# plot REAPR results

library(ggplot2)
library(reshape)
library(data.table)

# grab the info
grabN <- function(path, N){
  Ns <- fread(cmd=paste("grep 'N", N, " =' ", path, sep=''))
  N_orig <- as.numeric(strsplit(names(Ns)[1], ' = ')[[1]][2])
  N_broken <- as.numeric(strsplit(as.character(Ns[1,1]), split = ' = ')[[1]][2])
  return(c(N_orig, N_broken))
}

grabNseq <- function(path){
  Nseqs <- fread(cmd=paste("grep 'Number of sequences'", path, sep=' '))
  Nseqs_orig <- as.numeric(Nseqs$V4)[1]
  Nseqs_broken <- as.numeric(Nseqs$V4)[2]
  frame <- data.frame('values'=c(Nseqs_orig, Nseqs_broken), 'fac'=c('01_original', '02_broken'))
  return(frame)
}

grabERR <- function(path){
  err <- fread(cmd=paste("grep 'Error free bases:'", path, sep=' '))
  errval <- as.numeric(gsub(pattern='%', replacement='', x=err$V4))
  errval_frame <- data.frame('error_free'=errval)
  return(errval_frame)
}

# get all Ns
grabAllN <- function(path){
  Nvals <- seq(50, 100, by=10)
  mat <- t(sapply(Nvals, grabN, path=path))
  fac <- c(rep('01_original', nrow(mat)), rep('02_broken', nrow(mat)))
  frame <- data.frame('N_val' = c(mat[,1], mat[,2]), 'fac' = fac, 'Nx'=Nvals)
  return(frame)
}

# get info from one assembly
reapr_summary <- function(path, name){
  # Ns info
  frame_N <- grabAllN(path)
  frame_N$assembler <- rep(name, nrow(frame_N))
  # number of seqs
  frame_Nseq <- grabNseq(path)
  frame_Nseq$assembler <- name
  # error free bases
  frame_err <- grabERR(path)
  frame_err$assembler <- name
  
  return(list(frame_N, frame_Nseq, frame_err))
}


data <- c(
falcon = "~/cube_scratch/03_REAPR/reapr_falcon_out/05.summary.report.txt",
miniasm = "~/cube_scratch/03_REAPR/reapr_miniasm_out/05.summary.report.txt",
canu_original = "~/cube_scratch/03_REAPR/reapr_canu_out/05.summary.report.txt",
canu_pilon_PE = "~/cube_scratch/08_reapr_new_assemblies/canu_pilon_reapr/05.summary.report.txt",
canu_racon = "~/cube_scratch/08_reapr_new_assemblies/canu_racon_reapr/05.summary.report.txt",
canu_r0p3_c0p045 = "~/cube_scratch/08_reapr_new_assemblies/scol_r0p03_c0p045_reapr/05.summary.report.txt",
canu_both_spikeins = "~/cube_scratch/04_mergebam/reapr_canu_out_spikeins/05.summary.report.txt "
)

NGs <- data.frame()
Nseqs <- data.frame()
Err <- data.frame()

for (i in seq(1,length(data))){
  tmp <- reapr_summary(data[i], name=names(data[i]))
  NGs <- rbind(NGs, tmp[[1]])
  Nseqs <- rbind(Nseqs, tmp[[2]])
  Err <- rbind(Err, tmp[[3]])
}


# plot1: change of N values after breaking
p <- ggplot(data=NGs, aes(x=Nx, y=N_val, color=fac)) + 
  geom_line() +
  facet_wrap( ~ assembler) +
  scale_color_manual(values=c('skyblue4', 'skyblue1')) +
  labs(color='') +
  ylab("contig length") +
  theme_minimal() 
#p


# plot2: change in number of seqs
q <- ggplot(data=Nseqs, aes(x=assembler, y=values, fill=fac)) +
  geom_bar(stat='identity', position='dodge') +
  scale_fill_manual(values=c('skyblue4', 'skyblue1')) +
  labs(fill='') +
  ylab("# contigs") +
  xlab('') +
  coord_flip() +
  theme_minimal() +
  theme(legend.position = 'bottom')
#q


# plot3: rate of error free bases
r <- ggplot(data=Err, aes(x=assembler, y=error_free)) +
  geom_point() +
  geom_bar(stat='identity', fill='black', width=0.1) +
  coord_flip() +
  ylab('Error free bases (%)') +
  xlab('') +
  theme_minimal()
#r


combo <- plot_grid(p, plot_grid(q, r, nrow=1), nrow=2)
#combo

ggsave(combo, file="~/cube_home/results/reapr_combo.pdf", width=7,height=7)

