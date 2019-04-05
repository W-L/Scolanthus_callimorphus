library(ggplot2)
library(reshape)
library(dplyr)
library(data.table)
library(cowplot)
library(argparse)

parser <- ArgumentParser()
parser$add_argument('--assemblies', required=TRUE, help='')
parser$add_argument('--inputNG', required=TRUE, help='')
parser$add_argument('--inputCS', required=TRUE, help='')
parser$add_argument('--output', required=TRUE, help='')
args <- parser$parse_args()

read.tcsv = function(file, header=TRUE, sep=" ", ...) {
  
  n = max(count.fields(file, sep=sep), na.rm=TRUE)
  x = readLines(file)
  
  .splitvar = function(x, sep, n) {
    var = unlist(strsplit(x, split=sep))
    length(var) = n
    return(var)
  }
  
  x = do.call(cbind, lapply(x, .splitvar, sep=sep, n=n))
  x = apply(x, 1, paste, collapse=sep) 
  out = read.csv(text=x, sep=sep, header=header, ...)
  return(out)
  
}

assy_raw <- args$assemblies
inNG <- args$inputNG
inCS <- args$inputCS
output <- args$output

# assy_raw <- "canu_orig,miniasm,falcon,canu_pilon_PE,canu_racon,canu_r0p3_c0p045,masurca"
# inNG <- '/home/lukas/cube_scratch/09_assembly_plot/six_assemblies.NG'
# inCS <- '/home/lukas/cube_scratch/09_assembly_plot/six_assemblies.cumsum'
# output <- '/home/lukas/cube_scratch/09_assembly_plot/six_assemblies'


assy <- strsplit(assy_raw, split=',')[[1]]

# plot NGx graph
ng <- read.csv(inNG, sep="")
# reshape the datafram and rename the factors
ng_long <- melt(ng, id.vars='NGx', measure.vars = names(ng)[-1] )

assy_vec <- assy
names(assy_vec) <- names(ng)[-1]

ng_long$assembler <- recode(ng_long$variable, !!!assy_vec)

g <- ggplot(data=ng_long, aes(x=NGx, y=value, color=assembler)) +
  geom_line() +
  #geom_point(size=0.5) +
  scale_color_brewer(palette='Dark2') +
  #scale_y_continuous(breaks=c(0,250000,500000,750000,1000000), labels=c('0 Mb', '0.25 Mb', '0.5 Mb', '0.75 Mb', '1 Mb')) +
  ylab('Contig length') +
  theme_minimal() +
  theme(legend.position = 'bottom')
#g


# plot cumsums
cumsums <- read.tcsv(inCS)
assys <- names(cumsums)
cumsums$n_contig <- 1:nrow(cumsums)

cumsums_long <- melt(cumsums, id.vars='n_contig', measure.vars = assys)

assy_vec <- assy
names(assy_vec) <- assys

cumsums_long$assembler <- recode(cumsums_long$variable, !!!assy_vec)

h <- ggplot(data=cumsums_long, aes(x=n_contig, y=value, color=assembler)) +
  geom_line() +
  scale_color_brewer(palette="Dark2") +
  #scale_y_continuous(breaks=c(0,2e+8,4e+8,6e+8,8e+8), labels=c('0 Mb', '200 Mb', '400 Mb', '600 Mb', '800 Mb')) +
  ylab('Cumulative sequence length (bp)') +
  xlab('#Contig') +
  theme_minimal() +
  theme(legend.position = 'bottom')
#h

# common legend
leg <- get_legend(g)
g <- g + guides(color=FALSE)
h <- h + guides(color=FALSE)

combo <- plot_grid(plot_grid(g, h, nrow=1, ncol=2), leg, nrow=2, ncol=1, rel_heights = c(0.8,0.2))
ggsave(combo, filename = paste(output, '.pdf', sep=''), width=9, height=6)

