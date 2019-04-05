library(tidyverse)
library(cowplot)
library(ggthemes)
library(Biostrings)

assy_tbl_from_fasta <- function(fasta_filename, assembly) {
  bss <- readBStringSet(fasta_filename)
  tbl <- tibble(Length=lengths(bss), Name=word(names(bss), 1)) %>%
    arrange(desc(Length)) %>%
    mutate(Cumsum=cumsum(Length), Assembly=assembly, ContigNo=row_number())
  rm(bss)
  gc()
  tbl
}

tbl <- map2(list('scol_trim.contigs.fasta', 'scol_trim.unitigs.fasta'),
     list('contigs', 'unitigs'),
     assy_tbl_from_fasta) %>%
  bind_rows()

summaries <- tbl %>%
  group_by(Assembly) %>%
  summarize(Length=sum(Length),
            NumContigs=n(),
            Labels=paste(unique(Assembly),'\n',n(),"contigs\n",round(sum(Length)/1e6), "Mb"))

ggplot(tbl, aes(x=ContigNo,y=Cumsum,col=Assembly)) + geom_line() +
    geom_point(data=summaries,mapping=aes(NumContigs,Length),color='black') +
    geom_text(data=summaries,mapping=aes(NumContigs,Length,label=Labels),color='black',hjust = 0,vjust=1,nudge_x = 200, nudge_y = 2000) +
    scale_color_ptol('Assembly') + theme(legend.justification = c(1, 1), legend.position = c(1, 1))  +
    coord_cartesian(xlim=c(1,70000))

