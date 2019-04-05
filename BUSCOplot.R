# plot BUSCO results 

library(ggplot2)
library(data.table)
library(plyr)
library(forcats)
library(ggrepel)


BUSCOtable <- function(path, name){
  # load data
  data <- fread(cmd=paste('grep -v "#"', path))
  
  # add col for busco run name
  data$name <- rep(name, length(data$V2))
  
  # drop first col
  data <- data[,-1]
  
  # rename factors
  data$recat <- revalue(data$V3, replace=c("Complete BUSCOs (C)"="complete", 
                                           "Complete and single-copy BUSCOs (S)"="complete (single-copy)",
                                           "Complete and duplicated BUSCOs (D)"="complete (duplicated)",
                                           "Fragmented BUSCOs (F)"="fragmented",
                                           "Missing BUSCOs (M)"="missing",
                                           "Total BUSCO groups searched"="total"))
  
  # exclude total and comp
  data <- subset(data, !(recat %in% c('total', 'complete')))
  
  # as factor, explicit order
  data$recat <- factor(data$recat, levels = rev(c("complete (single-copy)",
                                                              "complete (duplicated)",
                                                              "fragmented",
                                                              "missing")))
  return(data)
}


canu <- BUSCOtable(path = '/home/lukas/cube_scratch/02_BUSCO/run_busco_canu/short_summary_busco_canu.txt',
                  name = 'canu')

falcon <- BUSCOtable(path = '/home/lukas/cube_scratch/02_BUSCO/run_busco_falcon/short_summary_busco_falcon.txt',
                   name = 'falcon')

miniasm <- BUSCOtable(path = '/home/lukas/cube_scratch/02_BUSCO/run_busco_miniasm2/short_summary_busco_miniasm2.txt',
                   name = 'miniasm racon')

canu_augustus <- BUSCOtable(path = '/home/lukas/cube_scratch/11_BUSCO_proteins/run_busco_sc1_canu_augustus/short_summary_busco_sc1_canu_augustus.txt',
                            name = 'canu_augustus_proteins')

# merge the BUSCOs
BUSCO <- rbind(canu, falcon, miniasm, canu_augustus)


# make some labels - calc y pos and center - painful
BUSCO <- plyr::arrange(BUSCO, name, forcats::fct_rev(recat))
BUSCO <- ddply(BUSCO, "name", transform, label_y=cumsum(V2) - 0.5*V2)

# stacked bar plot
p <- ggplot(data=BUSCO, aes(y=V2, x=name, fill=recat)) + 
  geom_bar(stat = 'identity') +
  scale_fill_brewer() +
  geom_text_repel(aes(y=label_y, label=paste(round(V2/978*100,2), '%')), size=2.5, direction='y', alpha=0.6) +
  coord_flip() +
  labs(fill='') +
  ggtitle('BUSCO scores for metazoan gene set (n = 978)') +
  ylab("genes") +
  xlab('') +
  theme_minimal()
p
#ggsave(p, file="~/cube_home/02_BUSCO/BUSCO_scores.pdf", width=7,height=3)




