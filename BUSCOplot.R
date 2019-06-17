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


#canu <- BUSCOtable(path = '/home/lukas/cube_scratch/02_BUSCO/run_busco_canu/short_summary_busco_canu.txt',
#                  name = 'canu')

#falcon <- BUSCOtable(path = '/home/lukas/cube_scratch/02_BUSCO/run_busco_falcon/short_summary_busco_falcon.txt',
#                   name = 'falcon')

#miniasm <- BUSCOtable(path = '/home/lukas/cube_scratch/02_BUSCO/run_busco_miniasm2/short_summary_busco_miniasm2.txt',
#                   name = 'miniasm racon')

#canu_augustus <- BUSCOtable(path = '/home/lukas/cube_scratch/11_BUSCO_proteins/run_busco_sc1_canu_augustus/short_summary_busco_sc1_canu_augustus.txt',
#                            name = 'canu_augustus_proteins')

#canu_r0p3_c0p1 <- BUSCOtable(path = '/home/lukas/cube_scratch/02_BUSCO/run_busco_canu_r0p3_c0p1/short_summary_busco_canu_r0p3_c0p1.txt',
#                            name = 'canu_r0p3_c0p1')

#2PB_canu_r0p3_c0p045 <- BUSCOtable(path = '/home/lukas/cube_scratch/processing/02_BUSCO/run_busco_canu_r0p3_c0p045/short_summary_busco_canu_r0p3_c0p045.txt',
#                            name = '2PB_canu_r0p3_c0p045')

canu <- BUSCOtable(path = '/scratch/weilguny/analysis/02_busco/run_busco_scol3PB_r0p3_c0p045/short_summary_busco_scol3PB_r0p3_c0p045.txt',
                            name = '3PB_canu')

#3PB_canu_r0p3_c0p06 <- BUSCOtable(path = '/home/lukas/cube_scratch/analysis/02_busco/run_busco_scol3PB_r0p3_c0p06/short_summary_busco_scol3PB_r0p3_c0p06.txt',
#                            name = '3PB_canu_r0p3_c0p06')

#3PB_falcon <- BUSCOtable(path = '/home/lukas/cube_scratch/analysis/02_busco/run_busco_scol3PB_falcon/short_summary_busco_scol3PB_falcon.txt',
#                            name = '3PB_falcon')

canu_redundans <- BUSCOtable(path = '/scratch/weilguny/analysis/02_busco/run_busco_redundans/short_summary_busco_redundans.txt',
					name = '3PB_canu_redundans') 



# merge the BUSCOs
#BUSCO <- rbind(canu, falcon, miniasm, canu_augustus, canu_r0p3_c0p045, canu_r0p3_c0p1)

#BUSCO <- rbind(2PB_canu_r0p3_c0p045, 3PB_canu_r0p3_c0p045, 3PB_canu_r0p3_c0p06, 3PB_falcon)

BUSCO <- rbind(canu, canu_redundans)


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
ggsave(p, file="./BUSCO_3PB_redundans.pdf", width=7,height=2)




