# plot BUSCO results for three assemblies


library(ggplot2)
library(data.table)
library(plyr)
library(forcats)
library(ggrepel)

setwd("/home/lukas/cube_home/02_BUSCO/")

# read in the Busco output directly
canu <- fread(cmd='grep -v "#" run_busco_canu/short_summary_busco_canu.txt')
falcon <- fread(cmd='grep -v "#" run_busco_falcon/short_summary_busco_falcon.txt')
miniasm <- fread(cmd='grep -v "#" run_busco_miniasm2/short_summary_busco_miniasm2.txt')

# add a col for the assembler
canu$assembler <- rep('canu', length(canu$V2))
falcon$assembler <- rep('falcon', length(falcon$V2))
miniasm$assembler <- rep('miniasm', length(miniasm$V2))

# merge the asemblies
assemblies <- rbind(canu, falcon, miniasm)

# rename the factors
assemblies$recat <- revalue(assemblies$V3, replace=c("Complete BUSCOs (C)"="complete", 
                                          "Complete and single-copy BUSCOs (S)"="complete (single-copy)",
                                          "Complete and duplicated BUSCOs (D)"="complete (duplicated)",
                                          "Fragmented BUSCOs (F)"="fragmented",
                                          "Missing BUSCOs (M)"="missing",
                                          "Total BUSCO groups searched"="total"))

# drop the first empty column
assemblies <- assemblies[,-1]

# exclude the total number and complete
assemblies <- subset(assemblies, !(recat %in% c('total', 'complete')))

assemblies$recat <- factor(assemblies$recat, levels = rev(c("complete (single-copy)",
                                                        "complete (duplicated)",
                                                        "fragmented",
                                                        "missing")))

# make some labels - calc y pos and center - painful
assemblies <- plyr::arrange(assemblies, assembler, forcats::fct_rev(recat))
assemblies <- ddply(assemblies, "assembler", transform, label_y=cumsum(V2) - 0.5*V2)

# stacked bar plot
p <- ggplot(data=assemblies, aes(y=V2, x=assembler, fill=recat)) + 
  geom_bar(stat = 'identity') +
  scale_fill_brewer() +
  geom_text_repel(aes(y=label_y, label=paste(round(V2/978*100,2), '%')), size=2.5, direction='y', alpha=0.6) +
  coord_flip() +
  labs(fill='') +
  ggtitle('BUSCO scores for metazoan gene set (n = 978)') +
  ylab("genes") +
  xlab('') +
  theme_minimal()
ggsave(p, file="~/cube_home/02_BUSCO/BUSCO_scores.pdf", width=7,height=3)



  
