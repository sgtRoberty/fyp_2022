library(ggplot2)
library(dplyr)
Ntip(ultra_aa_mitogenome_Malaysia_tree)
Ntip(ultra_aa_mitogenome_Panama_tree)

malaysia_family_df <- data.frame(table(Malaysia_metadata$family))
colnames(malaysia_family_df) <- c("family_name", "malaysia_count")

panama_family_df <- data.frame(table(Panama_metadata$family))
colnames(panama_family_df) <- c("family_name", "panama_count")

library(tidyr)
mp_family_top12_table <- merge(top_n(malaysia_family_df, n=12, malaysia_count), 
                               top_n(panama_family_df, n=12, panama_count), 
                               all = T)
mp_family_top12_table <- gather(mp_family_top12_table, 
                                country, count,
                                malaysia_count:panama_count,
                                factor_key=TRUE)
mp_family_top12_table[is.na(mp_family_top12_table)] <- 0
mp_family_top12_table$country <- factor(mp_family_top12_table$country, 
                                        levels = c('panama_count',
                                                   'malaysia_count'))

species_compo_top12_plot <-
  ggplot(mp_family_top12_table, aes(x=reorder(family_name, count), 
                                    y=count, fill=country)) +
  geom_bar(stat='identity', position=position_dodge(), 
           width=0.8) +
  labs(x = "Family name", y = "Count") +
  geom_text(aes(label=count), hjust=-0.2,
            position = position_dodge(0.8), 
            color="black", size=3) +
  scale_fill_manual(values=c('orange','forestgreen'),
                    name = "Country", 
                    labels = c("Panama", "Malaysia")) + 
  theme_minimal()+ 
  theme(axis.text=element_text(size=10),
        axis.title=element_text(size=12, face="bold"),
        axis.text.x = element_text(vjust = 2.5),
        axis.title.x = element_text(vjust = -0.5),
        axis.title.y = element_text(vjust = 2.5)) +
  coord_flip()
ggsave("../results/FYP_figs/species_compo_top12_plot.pdf", 
       plot = species_compo_top12_plot, width = 8, height = 6)


top_n(malaysia_family_df, n=10, malaysia_count) %>%
  ggplot(., aes(x=reorder(family_name, malaysia_count), y=malaysia_count)) +
  geom_bar(stat='identity', fill = "forestgreen") +
  labs(x = "Family name", y = "Count") +
  geom_text(aes(label=malaysia_count), hjust=-0.3, 
            color="black", size=4) + 
  theme_minimal()+ 
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14, face="bold"),) +
  coord_flip()

top_n(panama_family_df, n=10, panama_count) %>%
  ggplot(., aes(x=reorder(family_name, panama_count), y=panama_count))+
  geom_bar(stat='identity', fill = "orange") +
  labs(x = "Family name", y = "Count") + 
  geom_text(aes(label=panama_count), hjust=-0.3, 
            color="black", size=4) + 
  theme_minimal()+ 
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14, face="bold"),) +
  coord_flip()

