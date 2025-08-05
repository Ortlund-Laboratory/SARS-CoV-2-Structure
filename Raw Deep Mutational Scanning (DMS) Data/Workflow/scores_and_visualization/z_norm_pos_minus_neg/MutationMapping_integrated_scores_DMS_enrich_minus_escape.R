library(ggplot2)
#library(hrbrthemes)
library(RColorBrewer)
library(plyr)
library(dplyr)
library(ggplot2)
library(tidyr)
library(seqinr)
require(ggseqlogo)
library(tidyverse)
library(grid)
library(readr)

################################################################################
###Rename 
################################################################################
antibody <- "Integrated_Nucleocapsid_DMS_Enrich_minus_escape"

mutations <- read_csv("/home/adkeith@Eu.Emory.Edu/DMS_Workflow/nucleocapsid/R_N-05_Myc_Mar_2022_version_to_use_for_paper/z_norm_pos_minus_neg/MycPos_Minus_MycNeg_Scores.csv")

wildtype <- read.fasta("N_Wuhan.fasta")

#####################################################################
#Combine duplicated mutations (same mutation with multiple barcodes)
#####################################################################

################################################################################
###Convert mutation "names" into: wildtype, site, and mutations
################################################################################
####################################################################
##Renaming for escape
####################################################################


fun_first <- function(x) {
  substring(x, 1,1)
}

fun_last <- function(x){
  substr(x, nchar(x), nchar(x))
}

fun_site <- function(x){
  substr(x, 2, nchar(x)-1)
}

site <- data.frame(sapply(mutations[1], fun_site))
plotmutations <- data.frame(sapply(mutations[1], fun_last))
colnames(site) <- c("site")
colnames(plotmutations) <- c("plotmutations")
mutations$site <- as.numeric(site[[1]])
mutations$plotmutations <- plotmutations[[1]]
score <- data.frame(mutations[2])

###############################################################################
###Add missing residues
###############################################################################
Wuhan <- read.fasta("N_Wuhan.fasta")

###First reference data
### If there are residues without any data, these have to be added manually
#1. Define the range of amino acids present
seq_range <- min(mutations$site):max(mutations$site)
missing_aa <- seq_range[!seq_range %in% unique(mutations$site)]
#2. Add the missing residue(s)
missing_data <- data.frame(site = missing_aa, 
                           name = paste(toupper(Wuhan[["N_Wuhan"]][missing_aa]), missing_aa, toupper(Wuhan[["N_Wuhan"]][missing_aa])),
                           mutation = toupper(Wuhan[["N_Wuhan"]][missing_aa]), 
                           wildtype = toupper(Wuhan[["N_Wuhan"]][missing_aa]))
complete_reference <- rbind.fill(mutations, missing_data)

#Expand the dataset to include NA values for synonymous amino acids
all <- complete_reference %>% expand(site, mutation)
# join with all, n will be NA for obs. in all that are not present in v
reference_trimmed = complete_reference %>% group_by_at(vars(wildtype, site, mutation)) %>% 
  right_join(all)

reference_trimmed <- reference_trimmed[
  with(reference_trimmed, order(site, mutation)),
]

mutations <- reference_trimmed

###################################################################
### Calculate abundance: n/N
###################################################################


###############################################################################
### Transform data using arcsine of the squareroot
### This transforms the data into something closer to a normal distribution
###############################################################################


################################################################################
###Z Normalization
################################################################################


################################################################################
### Generate a matrix for plotting sequence logos
### This matrix is also used to calculate total escape scores
################################################################################
#remove the extra rwos that only contain the site number and nothing else 
mutations <- mutations[which(!is.na(mutations$name)),]

logo_matrix <- matrix(ncol = 418,nrow=20)
row.names(logo_matrix) <- reference_trimmed$name[1:20]
colnames(logo_matrix) <- seq(2,ncol(logo_matrix)+1, 1)


for(i in 0:ncol(logo_matrix)-1){
  logo_matrix[,i+1] <- reference_trimmed$scaledscore[seq(from = 20*i+1, to=20*i+20, by = 1)]
}

################################################################################
##Plot heat maps
################################################################################


#Generate a dataframe containing the wild type amino acids 
#This will be used to mark wild type in the tiled heat map
tmp <- data.frame(sapply(mutations[1], fun_first), 
                  sapply(mutations[1], fun_site), 
                  sapply(mutations[1], fun_first))
frames <- na.omit(distinct(tmp))
colnames(frames) <- c("mutation", "site", "wildtype")

#Change data type to integer for "site"
#Required for proper heatmap plotting
frames$site <- as.integer(frames$site)
mutations$site <- as.integer(mutations$site)

#Set the order for amino acids in the heatmap (by aa property)
polar <- c("H", "C", "S", "T", "N", "Q")
nonpolar <- c("G", "A", "V", "L", "I", "M", "P")
aromatic <- c("F", "Y", "W")
positive <- c("K", "R")
negative <- c("D", "E")

aa_order <- c(negative,positive, polar, nonpolar, aromatic)

# Positive and negative value

#Residues 1-209
start = 0
end = 209

mut_range <- subset(mutations, site>start & site<end)
#mut_range <- c(586, 587, 588, 589, 590, 591, 593, 594, 600, 601, 602, 603, 604, 605, 616, 622, 623, 663, 664, 665, 685, 689, 690, 691, 694, 695, 698, 699, 706, 708, 709, 710)
mut_range <-as.data.frame(mut_range)
#Reverse the order so plotting will happen top to bottom:
mut_range$site <- as.factor(mut_range$site)
frames_range <- subset(frames, site>start & site<end)
frames_range$site <- as.factor(frames_range$site)


ggsave(filename = paste(antibody,"_All_Mutations_Integrated_part1.png", sep=""), 
       ggplot(mut_range, aes(site, plotmutations, size = 1)) + 
         geom_tile(color = "white",
                   fill = '#FCF0F0', #azure1 for blueish grey background
                   lwd = 0.1,
                   linetype = 1,
                   alpha = 1) +
         #scale_y_discrete(limits=rev) +
         ylim(rev(aa_order)) +
         geom_tile(data=mut_range, 
                   size=0,
                   height = 1,
                   fill='white', 
                   colour="white") +
         geom_point(aes(colour = scaledscore), alpha=1) +
         scale_colour_distiller(palette = "Spectral", direction = +1,
                                na.value = '#FCF0F0',
                                limits=c(-4,4)) +
         scale_size(range = c(0, 2)) +
         #The next section adds lightgrey tiles and black points for 
         #the wild type residues
         geom_point(inherit.aes = FALSE, 
                    data = frames_range,
                    aes(site, mutation),
                    shape=16,
                    size = 1.25,
                    colour="black") +
         xlab("Site") + ylab("Mutation") +
         scale_x_discrete(breaks = levels(mut_range$site)[c(rep(F,3),T,F)]) +
         theme(# Hide panel borders and remove grid lines
           panel.border = element_blank(),
           panel.grid.major = element_blank(),
           panel.grid.minor = element_blank(),
           # Change axis line
           #axis.line = element_line(colour = "black"),
           legend.key.size = unit(0.15, 'inch')
         ),
       width = 16, height = 3, dpi = 300, units = "in", device='png')

#Residues 210-end
start = 210
end = 420

mut_range <- subset(mutations, site>start & site<end)
#mut_range <- c(586, 587, 588, 589, 590, 591, 593, 594, 600, 601, 602, 603, 604, 605, 616, 622, 623, 663, 664, 665, 685, 689, 690, 691, 694, 695, 698, 699, 706, 708, 709, 710)
mut_range <-as.data.frame(mut_range)
#Reverse the order so plotting will happen top to bottom:
mut_range$site <- as.factor(mut_range$site)
frames_range <- subset(frames, site>start & site<end)
frames_range$site <- as.factor(frames_range$site)

site_vals <- as.numeric(levels(mut_range$site))
desired_vals <- site_vals[site_vals %% 5 == 0 & site_vals >= 215]

ggsave(filename = paste(antibody,"_All_Mutations_Integrated_part2.png", sep=""), 
       ggplot(mut_range, aes(site, plotmutations, size = 1)) + 
         geom_tile(color = "white",
                   fill = '#FCF0F0', #azure1 for blueish grey background
                   lwd = 0.1,
                   linetype = 1,
                   alpha = 1) +
         #scale_y_discrete(limits=rev) +
         ylim(rev(aa_order)) +
         geom_tile(data=mut_range, 
                   size=0,
                   height = 1,
                   fill='white', 
                   colour="white") +
         geom_point(aes(colour = scaledscore), alpha=1) +
         scale_colour_distiller(palette = "Spectral", direction = +1,
                                na.value = '#FCF0F0',
                                limits=c(-4,4)) +
         scale_size(range = c(0, 2)) +
         #The next section adds lightgrey tiles and black points for 
         #the wild type residues
         geom_point(inherit.aes = FALSE, 
                    data = frames_range,
                    aes(site, mutation),
                    shape=16,
                    size = 1.25,
                    colour="black") +
         xlab("Site") + ylab("Mutation") +
#         scale_x_discrete(breaks = levels(mut_range$site)[c(rep(F,3),T,F)]) +
         scale_x_discrete(breaks = as.character(desired_vals)) +
         theme(# Hide panel borders and remove grid lines
           panel.border = element_blank(),
           panel.grid.major = element_blank(),
           panel.grid.minor = element_blank(),
           # Change axis line
           #axis.line = element_line(colour = "black"),
           legend.key.size = unit(0.15, 'inch')
         ),
       width = 16, height = 3, dpi = 300, units = "in", device='png')


################################################################################
###Plot a mini-heatmap 
################################################################################
#Residues 0-end
start = 0
end = 420

mut_range <- subset(mutations, site>start & site<end)
#mut_range <- c(586, 587, 588, 589, 590, 591, 593, 594, 600, 601, 602, 603, 604, 605, 616, 622, 623, 663, 664, 665, 685, 689, 690, 691, 694, 695, 698, 699, 706, 708, 709, 710)
mut_range <-as.data.frame(mut_range)
#Reverse the order so plotting will happen top to bottom:
mut_range$site <- as.factor(mut_range$site)
frames_range <- subset(frames, site>start & site<end)
frames_range$site <- as.factor(frames_range$site)



ggsave(filename = paste(antibody,"_heatmap_mini.png", sep=""), 
       ggplot(mut_range, aes(site, plotmutations, size = 1)) + 
         geom_tile(color = "white",
                   fill = '#f7fbff', #azure1 for blueish grey background
                   lwd = 0.1,
                   linetype = 1,
                   alpha = 0.5) +
         geom_tile(data=mut_range, 
                   size=0,
                   height = 1,
                   fill='white', 
                   colour="black") +
         #scale_y_discrete(limits=rev) +
         ylim(rev(aa_order)) +
         geom_point(aes(colour = scaledscore), alpha=1) +
         scale_colour_distiller(palette = "Spectral", direction = +1,
                                na.value = '#f7fbff',
                                limits=c(-4,4)) +
         scale_size(range = c(0, 1.5)) +
         #The next section adds lightgrey tiles and black points for 
         #the wild type residues
         geom_point(inherit.aes = FALSE, 
                    data = frames,
                    aes(site, mutation),
                    shape=16,
                    size = 0.75,
                    colour="black") +
         scale_x_discrete(breaks = levels(mut_range$site)[c(rep(F,8),T,F)]) +
         theme_void() +
         theme(legend.position = "none") +
         scale_fill_distiller(palette = "Spectral", direction = +1, 
                              na.value = '#f7fbff',
                              limits=c(-4,4)),
       width = 24, height = 1.6, dpi = 300, units = "in", device='png')


################################################################################
### Plot average Escape histogram
################################################################################
average_escape <- as.data.frame(colMeans(logo_matrix, na.rm=T))
average_escape$site <- seq(2,419,1)
colnames(average_escape) <- c("escape", "site")

ggsave(filename = paste(antibody,"_average_escape_znorm.png", sep=""), 
       ggplot(average_escape, aes(site, escape)) + 
         geom_bar(stat = "identity") + 
         xlim(2,419) +
         ylim(0, max(average_escape$escape+5)) +
         theme_minimal() + 
         theme_void(),
       width = 12, height = .4, dpi = 300, units = "in", device='png')


#Individual escapes:
escape_fractions_out <- na.omit(mutations[,c(1,2)])
#escape_fractions_out$name <- paste(escape_fractions_out$wildtype, escape_fractions_out$site, escape_fractions_out$mutation, sep = "")
write.csv(escape_fractions_out, file = paste(antibody,"_escape_fractions_znorm.csv", sep = ""), row.names = FALSE)

#average escapes:
average_escape_out <- as.data.frame(average_escape$site)
average_escape_out$'average Escape Score' <- average_escape$escape
colnames(average_escape_out) <- c("Site", "average Escape Score")
write.csv(average_escape_out, file = paste(antibody,"_average_escape_znorm.csv", sep=""), row.names = FALSE)

################################################################################
#Pymol file for mapping average scores onto the structure
################################################################################
#Escapees
m_ <- mean(average_escape$escape, na.rm=T)
sd_ <- sd(average_escape$escape, na.rm=T)
cut_1 <- (m_ + 1*sd_)
cut_2 <- (m_ + 2*sd_)
cut_3 <- (m_ - 1*sd_)
cut_4 <- (m_ - 2*sd_)
escapees_C01 <- subset(average_escape,escape > cut_1 & escape < cut_2)
escapees_C02 <- subset(average_escape,escape > cut_2)
escapees_C03 <- subset(average_escape,escape < cut_3 & escape > cut_4)
escapees_C04 <- subset(average_escape,escape < cut_4)
top_stab <- na.omit(escapees_C02[,c(1,2)])
top_destab <- na.omit(escapees_C04[,c(1,2)])
write.csv(top_stab, file = paste(antibody,"_top_stabs.csv", sep=""), row.names = FALSE)
write.csv(top_destab, file = paste(antibody,"_top_destabs.csv", sep=""), row.names = FALSE)

#average_escape$escape <- ifelse(average_escape$escape < 0, 0, average_escape$escape)
average_escape <- average_escape %>%
  mutate(level = case_when(
    escape < cut_1 & escape > cut_3 ~ "grey",
    escape < cut_2 & escape > cut_1 ~ "red_1",
    escape > cut_2 ~ "red_2",
    escape < cut_3 & escape > cut_4 ~ "red_3",
    TRUE ~ "red_4"
  ))
average_escape$level <- as.factor(average_escape$level)

#Residues 0-209
start = 0
end = 209
average_escape1 <- subset(average_escape, site>start & site<end)
ggsave(filename = paste(antibody, "_color_average_escape_znorm_part1.png", sep=""), 
       ggplot(average_escape, aes(site, escape, fill = level)) + 
         geom_bar(stat="identity") + 
         geom_segment(aes(x=1,xend=209,y=cut_1, yend = cut_1), linetype = 3, size = 0.15) +
         geom_segment(aes(x=1,xend=209,y=cut_2, yend = cut_2), linetype = 3, size = 0.15) +
         geom_segment(aes(x=1,xend=209,y=cut_3, yend = cut_3), linetype = 3, size = 0.15) +
         geom_segment(aes(x=1,xend=209,y=cut_4, yend = cut_4), linetype = 3, size = 0.15) +
         #geom_segment(aes(x=215,xend=448,y=m_, yend = m_), linetype = 1, size = 0.03, alpha = 0.5) +
         scale_fill_manual(values = c("grey" = "#d9d9d9",
                                      "red_1" ="#90EE90",
                                      "red_2"="#228B22",
                                      "red_3"="#FFA500",
                                      "red_4"="#FF0000")) +
         xlim(1,209) +
         coord_cartesian(ylim = c(-2.5, 2.5)) +
         theme_void() +
         theme(legend.position = "none"),
       width = 12, height = 0.4, dpi = 300, units = "in", device='png')

#Residues 210-end
start = 210
end = 420
average_escape1 <- subset(average_escape, site>start & site<end)
ggsave(filename = paste(antibody, "_color_average_escape_znorm_part2.png", sep=""), 
       ggplot(average_escape, aes(site, escape, fill = level)) + 
         geom_bar(stat="identity") + 
         geom_segment(aes(x=210,xend=420,y=cut_1, yend = cut_1), linetype = 3, size = 0.15) +
         geom_segment(aes(x=210,xend=420,y=cut_2, yend = cut_2), linetype = 3, size = 0.15) +
         geom_segment(aes(x=210,xend=420,y=cut_3, yend = cut_3), linetype = 3, size = 0.15) +
         geom_segment(aes(x=210,xend=420,y=cut_4, yend = cut_4), linetype = 3, size = 0.15) +
         #geom_segment(aes(x=215,xend=448,y=m_, yend = m_), linetype = 1, size = 0.03, alpha = 0.5) +
         scale_fill_manual(values = c("grey" = "#d9d9d9",
                                      "red_1" ="#90EE90",
                                      "red_2"="#228B22",
                                      "red_3"="#FFA500",
                                      "red_4"="#FF0000")) +
         xlim(210,420) +
         coord_cartesian(ylim = c(-2.5, 2.5)) +
         theme_void() +
         theme(legend.position = "none"),
       width = 12, height = 0.4, dpi = 300, units = "in", device='png')

################################################################################
### Write a script for coloring residues in pymol
################################################################################

escapees_C01$site <- escapees_C01$site+1
escapees_C01$resi <- paste("or resi", escapees_C01$site, sep = " ")
reds01 <- paste(escapees_C01$resi)
escapees_C02$site <- escapees_C02$site+1
escapees_C02$resi <- paste("or resi", escapees_C02$site, sep = " ")
reds02 <- paste(escapees_C02$resi)
escapees_C03$site <- escapees_C03$site+1
escapees_C03$resi <- paste("or resi", escapees_C03$site, sep = " ")
reds03 <- paste(escapees_C03$resi)
escapees_C04$site <- escapees_C04$site+1
escapees_C04$resi <- paste("or resi", escapees_C04$site, sep = " ")
reds04 <- paste(escapees_C04$resi)


