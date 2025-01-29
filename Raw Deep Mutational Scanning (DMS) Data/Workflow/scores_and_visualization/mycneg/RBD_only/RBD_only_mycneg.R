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
antibody <- "MycNeg"
ab_name <- antibody

escape_counts <- read_tsv("/home/adkeith@Eu.Emory.Edu/DMS_Workflow/nucleocapsid/R_N-05_Myc_Mar_2022_version_to_use_for_paper/exp3_mycneg/RBD_only/sp_nwt_escape.txt")
input_counts <- read_tsv("/home/adkeith@Eu.Emory.Edu/DMS_Workflow/nucleocapsid/R_N-05_Myc_Mar_2022_version_to_use_for_paper/exp3_mycneg/RBD_only/sp_nwt_input.txt")

escape_counts <- escape_counts[,2:3]
input_counts <- input_counts[,2:3]

################################################################################
### Histogram for read counts
################################################################################
ggsave(filename = paste(antibody, "_Number_of_Reads_escaped.png", sep=""), 
       ggplot(escape_counts, aes(x=count), title = "Mutations reads")+
         geom_histogram(color = "black", fill = "grey", ) +
         xlim(0,quantile(escape_counts$count, probs = .95, na.rm = T)) +
         ylim(0,5000)+
         # geom_vline(aes(xintercept=mu_bc_depth),
         #            linetype="dashed",
         #            color = "red") + 
         #geom_text(aes(label = mu), x = 3, y = 3, vjust = "inward", hjust = "inward") +
         ggtitle("Escape Read Counts") +
         xlab("# Reads") + 
         ylab("Count") +
         theme_bw(base_size = 10),
       width = 3, height = 2, dpi = 300, units = "in", device='png')

ggsave(filename = paste(antibody, "_Number_of_Reads_Myc.png", sep=""), 
       ggplot(input_counts, aes(x=count), title = "Mutations reads")+
         geom_histogram(color = "black", fill = "grey") +
         xlim(0,quantile(input_counts$count, probs = .99, na.rm = T)) +
         ylim(0,30000)+
         # geom_vline(aes(xintercept=mu_bc_depth),
         #            linetype="dashed",
         #            color = "red") + 
         #geom_text(aes(label = mu), x = 3, y = 3, vjust = "inward", hjust = "inward") +
         ggtitle("Reference Read Counts") +
         xlab("# Reads") + 
         ylab("Count") +
         theme_bw(base_size = 10),
       width = 3, height = 2, dpi = 300, units = "in", device='png')

#####################################################################
#Combine duplicated mutations (same mutation with multiple barcodes)
#####################################################################
escape_combined <- escape_counts %>% 
  group_by(mutation) %>% 
  summarise_all(funs(sum))

reference_combined <- input_counts %>% 
  group_by(mutation) %>% 
  summarise_all(funs(sum))

colnames(reference_combined) <- c("name", "depth")
colnames(escape_combined) <- c("name", "depth")

##Adjust all mutations in the reference set with less than the 5th percentile of reads
#depth_cutoff  <- quantile(reference_combined$depth, probs = .01, na.rm = T)
#reference_combined$depth <- ifelse(reference_combined$depth < depth_cutoff, depth_cutoff, reference_combined$depth)

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

wildtype <- data.frame(sapply(escape_combined[1], fun_first))
#colnames(wildtype) <- c("wildtype")
escape_combined$wildtype <- wildtype[[1]]

mutation <- data.frame(sapply(escape_combined[1], fun_last))
colnames(mutation) <- c("mutation")
escape_combined$mutation <- as.factor(mutation[[1]])

site <- data.frame(sapply(escape_combined[1], fun_site))
colnames(site) <- c("site")
escape_combined$site <- as.numeric(site[[1]])
#Reorder the columns (name, wild type, site, mutation, depth)
escape_combined <- escape_combined[,c(1,3,5,4,2)]

###Order reference_combined by site and mutation
escape_combined <- escape_combined[
  with(escape_combined, order(site, mutation)),
]



####################################################################
##Renaming for reference
####################################################################
wildtype <- data.frame(sapply(reference_combined[1], fun_first))
#colnames(wildtype) <- c("wildtype")
reference_combined$wildtype <- wildtype[[1]]

mutation <- data.frame(sapply(reference_combined[1], fun_last))
colnames(mutation) <- c("mutation")
reference_combined$mutation <- as.factor(mutation[[1]])

site <- data.frame(sapply(reference_combined[1], fun_site))
colnames(site) <- c("site")
reference_combined$site <- as.numeric(site[[1]])
#Reorder the columns (name, wild type, site, mutation, depth)
reference_combined <- reference_combined[,c(1,3,5,4,2)]

###Order reference_combined by site and mutation
reference_combined <- reference_combined[
  with(reference_combined, order(site, mutation)),
]

#Clean-up
remove(mutation)
remove(site)
remove(wildtype)




###Remove stop codons
escape_trimmed <- escape_combined[which(escape_combined$mutation != "*"),]
reference_trimmed <- reference_combined[which(reference_combined$mutation != "*"),]


###############################################################################
###Add missing residues
###############################################################################
Wuhan <- read.fasta("N_Wuhan.fasta")

###First reference data
### If there are residues without any data, these have to be added manually
#1. Define the range of amino acids present
seq_range <- min(reference_trimmed$site):max(reference_trimmed$site)
missing_aa <- seq_range[!seq_range %in% unique(reference_trimmed$site)]
#2. Add the missing residue(s)
missing_data <- data.frame(site = missing_aa, 
                           name = paste(toupper(Wuhan[["N_Wuhan"]][missing_aa]), missing_aa, toupper(Wuhan[["N_Wuhan"]][missing_aa])),
                           mutation = toupper(Wuhan[["N_Wuhan"]][missing_aa]), 
                           wildtype = toupper(Wuhan[["N_Wuhan"]][missing_aa]))
complete_reference <- rbind.fill(reference_trimmed, missing_data)

#Expand the dataset to include NA values for synonymous amino acids
all <- complete_reference %>% expand(site, mutation)
# join with all, n will be NA for obs. in all that are not present in v
reference_trimmed = complete_reference %>% group_by_at(vars(wildtype, site, mutation)) %>% 
  right_join(all)


###Next escape data
### If there are residues without any data, these have to be added manually
#1. Define the range of amino acids present
seq_range <- min(escape_trimmed$site):max(escape_trimmed$site)
missing_aa <- seq_range[!seq_range %in% unique(escape_trimmed$site)]
#2. Add the missing residue(s)
missing_data <- data.frame(site = missing_aa, 
                           name = paste(toupper(Wuhan[["N_Wuhan"]][missing_aa]), missing_aa, toupper(Wuhan[["N_Wuhan"]][missing_aa])),
                           mutation = toupper(Wuhan[["N_Wuhan"]][missing_aa]), 
                           wildtype = toupper(Wuhan[["N_Wuhan"]][missing_aa]))
complete_escape <- rbind.fill(escape_trimmed, missing_data)

#Expand the dataset to include NA values for synonymous amino acids
all <- complete_escape %>% expand(site, mutation)
# join with all, n will be NA for obs. in all that are not present in v
escape_trimmed = complete_escape %>% group_by_at(vars(wildtype, site, mutation)) %>% 
  right_join(all)



###Order reference_trimmed and escape_trimmed by site and mutation
escape_trimmed <- escape_trimmed[
  with(escape_trimmed, order(site, mutation)),
]
reference_trimmed <- reference_trimmed[
  with(reference_trimmed, order(site, mutation)),
]


###################################################################
### Calculate abundance: n/N
###################################################################
escape_trimmed$abundance <- escape_trimmed$depth/sum(escape_trimmed$depth, na.rm=TRUE)
reference_trimmed$abundance <- reference_trimmed$depth/sum(reference_trimmed$depth, na.rm=TRUE)

###Then adjust the lowest reference abundance values
#Use 95 percentile
#This helps remove exaggerated escape fractions caused by dividing by a very small number
cutoff  <- quantile(reference_trimmed$abundance, probs = .05, na.rm = T)
reference_trimmed$abundance <- ifelse(reference_trimmed$abundance<cutoff, cutoff, reference_trimmed$abundance)

##Fishers exact test
esc_fisher <- escape_trimmed[,1:5]
esc_fisher$ref_depth <- reference_trimmed$depth
esc_sum <- sum(na.omit(esc_fisher$depth))
ref_sum <- sum(na.omit(esc_fisher$ref_depth))

p_values <- data.frame(matrix(ncol = 1, nrow = nrow(esc_fisher)))
colnames(p_values) <- "p"
#Calculate p values for each comparison (using Fisher's exact test)
for (i in 1:nrow(esc_fisher)){
  if (!anyNA(esc_fisher[i,])){
    f_ <- matrix(unlist(c(esc_fisher[i,5], 
                          esc_fisher[i,6],
                          esc_sum-esc_fisher[i,5],
                          ref_sum-esc_fisher[i,6])),2)
    p_ <- fisher.test(f_)
    p_values$p[i] <-p_$p.value
  }
}

p_values$p_adj <- p.adjust(p_values$p, "bonferroni")
p_values$p_adj <- ifelse(p_values$p_adj == 0, 1e-300, p_values$p_adj)

###################################################################
### Calculate escape fractions
### And add adjusted p values
###################################################################
escape_trimmed$Escape <- (escape_trimmed$abundance/reference_trimmed$abundance)
escape_trimmed$p_adj <- p_values$p_adj
max(escape_trimmed$Escape, na.rm=T)
min(escape_trimmed$Escape, na.rm=T)
#Write out raw escape score before normalizing
escape_raw_out <- na.omit(escape_trimmed[,c(1, ncol(escape_trimmed)-2)])
colnames(escape_raw_out) <- c("Mutation", "EscapeScore")
write.csv(escape_raw_out, file = paste(antibody,"_raw_escape.csv", sep = ""), row.names = FALSE)

escape_trimmed$Non_Normalized_Enrich <- escape_trimmed$Escape

###Normalize the data between the 99 and 1 percentiles
max(escape_trimmed$Escape, na.rm=T)
upper_limit  <- quantile(escape_trimmed$Escape, probs = .99, na.rm = T)
lower_limit <- min(escape_trimmed$Escape, na.rm=T)
escape_trimmed$Escape <- ifelse(escape_trimmed$Escape>upper_limit, upper_limit, escape_trimmed$Escape)
escape_trimmed$Escape <- (escape_trimmed$Escape-lower_limit)/(upper_limit-lower_limit)

ggsave(filename = paste(antibody,"_EscapeHistogram_preArcSine.png", sep=""), 
       ggplot(escape_trimmed, aes(x=Escape), title = "Escape Fractions")+
         #geom_histogram(aes(y=..density..), color = "black", fill = "grey", binwidth = 0.025) +
         geom_histogram(color = "black", fill = "grey") +
         xlim(0,quantile(escape_trimmed$Escape, probs = .99, na.rm = T)) +
         xlab("Escape Fraction") + 
         ylab("Count") +
         theme_bw(base_size = 10),
       width = 3, height = 2, dpi = 300, units = "in", device='png')

# ##Calculate destabilization instead of escape!!!
# mx <- max(escape_trimmed$Escape, na.rm=T)
# mean_dest <- median(escape_trimmed$Escape, na.rm=T)
# escape_trimmed$destabilization <-abs(mx-escape_trimmed$Escape)
#escape_trimmed <- na.omit(escape_trimmed)


################################################################################
### Generate a matrix for plotting sequence logos
### This matrix is also used to calculate average escape scores
################################################################################
#remove the extra rwos that only contain the site number and nothing else 
escape_trimmed <- escape_trimmed[which(!is.na(escape_trimmed$mutation)),]

logo_matrix <- matrix(ncol = nrow(escape_trimmed)/20,nrow=20)
row.names(logo_matrix) <- escape_trimmed$mutation[1:20]
colnames(logo_matrix) <- seq(2,ncol(logo_matrix)+1, 1)


for(i in 0:ncol(logo_matrix)-1){
  logo_matrix[,i+1] <- escape_trimmed$Escape[seq(from = 20*i+1, to=20*i+20, by = 1)]
}

################################################################################
##Plot heat maps
################################################################################

#Generate a dataframe containing the wild type amino acids 
#This will be used to mark wild type in the tiled heat map
tmp <- data.frame(sapply(escape_trimmed[1], fun_first), 
                  sapply(escape_trimmed[1], fun_site), 
                  sapply(escape_trimmed[1], fun_first))
colnames(tmp) <- c("wildtype","site", "mutation")
frames <- distinct(tmp)
frames$site <- as.numeric(frames$site)

#Change data type to integer for "site"
#Required for proper heatmap plotting
frames$site <- as.integer(frames$site)
escape_trimmed$site <- as.integer(escape_trimmed$site)

write.table(escape_trimmed,file="differential_matrix.txt")

#Set the order for amino acids in the heatmap (by aa property)
polar <- c("H", "C", "S", "T", "N", "Q")
nonpolar <- c("G", "A", "V", "L", "I", "M", "P")
aromatic <- c("F", "Y", "W")
positive <- c("K", "R")
negative <- c("D", "E")

aa_order <- c(negative,positive, polar, nonpolar, aromatic)

#Residues 1-209
start = 46
end = 175
mut_range <- subset(escape_trimmed, site>start & site<end)
mut_range$site <- as.factor(mut_range$site)
frames_range <- subset(frames, site>start & site<end)
frames_range$site <- as.factor(frames_range$site)

ggsave(filename = paste(antibody,"_EscapeFraction_heatmap01.png", sep=""), 
       ggplot(mut_range, aes(site, mutation, size = -log(p_adj))) + 
         geom_tile(color = "white",
                   fill = '#FCF0F0', #azure1 for blueish grey background
                   lwd = 0.1,
                   linetype = 1,
                   alpha = 1) +
         #scale_y_discrete(limits=rev) +
         ylim(rev(aa_order)) +
         geom_tile(data=frames_range, 
                   size=0,
                   height = 1,
                   fill='white', 
                   colour="white") +
         geom_point(aes(colour = Escape), alpha=1) +
         scale_colour_distiller(palette = "RdPu", direction = +1,
                                na.value = '#FCF0F0',
                                limits=c(0,max(escape_trimmed$Escape, na.rm=T))) +
         scale_size(range = c(0, 2)) +
         #The next section adds lightgrey tiles and black points for 
         #the wild type residues
         geom_point(inherit.aes = FALSE, 
                    data = frames_range,
                    aes(site, mutation),
                    shape=16,
                    size = 1.25,
                    colour="black") +
         xlab("") + ylab("") +
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

################################################################################
### Plot average Escape histogram
################################################################################
average_escape <- as.data.frame(colMeans(logo_matrix, na.rm=T))
average_escape$site <- seq(2,419,1)
colnames(average_escape) <- c("escape", "site")

ggsave(filename = paste(antibody,"_average_escape.png", sep=""), 
       ggplot(average_escape, aes(site, escape)) + 
         geom_bar(stat = "identity") + 
         xlim(2,419) +
         ylim(0, max(average_escape$escape+5)) +
         theme_minimal() + 
         theme_void(),
       width = 12, height = .4, dpi = 300, units = "in", device='png')


#Individual escapes:
escape_fractions_out <- na.omit(escape_trimmed[,c(1,7,8)])
#escape_fractions_out$name <- paste(escape_fractions_out$wildtype, escape_fractions_out$site, escape_fractions_out$mutation, sep = "")
write.csv(escape_fractions_out, file = paste(antibody,"_escape_fractions.csv", sep = ""), row.names = FALSE)

#average escapes:
average_escape_out <- as.data.frame(average_escape$site)
average_escape_out$'average Escape Score' <- average_escape$escape
colnames(average_escape_out) <- c("Site", "average Escape Score")
write.csv(average_escape_out, file = paste(antibody,"_average_escape.csv", sep=""), row.names = FALSE)

################################################################################
#Pymol file for mapping average scores onto the structure
################################################################################
#Escapees
m_ <- mean(average_escape$escape, na.rm=T)
sd_ <- sd(average_escape$escape, na.rm=T)
cut_1 <- (m_ + 1*sd_)
cut_2 <- (m_ + 2*sd_)
escapees_C01 <- subset(average_escape,escape > cut_1 & escape < cut_2)
escapees_C02 <- subset(average_escape,escape > cut_2)
top_escapes <- na.omit(escapees_C02[,c(1,2)])
write.csv(top_escapes, file = paste(antibody,"_top_escapes.csv", sep=""), row.names = FALSE)

average_escape$escape <- ifelse(average_escape$escape < 0, 0, average_escape$escape)
average_escape$level <- ifelse(average_escape$escape < cut_1, "grey", 
                             ifelse(average_escape$escape < cut_2, "red_1", "red_2" ))
average_escape$level <- as.factor(average_escape$level)
#Residues 46-175
start = 46
end = 175
average_escape1 <- subset(average_escape, site>start & site<end)
ggsave(filename = paste(antibody, "_color_average_escape.png", sep=""), 
       ggplot(average_escape, aes(site, escape, fill = level)) + 
         geom_bar(stat="identity") + 
         geom_segment(aes(x=46,xend=175,y=cut_1, yend = cut_1), linetype = 3, size = 0.15) +
         geom_segment(aes(x=46,xend=175,y=cut_2, yend = cut_2), linetype = 3, size = 0.15) +
         #geom_segment(aes(x=215,xend=448,y=m_, yend = m_), linetype = 1, size = 0.03, alpha = 0.5) +
         scale_fill_manual(values = c("grey" = "#d9d9d9",
                                      "red_1" ="#dd3497",
                                      "red_2"="#4a2267")) +
         xlim(46,175) +
         theme_void() +
         theme(legend.position = "none"),
       width = 12, height = 0.4, dpi = 300, units = "in", device='png')

################################################################################
### Write a script for coloring residues in pymol
################################################################################
#Add "resi " in front of epitope residue names
#Redefine escapees_C01 and escapees_C02 to be 2 and 4 sd respectively, unlike 1
# and 2 above. This is to make the pymol structures less cluttered
cut_1 <- (m_ + 1*sd_)
cut_2 <- (m_ + 2*sd_)
escapees_C01 <- subset(average_escape,escape > cut_1 & escape < cut_2)
escapees_C02 <- subset(average_escape,escape > cut_2)

###################
# Need to add 1 to each site so that heatmap aligns with
# structures (since pymol resi numbers from 0)
escapees_C01$site <- escapees_C01$site+1
escapees_C02$site <- escapees_C02$site+1
###################

escapees_C01$resi <- paste("or resi", escapees_C01$site, sep = " ")
reds01 <- paste(escapees_C01$resi)
escapees_C02$resi <- paste("or resi", escapees_C02$site, sep = " ")
reds02 <- paste(escapees_C02$resi)

#RBD
sink(paste(ab_name, "_pymol_RBD.txt", sep = ""))
cat("#Set custom colors with RGB:")
cat("\n")
cat("set ray_shadow, 0")
cat("\n")
cat("color white, ChainD")
cat("\n")
cat("set ray_trace_mode, 1")
cat("\n")
cat("set_color blue_1, [192,192,192]")
cat("\n")
cat("set_color grey, [70,70,70]")
cat("\n")
cat("color blue_1, ChainD")
cat("\n")
cat("set_color red_1, [221,52,151]")
cat("\n")
cat("set_color red_2, [73,0,106]")
cat("\n")
cat("#Set sphere radius:")
cat("\n")
cat("set sphere_scale, 0.5")
cat("\n")
cat("#Select sites to highlight with spheres:")
cat("\n")
cat("create spheres01_", ab_name, ", ChainD and (resi 0 ", sep = "")
cat(paste(reds01, sep = ""), ")")
cat("\n")
cat("set sphere_color, red_1, spheres01_", ab_name, sep = "")
cat("\n")
cat("hide everything, spheres01_", ab_name, sep = "")
cat("\n")
cat("show spheres, spheres01_", ab_name, " and name CA", sep = "")
cat("\n")
cat("create spheres02_", ab_name, ", ChainD and (resi 0 ", sep = "")
cat(paste(reds02, sep = ""), ")")
cat("\n")
cat("hide everything, spheres02_", ab_name, sep = "")
cat("\n")
cat("show spheres, spheres02_", ab_name, " and name CA", sep = "")
cat("\n")
cat("set sphere_color, red_2, spheres02_", ab_name, sep = "")
cat("\n")
cat("create spheres03_CoreBottom, ChainD and (resi 51 or resi 52 or resi 53 or resi 71 or resi 72 or resi 73 or resi 84 or resi 85 or resi 86 or resi 87 or resi 110 or resi 111 or resi 112 or resi 113 or resi 114 or resi 115 or resi 116 or resi 130 or resi 132)
set sphere_color, grey, spheres03_CoreBottom")
cat("\n")
cat("set sphere_color, grey, spheres03_CoreBottom")
cat("\n")
cat("hide everything, spheres03_CoreBottom")
cat("\n")
cat("show spheres, spheres03_CoreBottom and name CA")
cat("\n")
cat("create surface_", ab_name, ", ChainD", sep = "")
cat("\n")
cat("select surface01_", ab_name, ", surface_", ab_name,  " and (resi 0 ", sep = "")
cat(paste(reds01, sep = ""), ")")
cat("\n")
cat("select surface02_", ab_name, ", surface_", ab_name,  " and (resi 0 ", sep = "")
cat(paste(reds02, sep = ""), ")")
cat("\n")
cat("select surface03_CoreBottom, surface_CoreBottom and (resi 51 or resi 52 or resi 53 or resi 71 or resi 72 or resi 73 or resi 84 or resi 85 or resi 86 or resi 87 or resi 110 or resi 111 or resi 112 or resi 113 or resi 114 or resi 115 or resi 116 or resi 130 or resi 132)")
cat("\n")
cat("color red_1, surface01_", ab_name, sep = "")
cat("\n")
cat("color red_2, surface02_", ab_name, sep = "")
cat("\n")
cat("color grey, surface03_CoreBottom")
cat("\n")
cat("show surface, surface_", ab_name, sep = "")
cat("\n")
sink()


