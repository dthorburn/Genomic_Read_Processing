## Parse an unphased VCF file for it's annotationsto generate a able/plots to inform hard filtering in GATK

suppressMessages(library(vcfR))
suppressMessages(library(dplyr))
suppressMessages(library(ggpubr))
suppressMessages(library(ggplot2))
suppressMessages(library(stringr))
suppressMessages(library(data.table))

Genome <- commandArgs(trailingOnly=TRUE)[1]
fil_name <- Genome %>% gsub(pattern = ".*/", "", .) %>% gsub(pattern = "\\..*", "", .)
temp_file <- read.vcfR(Genome, verbose = FALSE)

Ann_info <- temp_file@fix[,"INFO"] 

output <- data.table( 	File = fil_name,
						QUAL = temp_file@fix[,"QUAL"] %>% as.numeric(),
						QD 	= str_extract(pattern = "QD=[0-9]+\\.[0-9]+?", Ann_info) %>% gsub(pattern = "QD=", "", .) %>% as.numeric(),
						FS 	= str_extract(pattern = "FS=[0-9]+\\.[0-9]+?", Ann_info) %>% gsub(pattern = "FS=", "", .) %>% as.numeric(),
						MQ 	= str_extract(pattern = "MQ=[0-9]+\\.[0-9]+?", Ann_info) %>% gsub(pattern = "MQ=", "", .) %>% as.numeric(),
						SOR = str_extract(pattern = "SOR=[0-9]+\\.[0-9]+?", Ann_info) %>% gsub(pattern = "SOR=", "", .) %>% as.numeric())

#fwrite(file = paste0(Output_path, "/", fil_name, "_Annotation_Distribution.csv"), output)


p1 <- ggplot(output, aes(x = QD, fill = File, colour = File)) +
			geom_density(alpha = 0.2, size = 1) +
			theme_minimal(base_size = 20) +
			labs(x = "QD", y = "Density") +
			ggtitle(fil_name) +
			geom_vline(xintercept = mean(output$QD, na.rm = TRUE), size = 1, colour = "black", linetype = "dashed") +
			geom_vline(xintercept = median(output$QD, na.rm = TRUE), size = 1, colour = "black", linetype = "solid") +
			theme(legend.position = "none")

p2 <- ggplot(output, aes(x = FS, fill = File, colour = File)) +
			geom_density(alpha = 0.2, size = 1, colour = "blue") +
			theme_minimal(base_size = 20) +
			labs(x = paste0("FS (range:",min(output$FS),"-", max(output$FS), ")"), y = "Density") +
			geom_vline(xintercept = mean(output$FS, na.rm = TRUE), size = 1, colour = "black", linetype = "dashed") +
			geom_vline(xintercept = median(output$FS, na.rm = TRUE), size = 1, colour = "black", linetype = "solid") +
			theme(legend.position = "none") +
			scale_fill_manual(values = "blue") +
			xlim(0,10) 

p3 <- ggplot(output, aes(x = MQ, fill = File, colour = File)) +
			geom_density(alpha = 0.2, size = 1, colour = "green") +
			theme_minimal(base_size = 20) +
			labs(x = "MQ", y = "Density") +
			geom_vline(xintercept = mean(output$MQ, na.rm = TRUE), size = 1, colour = "black", linetype = "dashed") +
			geom_vline(xintercept = median(output$MQ, na.rm = TRUE), size = 1, colour = "black", linetype = "solid") +
			theme(legend.position = "none") +
			scale_fill_manual(values = "green")

p4 <- ggplot(output, aes(x = SOR, fill = File, colour = File)) +
			geom_density(alpha = 0.2, size = 1, colour = "purple") +
			theme_minimal(base_size = 20) +
			labs(x = "SOR", y = "Density") +
			geom_vline(xintercept = mean(output$SOR, na.rm = TRUE), size = 1, colour = "black", linetype = "dashed") +
			geom_vline(xintercept = median(output$SOR, na.rm = TRUE), size = 1, colour = "black", linetype = "solid") +
			theme(legend.position = "none") +
			scale_fill_manual(values = "purple")

figure <- ggarrange(p1, p2, p3, p4, ncol = 2, nrow = 2)
figure
ggsave(file = paste0("./06_SelectVariants/", fil_name, "_Annotation_Distribution.png"), height = 300, width = 400, units = "mm")
