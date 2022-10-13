# Dealing with the vcf files.
# Building background database
normal_list <- c('Ts11-P2N1','Ts11-P2N5','Ts15-P2N1','Ts15-P2N1','Ts6+8-P2N','Ts6+8-P2N7','Ts6+8-P2N11','Ts15+8-P2N2','Ts15+8-P2N5',
                 'Ts15+8-P2N11','Ts15+8-P2N12','Ts15+8-P2N13')
ident1 <- c()
for (i in normal_list) {
  temp1 <- read.table(paste0("./temp_data/",i,".vcf"), sep = '\t')
  ident <- paste(temp1$V1,temp1$V2,temp1$V4,temp1$V5,sep = '_')
  ident1 <- c(ident1, ident)
}
# Filtering the target mutations
tumor_list <- c('Ts11-P2M1','Ts11-P2M5','Ts15-P2M1','Ts15-P2M1-2','Ts6+8-P2M','Ts6+8-P2M7','Ts6+8-P2M11','Ts15+8-P2M2','Ts15+8-P2M2',
                'Ts15+8-P2M5','Ts15+8-P2M12','Ts15+8-P2M13')
for (i in tumor_list) {
  temp1 <- read.table(paste0("./temp_data/",i,".vcf"), sep = '\t')
  ident <- paste(temp1$V1,temp1$V2,temp1$V4,temp1$V5,sep = '_')
  temp1 <- temp1[!ident%in%ident1,]
  print(i)
  print(dim(temp1)[1]) 
  write.table(temp1, paste0("./temp_data/",i,"_filtered.vcf"), sep = '\t', row.names = F, quote = F, col.names = F)
}

# Filtering the target mutations
tumor_list <- c('WES_Ts11-1-M2P',
                'WES_Ts11-5-M2P',
                'WES_Ts15-1-M2P-1',
                'WES_Ts6-3-M2P',
                'WES_Ts8+15-2-M2P',
                'WES_Ts8-2-M2P',
                'WES_Ts6-4-M2P-2')
for (i in tumor_list) {
  temp1 <- read.table(paste0("./temp_data/",i,".vcf"), sep = '\t')
  temp1 <- temp1[temp1$V1!="MT",]
  ident <- paste(temp1$V1,temp1$V2,temp1$V4,temp1$V5,sep = '_')
  temp1 <- temp1[!ident%in%ident1,]
  temp1 <- temp1[temp1$V7=='PASS',]
  print(i)
  print(dim(temp1)[1]) 
  write.table(temp1, paste0("E:./temp_data/",i,"_filtered.vcf"), sep = '\t', row.names = F, quote = F, col.names = F)
  ident <- paste(temp1$V1,temp1$V2,temp1$V4,temp1$V5,sep = '_')
  ident1 <- c(ident1, ident)
}
