load('C:\\Users\\budel\\Desktop\\Multiple_SNP_Multiple_traits\\chap_4_real_data_analysis\\deliang_bu\\real_data_processed\\Tanaka\\data_matrix_final.Rdata')
load('C:\\Users\\budel\\Desktop\\Multiple_SNP_Multiple_traits\\chap_4_real_data_analysis\\deliang_bu\\real_data_processed\\Tanaka\\final_1kg_genotype_correct.Rdata')
load('C:\\Users\\budel\\Desktop\\Multiple_SNP_Multiple_traits\\chap_4_real_data_analysis\\deliang_bu\\real_data_processed\\Tanaka\\covariance_matrix_estimated.Rdata')
load('C:\\Users\\budel\\Desktop\\Multiple_SNP_Multiple_traits\\chap_4_real_data_analysis\\deliang_bu\\real_data_processed\\Tanaka\\gene_list.Rdata')
load('C:\\Users\\budel\\Desktop\\Multiple_SNP_Multiple_traits\\chap_4_real_data_analysis\\deliang_bu\\real_data_processed\\Tanaka\\gene_length_list.Rdata')
PCA_data<-read.table('C:\\Users\\budel\\Desktop\\Multiple_SNP_Multiple_traits\\chap_4_real_data_analysis\\deliang_bu\\real_data_processed\\Tanaka\\snpinfo_12.txt',head=FALSE)
library(MSKAT)
library(AssocTests)

used_PCA_SNP<-intersect(PCA_data$V1, final_1kg_genotype_correct$SNP)
PCA_data_1000_genome_index<-match(used_PCA_SNP,final_1kg_genotype_correct$SNP)
PCA_data_1000_genome<-final_1kg_genotype_correct[PCA_data_1000_genome_index,]
PCA_data_clean<-PCA_data_1000_genome[,7:2510]
write.table(PCA_data_clean, file = "pcocG.eg.txt", quote = FALSE,sep = "", row.names = FALSE, col.names = FALSE)
PCA_result<-pcoc("pcocG.eg.txt",outFile.txt = "pcoc.result.txt",n.MonteCarlo = 1000,num.splits = 10,miss.val = 9)
covariates<-PCA_result$principal.coordinates[,1:2]

n_sim=1000000
q=6
gamma_value<-1
phenotype_covariance<-cor(covariance_matrix_data[8:13])
select_gene<-gene_list[9592]
match_char=paste(select_gene,"\\(",sep = "")
match_line=which(grepl(match_char,data_matrix_final$ANNOT))
genotype_data<-t(final_1kg_genotype_correct[match_line,7:2510])
adj_geno<-lm(genotype_data~covariates)$residuals
genotype_covariance<-cor(adj_geno)
n_sample<-nrow(genotype_data)
geno_length<-ncol(genotype_data)
cutoff_value<-seq(from=0.1,to=1,by=0.1)
c<-as.matrix(rnorm(2504,0,1))
gamma_vec<-as.matrix(rep(gamma_value,q))
c_mat<-c%*%t(gamma_vec)

pheno_rho<-0.8
covariance_matrix<- true_cov_mat_generation_autoregressive(q, pheno_rho)
covariance_matrix<-cor(covariance_matrix_data[8:13])

B_mat<-matrix(0, nrow = geno_length, ncol = q)

p_value_MGAS<-rep(NA,n_sim)
p_value_TWT<-rep(NA,n_sim)
p_value_metaCCA<-rep(NA,n_sim)
p_value_MAT<-rep(NA,n_sim)
p_1<-rep(NA,n_sim)
p_2<-rep(NA,n_sim)
p_3<-rep(NA,n_sim)


est_pheno_rho_mat<-covariance_matrix
total_mat<-kronecker(genotype_covariance,est_pheno_rho_mat)
null_distribution_T3<-generate_null_distribution_T3(q*geno_length,10000,total_mat,cutoff_value)
coefficient_matrix_T3<-approximate_distribution_coefficient_estimate_T3(null_distribution_T3)
for (i in 1:n_sim) {
  eplison<- mvrnorm(n=n_sample,rep(0,q),covariance_matrix)
  Y_mat<-genotype_data%*%B_mat+c_mat+eplison
  #for (j in 1:ncol(independent_snp_mat)) {
  #  for (k in 1:q) {
  #    model<-lm(Y_mat[,k]~independent_snp_mat[,j])
  #    wald <- coef(model)[2] / sqrt((diag(vcov(model)))[2])
  #    wald_mat_independent[j,k]<-wald
  #  }
  #}
  #est_pheno_rho_mat<-cor(wald_mat_independent)
  Z_mat<-matrix(NA,ncol(genotype_data),ncol(Y_mat))
  std_beta_mat<-matrix(NA,ncol(genotype_data),ncol(Y_mat))
  for (j in 1:ncol(genotype_data)) {
    for (k in 1:ncol(Y_mat)) {
      model<- lm(Y_mat[,k]~genotype_data[,j]+covariates+c)
      wald <- coef(model)[2] / sqrt(diag(vcov(model)))[2]
      Z_mat[j,k] <- wald
      std_beta_mat[j,k]<- coef(model)[2]/sqrt(n_sample*(diag(vcov(model))[2]))
    }
  }
  z_vector<-as.vector(t(Z_mat))
  p_value_MGAS[i]<-MGAS(z_vector,genotype_covariance,est_pheno_rho_mat)
  cc<-new_method(Z_mat,genotype_covariance,est_pheno_rho_mat,cutoff_value,coefficient_matrix_T3)
  p_value_TWT[i]<-cc$p_value_final
  p_1[i]<-cc$p_1
  p_2[i]<-cc$p_2
  p_3[i]<-cc$p_3
  p_value_metaCCA[i]<-metaCCA(genotype_covariance,est_pheno_rho_mat,std_beta_mat,n_sample)
  p_value_MAT[i]<-MSATS(Z_mat,est_pheno_rho_mat,genotype_covariance)$p.value[1]
  print(i)
}

sum(as.numeric(p_value_MGAS<0.00001),na.rm = T)
sum(as.numeric(p_value_TWT<0.00001),na.rm = T)
sum(as.numeric(p_value_metaCCA<0.00001),na.rm = T)
sum(as.numeric(p_value_MAT<0.00001),na.rm = T)
sum(as.numeric(p_1<0.00001),na.rm = T)
sum(as.numeric(p_2<0.00001),na.rm = T)
sum(as.numeric(p_3<0.00001),na.rm = T)

sum(as.numeric(p_value_MGAS<0.0001),na.rm = T)
sum(as.numeric(p_value_TWT<0.0001),na.rm = T)
sum(as.numeric(p_value_metaCCA<0.0001),na.rm = T)
sum(as.numeric(p_value_MAT<0.0001),na.rm = T)
sum(as.numeric(p_1<0.0001),na.rm = T)
sum(as.numeric(p_2<0.0001),na.rm = T)
sum(as.numeric(p_3<0.0001),na.rm = T)

sum(as.numeric(p_value_MGAS<0.001),na.rm = T)
sum(as.numeric(p_value_TWT<0.001),na.rm = T)
sum(as.numeric(p_value_metaCCA<0.001),na.rm = T)
sum(as.numeric(p_value_MAT<0.001),na.rm = T)
sum(as.numeric(p_1<0.001),na.rm = T)
sum(as.numeric(p_2<0.001),na.rm = T)
sum(as.numeric(p_3<0.001),na.rm = T)

sum(as.numeric(p_value_MGAS<0.01),na.rm = T)
sum(as.numeric(p_value_TWT<0.01),na.rm = T)
sum(as.numeric(p_value_metaCCA<0.01),na.rm = T)
sum(as.numeric(p_value_MAT<0.01),na.rm = T)
sum(as.numeric(p_1<0.01),na.rm = T)
sum(as.numeric(p_2<0.01),na.rm = T)
sum(as.numeric(p_3<0.01),na.rm = T)