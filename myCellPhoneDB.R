callCellPhoneDB <- function(SeuratObject,
                            savePath,
                            pvalues.threshold = 0.01){
    
    sampleName <- SeuratObject@project.name
    
    savePath <- filePathCheck(savePath)
    
    savePath <- R.utils::getAbsolutePath(savePath)
    
    savePath_basic <- file.path(savePath, sampleName)
    savePath_interact <- file.path(savePath_basic, "interact")
    
    labels <- list.files(savePath_interact)
    files <- file.path(savePath_interact, labels)
    
    cellphone.db <- readRDS("./cellphonedb.RDS")
    
    for(i in 1 : length(files)){
        meta <- read.table(file.path(files[[i]], "meta.txt"), 
                           sep = "\t", 
                           header = T)
        counts <- read.table(file.path(files[[i]], "counts.txt"), 
                             sep = "\t", 
                             header = T, 
                             row.names = 1)
        myCellPhoneDB(meta, 
                      counts, 
                      files[[i]], 
                      cellphone.db, 
                      pvalues.threshold)
    }
}



myCellPhoneDB <- function(meta, 
                          counts, 
                          savePath, 
                          cellphone.db, 
                          pvalues.threshold = 0.01){
    
    ### transform mRNA counts to protein counts
    # read all tables
    complex_table <- cellphone.db[["complex_table"]]
    complex_composition_table <- cellphone.db[["complex_composition_table"]]
    gene_table <- cellphone.db[["gene_table"]]
    interaction_table <- cellphone.db[["interaction_table"]]
    multidata_table <- cellphone.db[["multidata_table"]]
    protein_table <- cellphone.db[["protein_table"]]
    
    # find common genes
    common_gene <- intersect(rownames(counts), gene_table$ensembl)
    
    # filter data and tables
    counts <- counts[common_gene, ]
    gene_table <- gene_table %>%
      subset(ensembl %in% common_gene)
    
    ## get the gene and protein table, and find out all the possible protein
    gene_protein_table <- merge(gene_table, 
                                protein_table, 
                                by.x = "protein_id", 
                                by.y = "id_protein")
    # find single protein
    protein_multidata_id <- gene_protein_table$protein_multidata_id
    
    # find complex
    filtered_complex_composition_table <- merge(complex_composition_table, 
                                                gene_protein_table, 
                                                by.x = "protein_multidata_id", 
                                                by.y = "protein_multidata_id")
    complex_multidata_id <- unique(filtered_complex_composition_table$complex_multidata_id)
    for(i in 1 : length(complex_multidata_id)){
        i_id <- complex_multidata_id[i]
        i_table <- filtered_complex_composition_table %>%
          subset(complex_multidata_id == i_id)
        if(nrow(i_table) != i_table$total_protein){
            complex_multidata_id[i] <- NA
        }
    }
    complex_multidata_id <- na.omit(complex_multidata_id)
    
    # combine single protein and complex
    multidata_id <- c(protein_multidata_id, complex_multidata_id)
    
    ## build protein counts
    protein_counts <- matrix(nrow = length(multidata_id), 
                             ncol = ncol(counts))
    rownames(protein_counts) <- as.character(multidata_id)
    colnames(protein_counts) <- colnames(counts)
    
    rownames(gene_protein_table) <- gene_protein_table$protein_multidata_id
    gene_protein_table <- gene_protein_table[as.character(protein_multidata_id), ]
    
    for(i in 1 : length(protein_multidata_id)){
        i_ensembl <- gene_protein_table$ensembl[i]
        protein_counts[i, ] <- as.matrix(counts[i_ensembl, ])
    }
    
    ## build complex counts
    for(i in 1 : length(complex_multidata_id)){
        i_table <- complex_composition_table[
          complex_composition_table$complex_multidata_id %in% 
            complex_multidata_id[i], ]
        i_protein_ids <- i_table$protein_multidata_id
        i_gene_protein_table <- gene_protein_table[as.character(i_protein_ids), ]
        i_ensembl <- i_gene_protein_table$ensembl
        protein_counts[as.character(complex_multidata_id[i]), ] <- 
          cal_complex(counts[i_ensembl, ])
        protein_counts[as.character(complex_multidata_id[i]), ] <- 
          cal_complex(counts[i_ensembl, ])
    }
    
    # find possible interaction
    filtered_interaction_table <- interaction_table %>% 
      subset(multidata_1_id %in% multidata_id) %>%
      subset(multidata_2_id %in% multidata_id)
    
    ## statistic method
    statistic_results <- statistic_method(meta, 
                                          protein_counts, 
                                          filtered_interaction_table)
    
    ## get pvalue
    clusters <- unique(meta[, 2])
    cluster_1 <- clusters[1]
    cluster_2 <- clusters[2]
    n_interaction <- nrow(filtered_interaction_table)
    
    protein_counts_1 <- protein_counts[, make.names(meta[meta[, 2] == cluster_1, ][, 1])]
    protein_counts_2 <- protein_counts[, make.names(meta[meta[, 2] == cluster_2, ][, 1])]
    mean_protein_counts_1 <- rowMeans(protein_counts_1)
    mean_protein_counts_2 <- rowMeans(protein_counts_2)
    
    # build real matrix
    real_mean <- matrix(nrow = n_interaction, 
                           ncol = 4)
    colnames(real_mean) <- c(paste0(cluster_1, "|", cluster_1), 
                             paste0(cluster_1, "|", cluster_2), 
                             paste0(cluster_2, "|", cluster_1), 
                             paste0(cluster_2, "|", cluster_2))
    rownames(real_mean) <- filtered_interaction_table$id_cp_interaction
    
    
    for(i in 1 : n_interaction){
        real_mean[i, paste0(cluster_1, "|", cluster_1)] <- 
          cal_mean(mean_protein_counts_1[as.character(filtered_interaction_table[i, ]$multidata_1_id)], 
                   mean_protein_counts_1[as.character(filtered_interaction_table[i, ]$multidata_2_id)])
        real_mean[i, paste0(cluster_1, "|", cluster_2)] <- 
          cal_mean(mean_protein_counts_1[as.character(filtered_interaction_table[i, ]$multidata_1_id)], 
                   mean_protein_counts_2[as.character(filtered_interaction_table[i, ]$multidata_2_id)])
        real_mean[i, paste0(cluster_2, "|", cluster_1)] <- 
          cal_mean(mean_protein_counts_2[as.character(filtered_interaction_table[i, ]$multidata_1_id)], 
                   mean_protein_counts_1[as.character(filtered_interaction_table[i, ]$multidata_2_id)])
        real_mean[i, paste0(cluster_2, "|", cluster_2)] <- 
          cal_mean(mean_protein_counts_2[as.character(filtered_interaction_table[i, ]$multidata_1_id)], 
                   mean_protein_counts_2[as.character(filtered_interaction_table[i, ]$multidata_2_id)])
    }
    
    # pvalue matrix
    pvalue_matrix <- matrix(nrow = n_interaction, 
                            ncol = 4)
    colnames(pvalue_matrix) <- c(paste0(cluster_1, "|", cluster_1), 
                                 paste0(cluster_1, "|", cluster_2), 
                                 paste0(cluster_2, "|", cluster_1), 
                                 paste0(cluster_2, "|", cluster_2))
    rownames(pvalue_matrix) <- filtered_interaction_table$id_cp_interaction
    
    for(i in 1 : 4){
        for(j in 1 : n_interaction){
            n_bigger <- 0
            for(iter in 1 : 1000){
                if(statistic_results[[i]][j, iter] > real_mean[j, i]){
                    n_bigger <- n_bigger + 1
                }
            }
            pvalue_matrix[j, i] <- n_bigger / 1000
        }
    }
    
    real_mean <- real_mean[which(rowSums(real_mean) > 0), ]
    pvalue_matrix <- pvalue_matrix[rownames(real_mean), ]
    
    # select interactions whose pvalue is less than pvalues.threshold
    pvalue_matrix <- pvalue_matrix[rowSums(pvalue_matrix <= pvalues.threshold) != 0, ]
    valid_interaction <- rownames(pvalue_matrix)
    
    mean_matrix <- real_mean[valid_interaction, ]
    pvalue_matrix <- pvalue_matrix[valid_interaction, ]
    
    rownames(filtered_interaction_table) <- filtered_interaction_table$id_cp_interaction
    filtered_interaction_table <- filtered_interaction_table[valid_interaction, ]
    
    ## get table base
    table_base <- build_table_base(filtered_interaction_table = filtered_interaction_table, 
                                   multidata_table = multidata_table, 
                                   gene_protein_table = gene_protein_table)
    
    ## means.txt
    means.txt <- cbind(table_base, 
                       as.data.frame(mean_matrix))
    pvalues.txt <- cbind(table_base, 
                        as.data.frame(pvalue_matrix))
    
    savePath <- file.path(savePath, "out")
    dir.create(savePath, recursive = T)
    
    write.table(means.txt,
                file.path(savePath,
                          "means.txt"),
                row.names = F,
                quote = F,
                sep = "\t")
    
    write.table(pvalues.txt,
                file.path(savePath,
                          "pvalues.txt"),
                row.names = F,
                quote = F,
                sep = "\t")
}


statistic_method <- function(meta,
                             counts,
                             interaction_table, 
                             iteration = 1000){
    
    clusters <- unique(meta[, 2])
    
    if(length(clusters) != 2){
        stop("The number of types of meta is not 2!")
    }
    
    cluster_1 <- clusters[1]
    cluster_2 <- clusters[2]
    
    n_cluster <- nrow(meta)
    n_cluster_1 <- sum(meta[, 2] == cluster_1)
    n_cluster_2 <- sum(meta[, 2] == cluster_2)
    n_interaction <- nrow(interaction_table)
    
    # "results" includes 4 matrices, and the name of matrix stands for "ligand_receptor"
    results <- list()
    results[[paste0(cluster_1, "|", cluster_1)]] <- matrix(nrow = n_interaction,
                                                           ncol = iteration)
    results[[paste0(cluster_1, "|", cluster_2)]] <- matrix(nrow = n_interaction,
                                                           ncol = iteration)
    results[[paste0(cluster_2, "|", cluster_1)]] <- matrix(nrow = n_interaction,
                                                           ncol = iteration)
    results[[paste0(cluster_2, "|", cluster_2)]] <- matrix(nrow = n_interaction,
                                                           ncol = iteration)
    
    for(i in 1 : iteration){
        rank <- sample(c(1:n_cluster), n_cluster)
        i_counts_cluster_1 <- counts[, rank[1:n_cluster_1]]
        i_counts_cluster_2 <- counts[, rank[(n_cluster_1+1):n_cluster]]
        
        for(j in 1 : nrow(interaction_table)){
            j_table <- interaction_table[j, ]
            interaction_1_l <- mean(i_counts_cluster_1[as.character(j_table$multidata_1_id), ])
            interaction_1_r <- mean(i_counts_cluster_1[as.character(j_table$multidata_2_id), ])
            interaction_2_l <- mean(i_counts_cluster_2[as.character(j_table$multidata_1_id), ])
            interaction_2_r <- mean(i_counts_cluster_2[as.character(j_table$multidata_2_id), ])
            
            results[[paste0(cluster_1, "|", cluster_1)]][j, i] <- cal_mean(interaction_1_l, interaction_1_r)
            results[[paste0(cluster_1, "|", cluster_2)]][j, i] <- cal_mean(interaction_1_l, interaction_2_r)
            results[[paste0(cluster_2, "|", cluster_1)]][j, i] <- cal_mean(interaction_2_l, interaction_1_r)
            results[[paste0(cluster_2, "|", cluster_2)]][j, i] <- cal_mean(interaction_2_l, interaction_2_r)
        }
    }
    
    return(results)
}

cal_mean <- function(mean_1, 
                     mean_2){
    if(mean_1 == 0 | mean_2 == 0){
        return(0)
    }else{
        return((mean_1 + mean_2) / 2)
    }
}

cal_complex <- function(counts){
    if(sum(counts == 0) != 0){
        return(0)
    }else{
        # return(colMeans(counts))
        return(apply(counts, 2, min))
    }
}


build_table_base <- function(filtered_interaction_table, 
                             multidata_table, 
                             gene_protein_table){
    ## id_cp_interaction
    id_cp_interaction <- filtered_interaction_table$id_cp_interaction
    ## interaction_pair, ensembl, partner, secreted and is_integrin
    ligand_id <- filtered_interaction_table$multidata_1_id
    receptor_id <- filtered_interaction_table$multidata_2_id
    interaction_pair <- list()
    partner_a <- list()
    partner_b <- list()
    ensembl_a <- list()
    ensembl_b <- list()
    secreted <- list()
    receptor_a <- list()
    receptor_b <- list()
    is_integrin <- list()
    
    for(i in 1 : length(ligand_id)){
        # secreted flag
        flag_secreted <- FALSE
        # is_integrin flag
        flag_integrin <- FALSE
        # ligand
        i_table <- subset(gene_protein_table, protein_multidata_id == as.character(ligand_id[i]))
        # complex
        if(nrow(i_table) != 1){
            i_table <- subset(multidata_table, id_multidata == as.character(ligand_id[i]))
            i_ligand <- i_table$name
            partner_a[[i]] <- paste0("complex:", i_table$name)
            ensembl_a[[i]] <- ""
            receptor_a[[i]] <- i_table$receptor == 1
            
            if(i_table$secreted){
                flag_secreted <- TRUE
            }
            flag_integrin <- TRUE
        }else{ # single protein
            i_ligand <- i_table$gene_name
            ensembl_a[[i]] <- i_table$ensembl
            i_table <- subset(multidata_table, id_multidata == as.character(ligand_id[i]))
            partner_a[[i]] <- paste0("single:", i_table$name)
            receptor_a[[i]] <- i_table$receptor == 1
            
            if(i_table$secreted){
                flag_secreted <- TRUE
            }
        }
        
        # receptor
        i_table <- subset(gene_protein_table, protein_multidata_id == as.character(receptor_id[i]))
        # complex
        if(nrow(i_table) != 1){
            i_table <- subset(multidata_table, id_multidata == as.character(receptor_id[i]))
            i_receptor <- i_table$name
            partner_b[[i]] <- paste0("complex:", i_table$name)
            ensembl_b[[i]] <- ""
            receptor_b[[i]] <- i_table$receptor == 1
            
            if(i_table$secreted){
                flag_secreted <- TRUE
            }
            flag_integrin <- TRUE
        }else{ # single protein
            i_receptor <- i_table$gene_name
            ensembl_b[[i]] <- i_table$ensembl
            i_table <- subset(multidata_table, id_multidata == as.character(receptor_id[i]))
            partner_b[[i]] <- paste0("single:", i_table$name)
            receptor_b[[i]] <- i_table$receptor == 1
            
            if(i_table$secreted){
                flag_secreted <- TRUE
            }
        }
        
        # interaction_pair
        interaction_pair[i] <- paste0(i_ligand, "_", i_receptor)
        # secreted
        secreted[[i]] <- flag_secreted
        # is_integrin
        is_integrin[[i]] <- flag_integrin
    }
    
    return(data.frame(id_cp_interaction = unlist(id_cp_interaction), 
                      interaction_pair = unlist(interaction_pair), 
                      partner_a = unlist(partner_a), 
                      partner_b = unlist(partner_b), 
                      gene_a = unlist(ensembl_a), 
                      gene_b = unlist(ensembl_b), 
                      secreted = unlist(secreted), 
                      receptor_a = unlist(receptor_a), 
                      receptor_b = unlist(receptor_b), 
                      annotation_strategy = filtered_interaction_table$annotation_strategy, 
                      is_integrin = unlist(is_integrin), 
                      check.names = FALSE))
}
