args <- commandArgs(trailingOnly = TRUE)
data_path <- args[1]

getCutSite <- function(cigar,start_pos, ref_pos = 98){
  
  len <- unlist(strsplit(cigar, "[MXSID]"))
  len <- as.numeric(len[len != ""])
  type <- unlist(strsplit(cigar, "[0-9]+"))
  type <- type[type != ""]
  pos_ref <- start_pos - 1
  pos_seq <- 0
  for(i in 1 : length(len)){
    if(i != 1 & type[[i]] == "S"){
      break
    }
    if(type[[i]] == "D" | type[[i]] == "M" | type[[i]] == "X"){
      pos_ref = pos_ref + len[i]
    }
    if(i == 1 & type[[i]] == "S"){
      pos_seq = pos_seq + len[i]
    } else if(type[[i]] == "S"){
      break
    }
    if(type[[i]] == "I" | type[[i]] == "M" | type[[i]] == "X"){
      pos_seq = pos_seq + len[i]
    }
    if(pos_ref == ref_pos){
      break
    }
    if(pos_ref > ref_pos){
      type <- type[[i]]
      if(type == "M" | type == "X"){
        pos_seq <- pos_seq - (pos_ref - ref_pos)
      }
      if(type == "I"){
        pos_seq <- pos_seq - len[i]
      }
      break
    }
  }
  pos_seq
}
fixSeq2 <- function(seq, cigar, start_pos, reference, apos, delLen = 6){
  fixseq <- fixSeq(seq, cigar, start_pos)
  index <- 1
  step <- 1
  ref <- fixseq$ref
  seq <- fixseq$seq
  while(index < apos){
    if(step > length(ref)){
      return(fixseq)
    }
    if(str_detect(ref[step], "[A-Z]")){
      index <- index + 1
    }
    step <- step + 1
  }
  step_apos <- step
  ref_pos <- apos
  delete_end <- 0
  #找到delete的第一个碱基
  while(step > 1){
    step <- step - 1
    if(ref[step] != "-"){
      ref_pos <- ref_pos - 1
    }
    if(ref[step] == "N"){
      delete_end <- step
      break
    } 
  }
  if(delete_end < (step_apos - delLen) | delete_end <= start_pos + 1 | 
     delete_end == (step_apos - 1) | delete_end == 0){
    return(fixseq)
  }
  delete_start <- 0
  while(step > 1){
    step <- step - 1
    if(ref[step] != "N"){
      delete_start = step + 1
      break
    }
    if(ref[step] != "-"){
      ref_pos <- ref_pos - 1
    }
  }
  if(delete_start == 0){
    return(fixseq)
  }
  reference <- unlist(strsplit(reference, "*"))
  while(delete_end < (step_apos - 1)){
    delete_end <- delete_end + 1
    while(ref_pos < delete_end){
      if(reference[ref_pos] == ref[delete_end] | reference[ref_pos] == "N"){
        tmp <- ref[delete_start]
        ref[delete_start] <- ref[delete_end]
        ref[delete_end] <- tmp
        tmp <- seq[delete_start]
        seq[delete_start] <- seq[delete_end]
        seq[delete_end] <- tmp
        break 
      } else {
        ref_pos <- ref_pos + 1
        delete_start <- delete_start + 1
      }
    }
    
  }
  fixseq$ref <- ref
  fixseq$seq <- seq
  fixseq
  
}
#这个函数的目的是为了获得ref和seq但是是经过利用cigar填充的
fixSeq <- function(seq, cigar, start_pos){
  len <- unlist(strsplit(cigar, "[MXSID]"))
  len <- as.numeric(len[len != ""])
  type <- unlist(strsplit(cigar, "[0-9]+"))
  type <- type[type != ""]
  ref_seq <- c()
  map_seq <- c()
  seq <- unlist(strsplit(seq, "*"))
  if(start_pos != 1){
    map_seq <- ref_seq <- rep("N", start_pos - 1)
  }
  pos_seq <- 1
  for(i in 1 : length(len)){
    if(type[[i]] == "M" | type[[i]] == "X"){
      ref_seq <- c(ref_seq, seq[pos_seq : (pos_seq + len[i] - 1)])
      map_seq <- c(map_seq, seq[pos_seq : (pos_seq + len[i] - 1)])
      pos_seq <- pos_seq + len[i]
    }
    if(type[[i]] == "D"){
      ref_seq <- c(ref_seq, rep("N", len[i]))
      map_seq <- c(map_seq, rep("-", len[i]))
    }
    if(type[[i]] == "I"){
      ref_seq <- c(ref_seq, rep("-", len[i]))
      map_seq <- c(map_seq, str_to_lower(seq[pos_seq : (pos_seq + len[i] - 1)]))
      pos_seq <- pos_seq + len[i]
    }
    if(type[[i]] == "S"){
      ref_seq <- c(ref_seq, str_to_lower(seq[pos_seq : (pos_seq + len[i] - 1)]))
      map_seq <- c(map_seq, str_to_lower(seq[pos_seq : (pos_seq + len[i] - 1)]))
      pos_seq <- pos_seq + len[i]
    }
  }
  return(list(ref = ref_seq, seq = map_seq))
}
getSub <- function(fixseq, cutpos, leftSize, rightSize){
  index <- 1
  step <- 1
  ref <- fixseq$ref
  seq <- fixseq$seq
  while(index < cutpos){
    if(step > length(ref)){
      return(list(left = "", right = ""))
    }
    if(str_detect(ref[step], "[A-Z]")){
      index <- index + 1
    }
    step <- step + 1
  }
  index <- step
  validAcc <- 0
  leftRes <- c()
  while(validAcc < leftSize){
    if(index < 1){
      break
    }
    leftRes <- c(seq[index], leftRes)
    if(str_detect(ref[index], "[A-Z]"))
      validAcc <- validAcc + 1
    index <- index - 1
  }
  index <- step + 1
  validAcc <- 0
  rightRes <- c()
  while(validAcc < rightSize){
    if(index > length(ref)){
      break
    }
    rightRes <- c(rightRes, seq[index])
    if(str_detect(ref[index], "[A-Z]"))
      validAcc <- validAcc + 1
    index <- index + 1
  }
  leftRes <- paste(leftRes, collapse = "")
  rightRes <- paste(rightRes, collapse = "")
  if(!str_detect(leftRes, "[A-Z]")){
    leftRes <- ""
  }
  if(!str_detect(rightRes, "[A-Z]")){
    rightRes <- ""
  }
  return(list(left = leftRes, right = rightRes))
}
calMinDel <- function(cigar, start_pos){
  len <- unlist(strsplit(cigar, "[MXSID]"))
  len <- as.numeric(len[len != ""])
  type <- unlist(strsplit(cigar, "[0-9]+"))
  type <- type[type != ""]
  possible_len <- 0
  mapping_last <- start_pos - 1
  for(i in 1 : length(len)){
    if(type[i] %in% c("I", "S")){
      possible_len <- possible_len + len[i]
    } else if(type[i] == "M"){
      possible_len <- possible_len + len[i]
      mapping_last <- mapping_last + len[i]
    } else if(type[i] == "D"){
      mapping_last <- mapping_last + len[i]
    }
  }
  return(mapping_last - possible_len)
}
getRegionDel <- function(fixseq, start, end, checkMap = F){
  index <- 1
  step <- 1
  ref <- fixseq$ref
  seq <- fixseq$seq
  while(index < start){
    if(step > length(ref)){
      return(-999)
    }
    if(str_detect(ref[step], "[A-Z]")){
      index <- index + 1
    }
    step <- step + 1
  }
  delLen <- 0
  checked <- T
  #check mapped before del
  if(checkMap){
    if(!(str_detect(seq[step - 1], "[A-Z]") & str_detect(ref[step - 1], "[A-Z]"))){
      checked <- F
    }
  }
  
  while(index <= end){
    if(step > length(ref)){
      return(-999)
    }
    if(str_detect(ref[step], "[A-Z]")){
      index <- index + 1
    }
    if(ref[step] == "N"){
      delLen <- delLen + 1
    }
    step <- step + 1
  }
  if(!checked){
    return(-delLen) 
  }else{
    return(delLen)
  }
  
}
getRegionDel2 <- function(fixseq, start, NLen = 6){
  index <- 1
  step <- 1
  ref <- fixseq$ref
  seq <- fixseq$seq
  while(index < start){
    if(step > length(ref)){
      return(-999)
    }
    if(str_detect(ref[step], "[A-Z]")){
      index <- index + 1
    }
    step <- step + 1
  }
  delLen <- 0
  acc <- 0
  #check mapped before del
  while(step > 1){
    step <- step - 1
    acc <- acc + 1
    if(ref[step] == "N"){
      delLen <- delLen + 1
    } else {
      if(str_detect(ref[step], "[A-Z]") && acc > NLen){
        break
      }
    }
    
    
  }
  return(delLen)
  
}
reverse_complement_column <- function(dna_sequences) {
  # 保存原始行数
  original_length <- length(dna_sequences)
  
  # 记录非空行的位置
  non_empty_indices <- which(dna_sequences != "")
  
  # 移除空行
  dna_sequences <- dna_sequences[non_empty_indices]
  
  # 如果移除空行后仍有NA值，返回NA
  if (any(is.na(dna_sequences))) return(rep(NA, original_length))
  
  # 处理DNA序列
  result <- tryCatch({
    dna_string_set <- Biostrings::DNAStringSet(dna_sequences)
    reverse_complement_set <- Biostrings::reverseComplement(dna_string_set)
    as.character(reverse_complement_set)
  }, error = function(e) {
    warning("Error processing DNA sequences: ", conditionMessage(e))
    return(rep(NA, length(dna_sequences)))
  })
  
  # 创建一个与原始长度一致的向量，填充NA
  full_result <- rep(NA, original_length)
  
  # 将处理后的结果填充回原始位置
  full_result[non_empty_indices] <- result
  
  return(full_result)
}

library(stringr)

# 获取所有子文件夹中的result.txt文件路径
file_paths <- list.files(path = data_path, pattern = "result.txt", recursive = TRUE, full.names = TRUE)

# 过滤出文件名为“result.txt”的文件路径
file_paths <- file_paths[basename(file_paths) == "result.txt"]

process_file <- function(file_path){cat("Processing file:", file_path, "\n")  # 添加调试信息
  merged_reads <- read.table(file_path, header = F, fill = T, sep = "\t")
  colnames(merged_reads) <- c("Seq", "Counts", "Cigar","Start Pos", "Length", "Indel Times", "Indel Len", "Cigar2","IndelInN")
  merged_reads <- merged_reads[order(merged_reads$Counts, decreasing = T),]
  
  merged_stat <- data.frame(type = c("Total Reads", ">3 Reads"), 
                            counts = c(sum(merged_reads$Counts), 
                                       sum(merged_reads$Counts[merged_reads$Counts > 3])))
  sample <- str_remove(file_path, "_result.txt")
  write.table(merged_stat, file = paste0(sample, "_stat.txt"), row.names = FALSE, col.names = TRUE, sep = "\t")
  filtered_reads <- merged_reads[merged_reads$Counts > 3,]
  filtered_reads <- lapply(split(filtered_reads, filtered_reads$Cigar2), function(x){
    x <- x[order(x$Counts, decreasing = T),]
    x$Counts <- sum(x$Counts)
    x[1,]
  })
  filtered_reads <- data.frame(do.call(rbind, filtered_reads))
  filtered_reads <- filtered_reads[order(filtered_reads$Counts, decreasing = T),]
  fixseqs <- lapply(1 : nrow(filtered_reads), function(index){
    fixSeq2(filtered_reads$Seq[index],filtered_reads$Cigar2[index], filtered_reads$Start.Pos[index], 
            reference = str_to_upper("CACGAGTGTCCATAGGTCACtgctcctcctctgggtttccatgtttacccgtttccagtcatcaggtgAGGCCGTATAGAGCGTCGTGTAGGGTTAGAGTGTCGATCTCGGAAGAGTCGCCTGATTC"), 
            apos = 69)})
  seqs <- lapply(1 : nrow(filtered_reads), function(index){
    getSub(fixseqs[[index]],  60, 10, 18)
  })
  filtered_reads$leftSeq <- unlist(lapply(seqs, function(x){
    x$left
  }))
  filtered_reads$rightSeq <- unlist(lapply(seqs, function(x){
    x$right
  }))
  
  filtered_reads$delLen <- unlist(lapply(1 : nrow(filtered_reads), function(index){
    getRegionDel2(fixseqs[[index]], 69, 1)
  }))
  
  write.table(filtered_reads, file = paste0(sample, "_after_repairmatch_filtered.txt"), row.names = FALSE, col.names = TRUE, sep = "\t")
  
  pattern <- "^(Seq|leftSeq|rightSeq)"
  dna_columns <- grep(pattern, names(filtered_reads), value = TRUE)
  
  filtered_reads_reverse_complement <- filtered_reads
  filtered_reads_reverse_complement[dna_columns] <- lapply(filtered_reads[dna_columns], reverse_complement_column)
  
  names(filtered_reads_reverse_complement)[names(filtered_reads_reverse_complement) == "leftSeq"] <- "rightseq"
  names(filtered_reads_reverse_complement)[names(filtered_reads_reverse_complement) == "rightSeq"] <- "leftseq"
  
  filtered_reads_reverse_complement <- filtered_reads_reverse_complement[, c("Seq", "Counts", "Cigar", "Start.Pos", "Length", "Indel.Times", "Indel.Len", "Cigar2", "IndelInN", "leftseq", "rightseq", "delLen")]
  
  write.table(filtered_reads_reverse_complement, file = paste0(sample, "_reverse_after_repairmatch_filtered.txt"), row.names = FALSE, col.names = TRUE, sep = "\t")
}

for (file_path in file_paths) {
  process_file(file_path)}
