
#These functions are from file: fold.crispr_version5_script.R
#Used to analyze single dataset (real data used by mageck paper).
get.median = function(x, y, func)
{
  return(tapply(y, x, func))    ############### mean vs median
}

#get permutated p.values. x is an element, y is a vector contains all permutated values.
get.pos.permu.p = function(x, y)
{
  pos.p = sum(y > x, na.rm = T)/length(y)
  return(pos.p)
}

get.neg.permu.p = function(x, y)
{
  neg.p = sum(y < x, na.rm = T)/length(y)
  return(neg.p)
}

#zs.gene: calculated stat for each gene according to the function
permu.pvalue = function(dat, sim.no = sim.no, zs.gene, func, seed = 7292016)
{

  zs = log(dat[,4]/dat[,3])
  set.seed(seed)
  index.mat = matrix(rep(dat$Gene, sim.no), ncol = sim.no)
  permu.index.mat = apply(index.mat, 2, sample)    #each column is a permutation of the original gene index.
  result.mat = apply(permu.index.mat, 2, get.median, y = zs, func)

  pos.p = sapply(zs.gene, get.pos.permu.p, as.vector(result.mat))
  neg.p = sapply(zs.gene, get.neg.permu.p, as.vector(result.mat))

  ##p.val = apply(xx, 1, sum, na.rm = T)/sim.no
  #pos.p.val = sapply(zs.gene, inter.pos.qvalue, result.mat)/sim.no/rank(-zs.gene)
  #neg.p.val = sapply(zs.gene, inter.neg.qvalue, result.mat)/sim.no/rank(-zs.gene)
  return(data.frame(pos.p, neg.p))
}



#dataset.id is the id number for simulated dataset name; off.ratio is the off target proportion in the simulated dataset

fold.crispr = function(dat,sim.no = 10, func = "median", alpha.threshold = .2)
{
  dat = dat[order(dat$Gene), ] # so that data set is in increasing order of gene
  dat[is.na(dat)] = 0
  dat[, 3:4] = dat[, 3:4] + .25
  datt = dat
  dat[, 3] = datt[, 3] * mean(c(sum(datt[, 3]), sum(datt[, 4])))/sum(datt[, 3])
  dat[, 4] = datt[, 4] * mean(c(sum(datt[, 3]), sum(datt[, 4])))/sum(datt[, 4])

  zs.gene = tapply(log(dat[,4]/dat[,3]), dat$Gene, func) #################### mean vs median

  initial.p.value = permu.pvalue(dat, sim.no = sim.no, zs.gene = zs.gene, func = func)
  initial.adj.pos.pvalue = initial.p.value$pos.p
  initial.adj.neg.pvalue = initial.p.value$neg.p

  #direction = sign(tapply(log(dat[,4]/dat[,3]), dat$Gene, median))
  initial.result = data.frame(Gene = sort(unique(dat$Gene)), initial.adj.pos.pvalue, initial.adj.neg.pvalue) #, direction)
  #hit = final.result[final.result$adj.p.value < threshold,]

  #summarize genes selected by the program
  initial.pos.gene = initial.result$Gene[initial.adj.pos.pvalue < alpha.threshold]         ############
  initial.neg.gene = initial.result$Gene[initial.adj.neg.pvalue < alpha.threshold]         ############
  update.dat = dat[!is.element(dat$Gene, c(initial.pos.gene, initial.neg.gene)),]

  p.value = permu.pvalue(update.dat, sim.no = sim.no, zs.gene = zs.gene, func = func)
#  adj.pos.pvalue = p.adjust(p.value$pos.p, method = method)
#  adj.neg.pvalue = p.adjust(p.value$neg.p, method = method)

  #direction = sign(tapply(log(dat[,4]/dat[,3]), dat$Gene, median))
  final.result = data.frame(Gene = sort(unique(dat$Gene)), pos.pvalue = p.value$pos.p, neg.pvalue = p.value$neg.p)
                            #pos.fdr = adj.pos.pvalue, neg.fdr = adj.neg.pvalue)

#  pos.gene = final.result$Gene[adj.pos.pvalue < fdr]  #genes that are selected as positive genes.
#  neg.gene = final.result$Gene[adj.neg.pvalue < fdr]  #genes that are selected as negative genes.
#  pos.no = length(pos.gene)
#  neg.no = length(neg.gene)


  #output = data.frame(id = 1:length(unique(dat$Gene)), pos.p.value = p.value$pos.p, pos.fdr = adj.pos.pvalue, neg.p.value = p.value$pos.p, neg.fdr = adj.neg.pvalue)

  #write.table(output, file = paste("output_offratio", off.ratio, "_simdata", dataset.id, ".txt", sep = ""), quote = F, row.names = F )

  #return.value = list(pos.gene, pos.no, neg.gene, neg.no, final.result)

  return(final.result)
}



#' @title Permutation Based Non-Parametric Analysis of CRISPR Screen Data
#'
#' @description
#' This function reads the raw read count data and conducts statistical
#' analysis for permutation based non-parametric analysis of CRISPR screen data.
#'
#' @param dat List type with each element being the raw read count data for one replicate.
#' Each element should be a dataframe with four columns. The first column is named
#' 'sgRNA' which is the sgRNA index; the second column is named 'Gene' which is the
#' gene index; the third column should be the initial read count or control read
#' count and the fourth column should be the final read count or treatment read count.
#' Missing values in the read count are replaced with 0.
#' @param sim.no Number of permutations used to get the un-adjusted p-value.Set to 10 by default.
#' @param alpha.threshold Threshold to remove genes with significant p-values. Set to 0.2 by default.
#' @param fdr The FDR threshold to determine the selected genes. Set to 0.05 by default.
#' @details PBNPA implements permutation based non-parametric analysis of CRISPR screen data. Details
#' about this algorithm are published in the following paper published on BMC genomics, Jia et al. (2017) <doi:10.1186/s12864-017-3938-5>: A permutation-based non-parametric analysis of CRISPR screen data.
#' Please cite this paper if you use this algorithm for your paper.
#'
#' @return A list of 5 elements will be returned. The first element is pos.gene, which is the index of
#' genes identified as hits for positive screen by controlling FDR at the selected level; the second
#' element is pos.number, which is the number of genes identified as hits for positive screen; The
#' third element is neg.gene, which is the index of genes identified as hits for negative screen by
#' controlling FDR at the selected level; the fourth element is neg.number, which is the number of genes
#' identified as hits for negative screen; the fifth element is a dataframe which contains unadjusted
#' p-values and FDR adjusted p-values for all the genes (for both negative selection and positive selection).
#' @export
#' @examples
#' dat11 = system.file('extdata','simdata_20per_off50.csv', package='PBNPA')
#' dat22 = system.file('extdata','simdata_20per_off49.csv', package='PBNPA')
#' dat33 = system.file('extdata','simdata_20per_off48.csv', package='PBNPA')
#' dat1 = read.csv(dat11, header = TRUE)
#' dat2 = read.csv(dat22, header = TRUE)
#' dat3 = read.csv(dat33, header = TRUE)
#' datlist = list(dat1, dat2, dat3)
#' result = PBNPA(datlist)
#' @import metaRNASeq
PBNPA = function(dat, sim.no = 10, alpha.threshold = .2, fdr = .05)
{
  nrep = length(dat)
  combine.pos = list()
  combine.neg = list()
  for (i in 1:nrep)
  {
    result = fold.crispr(dat[[i]], sim.no = sim.no, alpha.threshold = alpha.threshold)
    combine.pos[[i]] = result$pos.pvalue
    combine.neg[[i]] = result$neg.pvalue
  }

  #library(metaRNASeq)

  combined.pos = metaRNASeq::fishercomb(combine.pos, BHth = fdr)
  combined.neg = metaRNASeq::fishercomb(combine.neg, BHth = fdr)

  final.result = data.frame(Gene = sort(unique(dat[[1]]$Gene)), pos.pvalue = combined.pos$rawpval, pos.fdr = combined.pos$adjpval,
                            neg.pvalue = combined.neg$rawpval, neg.fdr = combined.neg$adjpval)
  pos.gene = final.result$Gene[combined.pos$DEindices]  #genes that are selected as positive genes.
  neg.gene = final.result$Gene[combined.neg$DEindices]  #genes that are selected as negative genes.
  pos.no = length(pos.gene)
  neg.no = length(neg.gene)

  return.value = list(pos.gene = pos.gene, pos.no = pos.no, neg.gene = neg.gene, neg.no = neg.no, final.result = final.result)
  return(return.value)
}


