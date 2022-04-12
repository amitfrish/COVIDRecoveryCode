enrichmentTest = function(corsForGSEA,geneSetsListSelectedLarge,direction){
  sapply(geneSetsListSelectedLarge, function(x){
    mutualGenes = intersect(x,names(corsForGSEA))
    inSetValues = corsForGSEA[mutualGenes]
    outSetValues = corsForGSEA[!(names(corsForGSEA) %in% mutualGenes)]
    outSetValues = outSetValues[!is.na(outSetValues)]
    if(length(mutualGenes)<10){
      NA
    }else{
      ks.test(outSetValues,inSetValues,alternative = direction)$p.value  
    }
  })
}

calculateMeanCor = function(corsForGSEA, geneSetsList){
  selectedGeneSetsLargeIndexes = c(grep("^HALLMARK",names(geneSetsList)),grep("^REACTOME",names(geneSetsList)),grep("^KEGG",names(geneSetsList)),grep("^GO",names(geneSetsList)))
  geneSetsListSelectedLarge = lapply(selectedGeneSetsLargeIndexes,function(currSet){toupper(geneSetsList[[currSet]])})
  names(geneSetsListSelectedLarge) = names(geneSetsList)[selectedGeneSetsLargeIndexes]
  sapply(geneSetsListSelectedLarge, function(x){
    mutualGenes = intersect(x,names(corsForGSEA))
    inSetValues = corsForGSEA[mutualGenes]
    if(length(mutualGenes)<10){
      NA
    }else{
      median(inSetValues)
    }
  })
}

runEnrichment = function(corsForGSEA, geneSetsList, singCutOff = 0.01, justHallmark = T){
  selectedGeneSetsLargeIndexes = grep("^HALLMARK",names(geneSetsList))
  if(!justHallmark){
    selectedGeneSetsLargeIndexes = c(grep("^HALLMARK",names(geneSetsList)),grep("^REACTOME",names(geneSetsList)),grep("^KEGG",names(geneSetsList)),grep("^GO",names(geneSetsList)))
  }
  geneSetsListSelectedLarge = lapply(selectedGeneSetsLargeIndexes,function(currSet){toupper(geneSetsList[[currSet]])})
  names(geneSetsListSelectedLarge) = names(geneSetsList)[selectedGeneSetsLargeIndexes]
  posEnriched = enrichmentTest(corsForGSEA,geneSetsListSelectedLarge,"greater")
  negEnriched = enrichmentTest(corsForGSEA,geneSetsListSelectedLarge,"less")
  totalEnriched = cbind(posEnriched,negEnriched)
  totalEnrichedAdjusted = p.adjust(totalEnriched, method = "fdr")
  totalEnrichedLog = as.matrix(cbind(-log(totalEnrichedAdjusted[1:length(posEnriched)],base = 10), log(totalEnrichedAdjusted[(length(posEnriched)+1):length(totalEnrichedAdjusted)],base=10)))
  bestSelection = unlist(lapply(apply(abs(totalEnrichedLog),1,which.max),function(currScore){
    if(length(currScore)==0){
      0
    }else{
      currScore
    }
  }))
  enrichemtntDF = as.matrix(totalEnrichedLog[,1])
  enrichemtntDF[which(bestSelection==2),] = totalEnrichedLog[which(bestSelection==2),2]
  colnames(enrichemtntDF) = "Score"
  row.names(enrichemtntDF) = names(geneSetsList)[selectedGeneSetsLargeIndexes]
  enrichemtntDF
}
