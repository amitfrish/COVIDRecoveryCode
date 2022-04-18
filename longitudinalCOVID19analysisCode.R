#folder = "<<Code location>>"
folder = "C:/Users/Shennor/Dropbox (Irit Gat Viks)/phD_documents/Postdoc/Joachim/COVIDLongitudinal/uploadedCode"
folder = file.path(folder,"FileFolder")
library(ggplot2)
library(sva)
library(TimeAx)
library(NMF)
library(scBio)
dataPlotter = function(currMetaData, currTraj, currMeta, additionalCol = NULL, box = T, yValues = NULL, xValues = NULL, regLine = F, DoF = 4){
  x = currMetaData[,currMeta]
  indsToUse = which(!is.na(x) & !is.na(currTraj))
  metaDataLevels = x[indsToUse]
  currTrajLevels = currTraj[indsToUse]
  currCorValue = cor(metaDataLevels, currTrajLevels)
  if(length(unique(metaDataLevels))<6 & box){
    if(length(unique(metaDataLevels)) == 2 & min(table(metaDataLevels))>3){
      currCorValue = t.test(currTrajLevels[metaDataLevels == unique(metaDataLevels)[1]],
                            currTrajLevels[metaDataLevels == unique(metaDataLevels)[2]])$p.value
    }
    ggplot(data = NULL, aes(x = factor(metaDataLevels), y = currTrajLevels))+geom_boxplot()+ggtitle(paste(currMeta, currCorValue,sep = " _ "))+
      theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
  }else{
    if(is.null(xValues)){
      xValues = c(min(currTrajLevels),max(currTrajLevels))
    }
    if(is.null(yValues)){
      yValues = c(min(metaDataLevels),max(metaDataLevels))
    }
    regLineFunc = NULL
    if(regLine){
      regLineFunc = geom_smooth(method = "lm", formula = y ~ poly(x, DoF), se = F)
    }
    if(!is.null(additionalCol)){
      additionalCol = factor(additionalCol[indsToUse])
      ggplot(data = NULL, aes(x = currTrajLevels, y = metaDataLevels))+geom_point(size = 3, aes(color = additionalCol))+ggtitle(paste(currMeta, currCorValue,sep = " _ "))+xlim(xValues)+ylim(yValues)+regLineFunc+
        theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
    }else{
      ggplot(data = NULL, aes(x = currTrajLevels, y = metaDataLevels))+geom_point(size = 3)+ggtitle(paste(currMeta, currCorValue,sep = " _ "))+xlim(xValues)+ylim(yValues)+regLineFunc+
        theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
    }
  }
}

#### Load and norm cohort 1 data ####
GEData = read.table(file.path(folder, "countMatrix.txt"),sep="\t", header = T, row.names = 1, check.names = F)
sampleDict = read.table(file.path(folder, "sampleDict.txt"),sep="\t", header = T, row.names = 1, check.names = F)
trajectoryBase = as.numeric(sapply(colnames(GEData),function(x){unlist(strsplit(x,"_"))[2]}))-1
sampleNames = sapply(colnames(GEData),function(x){unlist(strsplit(x,"_"))[1]})
batchInformation = sampleDict$Batch_num
dateNumbers <- as.numeric(as.Date(sampleDict$Sample_Name_02, format = "%m/%d/%Y", origin = "2020-03-20"))
metaDataNum = readRDS(file.path(folder,"metaDataSmall.rds"))

y <- edgeR::DGEList(counts=GEData) # Put counts into DGEList.
y <- edgeR::calcNormFactors(y) # TMM-normalization.
logCPM <- edgeR::cpm(y, log=F, prior.count=0.1) # Convert to CPM and log2 transformation.
GEDataCombat = ComBat(logCPM, batchInformation)

#### Create train and test ####
chosenSampleNamesTrain = names(which(table(sampleNames)>4))
deathInChosen = sapply(chosenSampleNamesTrain, function(currName){max(metaDataNum[sampleNames==currName,"Pat_deceased_in_hosp"])})
chosenSampleNamesTrain = chosenSampleNamesTrain[deathInChosen == 1]
indexesForTrain = which((sampleNames %in% chosenSampleNamesTrain))
sampleNamesTrain = sampleNames[indexesForTrain]
GEDataCombatTrain = GEDataCombat[,indexesForTrain]
indexesForTest = which(!(sampleNames %in% chosenSampleNamesTrain))
sampleNamesTest = sampleNames[indexesForTest]
GEDataCombatTest = GEDataCombat[,indexesForTest]

indexOrderAll = order(c(indexesForTrain,indexesForTest))
GEDataCombatAll = cbind(GEDataCombatTrain,GEDataCombatTest)[,indexOrderAll]
sampleNamesAll = c(sampleNamesTrain,sampleNamesTest)[indexOrderAll]
dateNumbersAll = dateNumbers[c(indexesForTrain,indexesForTest)][indexOrderAll]

dateListAll = lapply(unique(sampleNames),function(currName){
  dateNumbersAll[currName == sampleNamesAll]
})

#### Create original time scales ####
daysFromInclusion = metaDataNum[,"daysFromInclusion"]
daysFromSymph = metaDataNum[,"daysFromSymph"]

#### Alignment ####
proteinCodingGenes = as.character(as.matrix(read.table(file.path(folder,"proteinCodingGenes.txt"),sep = "\t", row.names = NULL, col.names = F)))
proteinCodingGenes = intersect(proteinCodingGenes, row.names(GEDataCombatAll))
multiAlignment = TimeAx::modelCreation(GEDataCombatTrain[proteinCodingGenes,], sampleNamesTrain)
profilePredictionsFull = predictByConsensus(multiAlignment,GEDataCombatAll)
trajToUseProfileToUse = profilePredictionsFull$predictions
originalTraj = daysFromSymph
model1Robustness = TimeAx::robustness(multiAlignment,GEDataCombatAll,sampleNamesAll,profilePredictionsFull)

#### Figure 2 ####
#Figure 2B
prePostData = as.data.frame(cbind(originalTraj,trajToUseProfileToUse))
colnames(prePostData) = c("Pre","Post")
prePostData$Sample = as.factor(sampleNamesAll)
prePostData = prePostData[prePostData$Pre>0,]
samplesToShow = c("COR031","COR027","COR062","COR069")
ggplot(prePostData[prePostData$Sample %in% samplesToShow,], aes(x = Pre, y = Post, colour = Sample))+geom_line()+
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

#Figure 2C
metasForP = c("CRP","Creatinine","IL6","GCSF")
metaP = sapply(metasForP, function(currMeta){
  currX = metaDataNum[,currMeta]
  summary(lm(currX[!is.na(currX)]~trajToUseProfileToUse[!is.na(currX)]))[[4]][2,4]
})
metaPOrig = sapply(metasForP, function(currMeta){
  currX = metaDataNum[,currMeta]
  summary(lm(currX[!is.na(currX) & daysFromSymph>0]~daysFromSymph[!is.na(currX)  & daysFromSymph>0]))[[4]][2,4]
})
metaDataCor = apply(metaDataNum[,metasForP],2,function(x){
  indsToUse = which(!is.na(x))
  if(length(indsToUse)>40){
    cor(trajToUseProfileToUse[indsToUse], x[indsToUse])
  }else{
    NA
  }
})
metaDataCorOrig = apply(metaDataNum[,metasForP],2,function(x){
  indsToUse = which(!is.na(x) & daysFromSymph>0)
  if(length(indsToUse)>40){
    cor(daysFromSymph[indsToUse], x[indsToUse])
  }else{
    NA
  }
})
dataPlotter(metaDataNum, trajToUseProfileToUse,"CRP",regLine = T)
dataPlotter(metaDataNum[originalTraj>0,], originalTraj[originalTraj>0],"CRP",regLine = T)

#Figure 2D
MVSGenes = read.table(file.path(folder,"MVSGenes.txt"),sep="\t", header = T, row.names = 1, check.names = F)
MVSGenes = MVSGenes[intersect(row.names(MVSGenes),row.names(GEDataCombatAll)),,drop=F]
MVSScores = apply(GEDataCombatAll[row.names(MVSGenes),],2,function(x){
  overExp = x[MVSGenes$EffectSize>0]
  underExp = x[MVSGenes$EffectSize<0]
  overExp = overExp[overExp>0]
  underExp = underExp[underExp>0]
  exp(mean(log(overExp)))-exp(mean(log(underExp)))
})

dataForMVSNMF = t(GEDataCombatAll[row.names(MVSGenes[MVSGenes$EffectSize>0,,drop=F]),])
dataForMVSNMF = dataForMVSNMF-min(dataForMVSNMF)
MVSModel = NMF::nmf(dataForMVSNMF, 1)@fit
MVSScores = MVSModel@W[,1]

MVSP = summary(lm(MVSScores~trajToUseProfileToUse))[[4]][2,4]
MVSPOrig = summary(lm(MVSScores[daysFromSymph>0]~daysFromSymph[daysFromSymph>0]))[[4]][2,4]

ggplot(data = NULL, aes(x = trajToUseProfileToUse, y = MVSScores))+geom_point(size = 3)+geom_smooth(method = "lm", formula = y ~ poly(x, 4), se = F)+ggtitle(cor(MVSScores,trajToUseProfileToUse))+
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
ggplot(data = NULL, aes(x = originalTraj[originalTraj>0], y = MVSScores[originalTraj>0]))+geom_point(size = 3)+geom_smooth(method = "lm", formula = y ~ poly(x, 4), se = F)+ggtitle(cor(MVSScores[originalTraj>0],originalTraj[originalTraj>0]))+
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

#Figure 2F
corGeneTraj = cor(t(GEDataCombatAll),trajToUseProfileToUse)
corGeneTrajPre = cor(t(GEDataCombatAll)[originalTraj>0,],originalTraj[originalTraj>0])
nonNaValues = which(!is.na(corGeneTraj))
corGeneTraj = corGeneTraj[names(corGeneTraj[nonNaValues,]),]
corGeneTrajPre = corGeneTrajPre[names(corGeneTrajPre[nonNaValues,]),]

cutOff = 0.3
geneWithCorOf0.5 = which(abs(corGeneTraj)>cutOff | abs(corGeneTrajPre)>cutOff)
diffInCor = (abs(corGeneTraj)-abs(corGeneTrajPre))[geneWithCorOf0.5]
diffInCorPlot = ggplot(data = NULL, aes(x = diffInCor)) + geom_density()+geom_vline(xintercept = 0, linetype="dashed")+
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

#Figure 2G
regPre = summary(lm(t(GEDataCombatAll)[daysFromSymph>0,]~daysFromSymph[daysFromSymph>0]))
genePValuesPre = sapply(regPre, function(x){
  x[[4]][2,4]
})

regAligned = summary(lm(t(GEDataCombatAll)~trajToUseProfileToUse))
genePValuesAligned = sapply(regAligned, function(x){
  x[[4]][2,4]
})

geneRegPre = p.adjust(genePValuesPre, method = "fdr")
geneRegAligned = p.adjust(genePValuesAligned, method = "fdr")

geneRegPreLog = -log(geneRegPre,base=10)
geneRegAlignedLog = -log(geneRegAligned,base=10)
newCutOff = -log(0.01, base=10)
#newCutOff = 5
geneColourProVsAligned = rep(0, length(geneRegPreLog))
geneColourProVsAligned[geneRegPreLog>newCutOff & geneRegAlignedLog<newCutOff] = 1
geneColourProVsAligned[geneRegPreLog<newCutOff & geneRegAlignedLog>newCutOff] = 2
geneColourProVsAligned[geneRegPreLog>newCutOff & geneRegAlignedLog>newCutOff] = 3
ggplot(data = NULL, aes(x = geneRegPreLog, y = geneRegAlignedLog, colour = factor(geneColourProVsAligned)))+geom_point()+
  geom_hline(yintercept = newCutOff, linetype = "dashed")+geom_vline(xintercept = newCutOff, linetype = "dashed")+
  scale_color_manual(values = c("0" = "lightgray","1" = "darkgray","2" = "brown", "3" = "green"))+
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

#Figure 2H
geneMatrixForComp = cbind(corGeneTraj,corGeneTrajPre)
geneMatrixForComp = geneMatrixForComp[abs(geneMatrixForComp[,1])>0.4 & abs(geneMatrixForComp[,2])<0.2,]
colnames(geneMatrixForComp) = c("Pseudo","Time from admission")

geneDiff = sort(geneMatrixForComp[,1]-geneMatrixForComp[,2])
geneMatrixForCompToShow = geneMatrixForComp[names(c(geneDiff[1:10],rev(geneDiff)[1:10])),]

geneToPlot = "IL17RA"
dataPlotter(t(GEDataCombatAll), trajToUseProfileToUse,geneToPlot,regLine = T)
dataPlotter(t(GEDataCombatAll)[daysFromSymph>0,], daysFromSymph[daysFromSymph>0],geneToPlot,regLine = T)

#### Figure 3 ####
caculateModelFitScore = function(sampleList, pseudo, pre){
  fitScores = sapply(unique(sampleList),function(currSample){
    sd(pseudo[sampleList == currSample])
  })
  names(fitScores) = unique(sampleList)
  fitScores
}

multiAlignmentRev = modelCreation(GEDataCombatTrain[proteinCodingGenes,dim(GEDataCombatTrain)[2]:1], rev(sampleNamesTrain))
model2Robustness = TimeAx::robustness(multiAlignmentRev,GEDataCombatAll,sampleNamesAll)
samplePseuSD = caculateModelFitScore(sampleNamesAll,model1Robustness$robustnessPseudo,daysFromSymph)
samplePseuSDRev = caculateModelFitScore(sampleNamesAll,model2Robustness$robustnessPseudo,daysFromSymph)

#Figure 3B
pseudoCutOff = 0.1
recoverGroup = names(which(samplePseuSD>pseudoCutOff & samplePseuSDRev<pseudoCutOff))
worsenGroup = names(which(samplePseuSD<pseudoCutOff & samplePseuSDRev>pseudoCutOff))
sameGroup = names(which(samplePseuSD<pseudoCutOff & samplePseuSDRev<pseudoCutOff))
strangeGroup = names(which(samplePseuSD>pseudoCutOff & samplePseuSDRev>pseudoCutOff))

stateVector = c(rep(-1, length(worsenGroup)),rep(0, length(sameGroup)),rep(1, length(recoverGroup)))
sampleOrderForGroups = c(worsenGroup,sameGroup,recoverGroup)
names(stateVector) = sampleOrderForGroups
ggplot(data = NULL, aes(x = samplePseuSD, y = samplePseuSDRev,colour = factor(stateVector[names(samplePseuSD)])))+geom_point(size = 3)+
  geom_vline(xintercept = pseudoCutOff)+geom_hline(yintercept = pseudoCutOff)+
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

#Figure 3C-E
AnnotationsToUseForSamplesFirst = do.call(rbind,lapply(unique(sampleNamesAll), function(currSample){
  metaDataNum[sampleNamesAll==currSample,][1,]
}))
row.names(AnnotationsToUseForSamplesFirst) = unique(sampleNamesAll)

getMetaGroupPlot = function(metaName){
  tableOfMetaPerGroup = table(stateVector,AnnotationsToUseForSamplesFirst[sampleOrderForGroups,metaName])
  currP = fisher.test(tableOfMetaPerGroup)$p.value
  chanceForMetaPerGroup = tableOfMetaPerGroup[,2]/rowSums(tableOfMetaPerGroup)
  ggplot(data = NULL, aes(x = factor(unique(stateVector)),y = chanceForMetaPerGroup))+geom_col()+ggtitle(paste(metaName,currP,sep=" _ "))+
    theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
}

deathGroupPlot = getMetaGroupPlot("Pat_deceased_in_hosp") # Provided only clinical outcome and not age and gender due to privacy (GDPR) issues.

#Figure 3F
diffIn2First = sapply(unique(sampleNamesAll),function(currSample){
  currTraj = trajToUseProfileToUse[currSample==sampleNamesAll]
  currTraj[2]-currTraj[1]
})
diffIn2FirstDF = data.frame(diff = diffIn2First, group = rep("Recovery", length(diffIn2First)))
diffIn2FirstDF[worsenGroup,"group"] = "Worsening"
diffIn2FirstDF[sameGroup,"group"] = "Same"
diffIn2FirstDF[strangeGroup,"group"] = "Strange"

diffIn2FirstNoStrangeDF = diffIn2FirstDF[diffIn2FirstDF$group!="Strange",]
diffIn2FirstP = summary(lm(diffIn2FirstNoStrangeDF$diff~as.numeric(as.factor(diffIn2FirstNoStrangeDF$group))))[[4]][2,4]
ggplot(data = diffIn2FirstNoStrangeDF, aes(x = factor(group), y = diff))+geom_boxplot()+ggtitle(diffIn2FirstP)+
  geom_hline(yintercept = 0, linetype = "dashed")+
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

#Figure 3G
dayOfMinPseudo = sapply(sampleOrderForGroups, function(currSample){
  #(daysFromSymph[sampleNamesAll==currSample])[which.min(trajToUseProfileToUse[sampleNamesAll==currSample])]
  (daysFromInclusion[sampleNamesAll==currSample])[which.min(trajToUseProfileToUse[sampleNamesAll==currSample])]
})
dayOfMaxPseudo = sapply(sampleOrderForGroups, function(currSample){
  (daysFromInclusion[sampleNamesAll==currSample])[which.max(trajToUseProfileToUse[sampleNamesAll==currSample])]
})
ggplot(data = NULL, aes(x = dayOfMinPseudo[stateVector!=0], y=dayOfMaxPseudo[stateVector!=0], colour = factor(stateVector[stateVector!=0])))+geom_point(size = 3)+
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))


#### Figure 4 ####
singleCellData = readRDS(file.path(folder,"SCDataForCPM.rds"))
trajGenesForAnalysis = multiAlignment$seed
mutualTrajGenes = intersect(intersect(trajGenesForAnalysis, row.names(singleCellData)),names(corGeneTraj))
cellsToUse = which(colSums(singleCellData[mutualTrajGenes,])!=0)
singleCellDataFinal = singleCellData[,cellsToUse]
singleCellDataFinal = singleCellDataFinal[,colnames(singleCellDataFinal)!="Mixed_cells"]
cellTypes = colnames(singleCellDataFinal)
cellTypeToUse = names(which(table(cellTypes)>10))
singleCellDataFinalRefined = singleCellDataFinal[,colnames(singleCellDataFinal) %in% cellTypeToUse]
cellTypesRefined = colnames(singleCellDataFinalRefined)

cellSpaceInfo = lapply(unique(cellTypesRefined),function(currCellType){
  selectedCellsInCellType = which(cellTypesRefined==currCellType)
  currCells = singleCellDataFinalRefined[,selectedCellsInCellType]
  otherCells = singleCellDataFinalRefined[,cellTypesRefined!=currCellType]
  numOfNonZeros = apply(currCells,1,function(x){length(x[x!=0])})
  selectedCellsToCheck = names(which(numOfNonZeros>0.4*dim(currCells)[2]))
  currFC = rowMeans(currCells[selectedCellsToCheck,])/rowMeans(otherCells[selectedCellsToCheck,])
  currData = singleCellDataFinalRefined[names(currFC)[currFC>2],selectedCellsInCellType]
  currData = currData[rowSums(currData)>0,]
  pcRes = prcomp(t(currData))$x[,1:2]
  selectedCells = which(pnorm(pcRes[,1],mean(pcRes[,1]),sd(pcRes[,1]))>0.1)
  if(length(selectedCells)>50){
    selectedCells = sample(selectedCells,50)
  }
  list(space = pcRes[selectedCells,],indexes = selectedCellsInCellType[selectedCells])
})

selectedCellsOverall = unlist(lapply(cellSpaceInfo,function(x){x$indexes}))
singleCellDataFinalForCPMRefined = singleCellDataFinalRefined[,selectedCellsOverall]
cellSpaceCPM = do.call(rbind,lapply(cellSpaceInfo,function(x){x$space}))
finalNamesCPM = names(which(table(colnames(singleCellDataFinalForCPMRefined))>20))
tempInds = which(colnames(singleCellDataFinalForCPMRefined) %in% finalNamesCPM)
singleCellDataFinalForCPMRefined = singleCellDataFinalForCPMRefined[,tempInds]
cellSpaceCPM = cellSpaceCPM[tempInds,]
cellNamesCPM = colnames(singleCellDataFinalForCPMRefined)

cpmResDeconv = CPM(singleCellDataFinalForCPMRefined, cellNamesCPM, GEDataCombatAll, cellSpaceCPM, modelSize = 100,quantifyTypes = T,typeTransformation = T)

#Figure 4A
lymph1 = trajToUseProfileToUse[which(metaDataNum[,"Lymphocytes"]<0.8)]
lymph2 = trajToUseProfileToUse[which(metaDataNum[,"Lymphocytes"]>0.8 | metaDataNum[,"Lymphocytes"]<5)]
lymph3 = trajToUseProfileToUse[which(metaDataNum[,"Lymphocytes"]>5)]

lymphDF = as.data.frame(cbind(c(lymph1,lymph2,lymph3),c(rep(1,length(lymph1)),rep(2,length(lymph2)),rep(3,length(lymph3)))))
colnames(lymphDF) = c("Lymph","Type")
ggplot(data = lymphDF,aes(x = factor(Type), y=Lymph))+geom_boxplot()+ggtitle(summary(lm(Lymph~Type,data = lymphDF))[[4]][2,4])+
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

#Figure 4B
NeuToLymphRatio = metaDataNum[,"Neutrophils"]/metaDataNum[,"Lymphocytes"]
NeuToLymphRatioLog = log(NeuToLymphRatio,base = 10)
NeuToLymphRatioDF = as.data.frame(cbind(trajToUseProfileToUse,NeuToLymphRatioLog))
NeuToLymphRatioDF = NeuToLymphRatioDF[!is.na(rowSums(NeuToLymphRatioDF)),]
colnames(NeuToLymphRatioDF) = c("Pseudo","NeuToLymphRatio")
summary(lm(Pseudo~NeuToLymphRatio, NeuToLymphRatioDF))[[4]][2,4]
ggplot(data = NeuToLymphRatioDF, aes(x = Pseudo, y= NeuToLymphRatio))+geom_point(size = 3)+ggtitle(cor(NeuToLymphRatioDF)[1,2])+
  geom_smooth(method = "lm", formula = y ~ poly(x, 4), se = F)+
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

#Figure 4C
cpmResFrac = cpmResDeconv$cellTypePredictions
cpmCorWithTraj = cor(cpmResFrac,trajToUseProfileToUse)[,1]
dataPlotter(cpmResFrac, trajToUseProfileToUse,"MME_Neutrophils",regLine = T)
dataPlotter(cpmResFrac, trajToUseProfileToUse,"CD14 Mono",regLine = T)
dataPlotter(cpmResFrac, trajToUseProfileToUse,"MAIT",regLine = T)

#Figure 4E
deadInds = names(which(AnnotationsToUseForSamplesFirst[,"Pat_deceased_in_hosp"]==2))
deadIndsForRetrain = names(which(table(sampleNamesAll)[deadInds]>3))
samplesWithDeath = c(unique(sampleNamesTrain), deadIndsForRetrain)
multiAlignmentDisrupted = modelCreation(GEDataCombatAll[proteinCodingGenes,sampleNamesAll %in% samplesWithDeath], sampleNamesAll[sampleNamesAll %in% samplesWithDeath])
profilePredictionsDisrupted = predictByConsensus(multiAlignmentDisrupted, GEDataCombatAll)
cpmCorWithTrajDisrupted = cor(cpmResFrac,profilePredictionsDisrupted$predictions)[,1]

namesOfCells = names(cpmCorWithTraj)
namesOfCells[which(abs(cpmCorWithTraj - cpmCorWithTrajDisrupted)<0.2)] = ""
minLevel = min(c(cpmCorWithTraj,cpmCorWithTrajDisrupted))
maxLevel = max(c(cpmCorWithTraj,cpmCorWithTrajDisrupted))
ggplot(data = NULL, aes(x = cpmCorWithTraj, y = cpmCorWithTrajDisrupted))+geom_point(size = 3)+
  xlim(minLevel,maxLevel)+ylim(minLevel,maxLevel)+ geom_line(aes(x = cpmCorWithTraj, y = cpmCorWithTraj))+
  geom_text(aes(label = namesOfCells), vjust = -1, size = 3)+ geom_hline(yintercept = 0, linetype = "dashed")+ geom_vline(xintercept = 0, linetype = "dashed")+
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

#Figure 4F-G
getCellFractions = function(ciberForSum){
  NeutSum = rowSums(ciberForSum[,grep("neutro",tolower(colnames(ciberForSum)))])
  CD4TSum = rowSums(ciberForSum[,grep("CD4",colnames(ciberForSum))])
  CD8TSum = rowSums(ciberForSum[,grep("CD8",colnames(ciberForSum))])
  RestTSum = rowSums(ciberForSum[,match(c("gdT","Treg","MAIT"),colnames(ciberForSum))])
  bloodSum = rowSums(ciberForSum[,match(c("Platelet","Eryth"),colnames(ciberForSum))])
  EosinophilsSum = ciberForSum[,match(c("Eosinophils"),colnames(ciberForSum))]
  DCSum = ciberForSum[,match(c("cDC2"),colnames(ciberForSum))]
  NKSum = rowSums(ciberForSum[,grep("NK",colnames(ciberForSum))])# - ciberForSum[,"NK Proliferating"]
  NKBrightSum = ciberForSum[,grep("NK_CD56bright",colnames(ciberForSum))]
  DNSum = ciberForSum[,grep("dnT",colnames(ciberForSum))]
  monoSum = rowSums(ciberForSum[,grep("Mono",colnames(ciberForSum))])
  PlasmaSum = ciberForSum[,grep("Plasmablast",colnames(ciberForSum))]
  TRegSum = ciberForSum[,grep("Treg",colnames(ciberForSum))]
  #TSum = CD4TSum+CD8TSum+RestTSum+DNSum
  TSum = CD4TSum+CD8TSum+RestTSum
  BSum = rowSums(ciberForSum[,match(c("B memory","B naive","B intermediate","Plasmablast"),colnames(ciberForSum))])
  Lymph = rowSums(ciberForSum)-bloodSum-NeutSum-EosinophilsSum-DCSum-monoSum
  Myelo = rowSums(ciberForSum)-bloodSum-Lymph
  lymphNeutRatio = Lymph/NeutSum
  NeutLymphRatio = 1/lymphNeutRatio
  immatureNeut = rowSums(ciberForSum[,grep("-Neut",colnames(ciberForSum))])
  #immatureNeut = ciberForSum[,grep("Pre-Neut",colnames(ciberForSum))]
  totalRealNeut = rowSums(ciberForSum[,grep("_Neut",colnames(ciberForSum))])
  neutRatio = totalRealNeut/immatureNeut
  list(Lymph = Lymph,Myelo = Myelo,totalRealNeut = totalRealNeut,TSum = TSum,
       BSum = BSum,immatureNeut = immatureNeut,neutRatio = neutRatio,NeutLymphRatio = NeutLymphRatio)
}
getFractionStats = function(generalMatrix,currSampleNames, samplesToUse, currState, diff = F){
  cellQuantFirstTime = do.call(rbind,lapply(samplesToUse, function(currSample){
    generalMatrix[currSample==currSampleNames,][1,]
  }))
  row.names(cellQuantFirstTime) = samplesToUse
  
  cellQuantSecondMinusFirst = do.call(rbind,lapply(samplesToUse, function(currSample){
    generalMatrix[currSample==currSampleNames,][2,]-generalMatrix[currSample==currSampleNames,][1,]
  }))
  row.names(cellQuantSecondMinusFirst) = samplesToUse
  
  dataToUseFirst = cellQuantFirstTime
  if(diff){
    dataToUseFirst = cellQuantSecondMinusFirst
  }
  firstCellStats = lapply(colnames(dataToUseFirst),function(currCell){
    currCellDF = as.data.frame(cbind(dataToUseFirst[,currCell],currState))
    colnames(currCellDF) = c("Quantity", "Type")
    currP = summary(lm(Quantity~Type, currCellDF))[[4]][2,4]
    currPlot = ggplot(data = currCellDF, aes(x = factor(Type), y=Quantity))+geom_boxplot()+ggtitle(paste(currCell,currP,sep=" _ "))+
      theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
    pToPrint = -log(currP,base=10)*
      sign(mean(currCellDF$Quantity[currCellDF$Type==1])-mean(currCellDF$Quantity[currCellDF$Type==-1]))
    list(p = pToPrint, plot = list(currPlot))
  })
  
  firstCellP = sapply(firstCellStats, function(x){x$p})
  names(firstCellP) = colnames(dataToUseFirst)
  
  firstCellPlot = sapply(firstCellStats, function(x){x$plot})
  
  firstCellPlot
}
pseudoForRatio = trajToUseProfileToUse
currFractions = getCellFractions(cpmResFrac)
generalMatrix = cbind(currFractions$Lymph,currFractions$Myelo,currFractions$totalRealNeut,
                      currFractions$TSum,currFractions$BSum,currFractions$immatureNeut,
                      currFractions$neutRatio,currFractions$NeutLymphRatio)
colnames(generalMatrix) = c("Lymph","Myelo","MatureNeut","T","B","ImmatureNeut","NeutRatio","NeutLymph")
fractionPlots = getFractionStats(generalMatrix, sampleNamesAll, sampleOrderForGroups, stateVector,diff = T)


#### Figure 5 ####
calculateResiduals = function(geneExpMatrix, cpmModel, cellTypeOnly = T){
  geneMarkersMutualAll = intersect(row.names(singleCellDataFinalRefined), row.names(geneExpMatrix))
  signatureMatrixMutual = signatureMatrixAll[geneMarkersMutualAll,]
  signatureMatrixMutual = t(apply(signatureMatrixMutual,1,function(x){(x-mean(x))/sd(x)}))
  
  predictedCellStates = cpmModel$predicted
  cellTypeTable = table(colnames(predictedCellStates))
  for(currCell in names(cellTypeTable)){
    predictedCellStates[,colnames(predictedCellStates) == currCell]= predictedCellStates[,colnames(predictedCellStates) == currCell] / cellTypeTable[currCell]
  }
  
  cellTypePart = singleCellDataFinalForCPMRefined %*% t(predictedCellStates)
  
  #cellTypePart = signatureMatrixMutual %*% t(cpmModel$cellTypePredictions[,colnames(signatureMatrixMutual)])
  cellTypePart = cellTypePart[geneMarkersMutualAll,]
  
  if(cellTypeOnly){
    return(cellTypePart)
  }
  dataRefined = geneExpMatrix[geneMarkersMutualAll,]
  cellTypePart = cellTypePart-min(cellTypePart)
  dataRefined = dataRefined * (rowMeans(cellTypePart) / rowMeans(dataRefined))
  #dataRefined = t(apply(dataRefined,1,function(x){(x-mean(x))/sd(x)}))
  #cellTypePart = t(apply(cellTypePart,1,function(x){(x-mean(x))/sd(x)}))
  dataResiduals = dataRefined-cellTypePart
  dataResiduals = dataResiduals[which(!is.nan(rowSums(dataResiduals))),]
  dataResiduals-min(dataResiduals)
  #dataResiduals
}
signatureMatrixAll = do.call(cbind,lapply(finalNamesCPM, function(currCellType){
  apply((singleCellDataFinalRefined[,selectedCellsOverall])[,cellTypesRefined[selectedCellsOverall]==currCellType],1,mean)
}))
colnames(signatureMatrixAll) = finalNamesCPM
signatureMatrixAll = signatureMatrixAll[,colnames(cpmResFrac)]
dataComposition = calculateResiduals(GEDataCombatAll, cpmResDeconv)

#Figure 5B
middlePoint = 0.5
createVolcanoPlot = function(data, currTraj, middlePoint, cutOff,textSize, volcModel = NULL,toShowText=T,minValueX = NULL, maxValueX = NULL){
  dataForVolcanoCell = as.data.frame(t(apply(data, 2, function(x){
    if(length(unique(x[currTraj<middlePoint]))==1 & length(unique(x[currTraj>middlePoint]))==1){
      return(c(0,0))
    }
    p = -log(t.test(x[currTraj<middlePoint],x[currTraj>middlePoint])$p.value,base = 10)
    fold = log(mean(x[currTraj>middlePoint])/mean(x[currTraj<middlePoint]))
    c(p,fold)
  })))
  colnames(dataForVolcanoCell) = c("minusLogP", "logFC")
  signCells = p.adjust(10^-dataForVolcanoCell$minusLogP,method = "fdr")<cutOff
  if(is.null(volcModel)){
    cellVolcanoColour = rep(0, dim(dataForVolcanoCell)[1])
    cellVolcanoColour[signCells & dataForVolcanoCell$logFC>0] = 1
    cellVolcanoColour[signCells & dataForVolcanoCell$logFC<0] = 2
  }else{
    cellVolcanoColour = volcModel$plot$plot_env$cellVolcanoColour
  }
  dataForVolcanoCell$col = cellVolcanoColour
  cellNames = row.names(dataForVolcanoCell)
  cellNames[!signCells] = ""
  if(!toShowText){
    cellNames = rep("",length(cellNames))
  }
  dataForVolcanoCell$name = cellNames
  minValueForPlotX = min(dataForVolcanoCell$logFC)
  maxValueForPlotX = max(dataForVolcanoCell$logFC)
  if(!is.null(minValueX)){
    minValueForPlotX = minValueX
  }
  if(!is.null(maxValueX)){
    maxValueForPlotX = maxValueX
  }
  maxP = max(dataForVolcanoCell$minusLogP)
  if(!is.null(volcModel)){
    maxP = max(volcModel$plot$plot_env$dataForVolcanoCell$minusLogP)
  }
  dataForVolcanoCell
  currPlot = ggplot(data = dataForVolcanoCell, aes(x = logFC, y = minusLogP, colour = factor(col)))+geom_point(size = 2)+geom_vline(xintercept = 0, linetype = "dashed")+
    geom_text(aes(label = name), vjust = -1, size = textSize)+ xlim(minValueForPlotX,maxValueForPlotX)+ylim(0,maxP)+
    theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
  list(plot = currPlot, data = dataForVolcanoCell)
}
volcMin = -1
volcMax = 1
geneVolcano = createVolcanoPlot(t(GEDataCombatAll[row.names(dataComposition),]),trajToUseProfileToUse,middlePoint, 0.01,2,toShowText = F,minValueX = volcMin, maxValueX = volcMax)
geneVolcanoComposition = createVolcanoPlot(t(dataComposition),trajToUseProfileToUse,middlePoint, 0.01,2,volcModel=geneVolcano,toShowText = F,minValueX = -1, maxValueX = 1)

newCutOff = -log(0.05/dim(geneVolcano$data)[1],base=10)
geneColour = rep(0, dim(geneVolcano$data)[1])
geneColour[geneVolcano$data$minusLogP>newCutOff & geneVolcanoComposition$data$minusLogP<newCutOff] = 1
geneColour[geneVolcano$data$minusLogP<newCutOff & geneVolcanoComposition$data$minusLogP>newCutOff] = 2
geneColour[geneVolcano$data$minusLogP>newCutOff & geneVolcanoComposition$data$minusLogP>newCutOff] = 3
ggplot(data = NULL, aes(x = geneVolcano$data$minusLogP, y = geneVolcanoComposition$data$minusLogP, colour = factor(geneColour)))+geom_point()+
  geom_hline(yintercept = newCutOff, linetype = "dashed")+geom_vline(xintercept = newCutOff, linetype = "dashed")+
  scale_color_manual(values = c("0" = "lightgray","1" = "darkgray","2" = "purple", "3" = "orange"))+
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

#Figure 5C
corGeneTrajComposition = cor(t(dataComposition),trajToUseProfileToUse)[,1]
geneCorMatrixWithComp = cbind(corGeneTraj[names(corGeneTrajComposition)],corGeneTrajComposition)
colnames(geneCorMatrixWithComp) = c("Cohort1", "Cohort1_Cell")

finalGeneCorJustCell = geneCorMatrixWithComp[which(abs(geneCorMatrixWithComp[,1])<0.1 & abs(geneCorMatrixWithComp[,2])>0.5),]
finalGeneCorNotCell = geneCorMatrixWithComp[which(abs(geneCorMatrixWithComp[,1])>0.5 & abs(geneCorMatrixWithComp[,2])<0.1),]
finalGeneCorOppPos = geneCorMatrixWithComp[which(geneCorMatrixWithComp[,1]>0.4 & geneCorMatrixWithComp[,2]<0),]
finalGeneCorOppNeg = geneCorMatrixWithComp[which(geneCorMatrixWithComp[,1]<(-0.4) & geneCorMatrixWithComp[,2]>0),]

geneCorForFigure = rbind(finalGeneCorNotCell[order(finalGeneCorNotCell[,1]),],finalGeneCorOppNeg[order(finalGeneCorOppNeg[,2]),])

#Figure 5D
source(file.path(folder, "enrichmentTest.R"))
geneSets = readRDS(file.path(folder, "geneSets.rds"))
geneSetsList = apply(geneSets,1,function(x){unlist(strsplit(as.character(as.matrix(x)),","))})

enrichmentOriginal = runEnrichment(corGeneTraj[names(corGeneTrajComposition)],geneSetsList,justHallmark = F)
enrichmentComp = runEnrichment(corGeneTrajComposition,geneSetsList,justHallmark = F)

enrichScoresBefore = abs(enrichmentOriginal)
enrichScoresAfter = abs(enrichmentComp)

enrichmentCompDF = as.data.frame(cbind(enrichScoresBefore,enrichScoresAfter))
colnames(enrichmentCompDF) = c("All","CellEffect")
enrichmentCompDF$name = row.names(enrichmentCompDF)

finalEnrichedMatrix = enrichmentCompDF[which(!is.na(enrichmentCompDF$CellEffect)),]

#Figure 5E
genesForIFNNFkappaComp = c("NFKB1","NFKB2","REL","RELA","RELB","IFIT1","IFIT2","IFIT3","IFIT5")
matrixForIFNNFkappaComp = geneCorMatrixWithComp[genesForIFNNFkappaComp,]

#Figure 5F
corGeneTrajDisrupted = cor(t(GEDataCombatAll),profilePredictionsDisrupted$predictions)[,1]
geneCorresidualMatrixWithDisrupted = cbind(corGeneTraj, corGeneTrajDisrupted[names(corGeneTraj)])
colnames(geneCorresidualMatrixWithDisrupted) = c("Cohort1","Cohort1_disrupted")

cutOffBadGroup = 0.2
cutOffBaseLvl = 0.4
lostAssociation = names(which(abs(geneCorresidualMatrixWithDisrupted[,1])>cutOffBaseLvl & abs(geneCorresidualMatrixWithDisrupted[,2])<cutOffBadGroup))
gainAssociation = names(which(abs(geneCorresidualMatrixWithDisrupted[,1])<cutOffBadGroup & abs(geneCorresidualMatrixWithDisrupted[,2])>cutOffBaseLvl))