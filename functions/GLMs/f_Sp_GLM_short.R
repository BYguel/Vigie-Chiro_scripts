test=F



#' a relatively generic function to fit multiple regression with glmmTMB, scaling numeric variables, supporting interactions and conversion of numeric variables to factors  

#' @param dataFile data file location

#' @param varInterest the name of the variable to be regressed (a character value)

#' @param listEffects the names of the explanatory variables (a character vector). Polynomial effects should be written as poly(x,n)

#' @param interactions a list of numeric vectors giving the position of the explanatory variables composing interactions (TO DO : adapt script to give variable names instead), default to NA (= no interactions)

#' @param formulaRandom the random part of the formula starting by "+" (default to ="+1" = no random effects)

#' @param selSample numeric, to downsample data (for testing)

#' @param tagModel a character tag identifying model outputs

#' @param family distrubution family of varInterest (default to "nbinom2", probably the better choice for abundance data)

#' @param asfactor a character vector giving the numeric variables to be treated as factor in the modelling, default to NA

#' @return write 5 files: (1) a .glm to save model fit in R format; (2) a "XXX_coefs.csv" table giving estimates and vif coefficients of the model; (3) a "XXX.log" to keep track of the formula of the model; (4) a "XXX_Res.csv" a table giving the residuals value; (5) a "forBackTransform_XXX.csv" table giving the mean and standard deviation value of numeric explanatory variables, to allow back transformation to real values when predicting

#' @example see at the end of this code

#' @param AutoCor: if TRUE, adds the autocorrelation structure starting by "+" with this form "ar1(varAutoCor + 1|groupingFact)" , by default is FALSE and TRUE only with no groupingFactor (i.e. grouping factor is the same for all data) 

#' @param  varAutoCor: the variable with autocorrelation by default is NA

Sp_GLM_short=function(dataFile,varInterest,listEffects,interactions=NA

                      ,formulaRandom="+1",AutoCor=FALSE,varAutoCor=NA, selSample=1e10,tagModel=""   ########## AJOUT DE AutoCor et varAutoCor

                      ,family="nbinom2",asfactor=NA)

{

  

  library(data.table)

  library(glmmTMB)

  library(plyr)

  library(beepr)

  library(corrplot)

  FAct=dataFile  # Variables à sélectionner et à tester en interaction  #### Emplacement du fichier avec les données

  VarAbondance=varInterest

  VarSimple=listEffects

  #Interactions=list(c(6,7,8),c(6,7,9),c(6,8,9),c(6,4))

  Interactions=interactions

  FormulaRandom=formulaRandom
  
  #FormulaRandom="+(1|espece)+(1|site)"

  #FormulaRandom="+(1|espece)"
  
  VarAutoCor=varAutoCor ### BEN MODIF
  
  #varAutoCor=annee 

  SelSample=selSample #for testing computing time

  #variables à rajouter : bioclim1 et 11, type de détecteur

  TagModel=tagModel

  # Famille

  familyMod=family

  # Modèle minimal

  #FormulaFix_TtSp="nb_contacts~(Jour+I(Jour^2)+I(Jour^3)+I(Jour^4)+I(Jour^5))*DecOT+(AT81+I(AT81^2))+((AT1+I(AT1^2))+(AT9+I(AT9^2)))+SpBioc12+SpHO1S+SpHO2S+SpHO4S+SpWS_S+SpWC_S"

  FormulaY=paste0(VarAbondance,"~1")  ##### ajoute 1 comme var explicative du modèle = modèle nul (si aucune variable ne sont ajoutés ensuite)

  FormulaXList=VarSimple

  

  #pour afficher les milisecondes

  op <- options(digits.secs=3)

  



  FormulaFix_TtSp=FormulaY #### Prends comme base la formule du modèle nul 

  for (i in 1:length(FormulaXList))  ####### Ajoute toutes les variables en effet additif (+) de la liste FormulaXList=VarSimple (qui est =listEffects dans les paramètres)

  {

    FormulaFix_TtSp=paste(FormulaFix_TtSp,FormulaXList[i],sep="+")

  }

  if(!is.na(Interactions)) ##### Ajoute les interactions à la formule 

  {

    for (i in 1:length(Interactions))

    {

      #    Intemp=paste(FormulaXList[Interactions[[i]][1]]

      #                ,FormulaXList[Interactions[[i]][2]],sep="*")

      Intemp=paste(FormulaXList[Interactions[[i]]],collapse="*") ### prends dans la liste de variable les interactions indiqués en paramètres

      

      FormulaFix_TtSp=paste(FormulaFix_TtSp,Intemp,sep="+")

    }

  }

  

  

  

  

  

  

  SpNuit=fread(FAct)  ########## charge les données indiqués en paramètres POTENTIELLEMENT A MODIF PARTOUT JUSTE POUR SAVOIR DE QUOI ON PARLE "spNuit" modif en "AbondData" donc les données d'abondance

  

  if(!is.na(asfactor))  #### transforme les variables indiqués en paramètres en facteur

  {

    SpNuit=as.data.frame(SpNuit)

    for (i in 1:length(asfactor))

    {

      test=match(asfactor[i],names(SpNuit)) #### match donne la position de la colonne de la variable i dans le datatable spNuit

      SpNuit[,test]=as.factor(SpNuit[,test])

    }

    SpNuit=as.data.table(SpNuit)

    

  }

  #compute summaries of activity

  ColA=match(VarAbondance,names(SpNuit))  ####  recupère le num de colonne de la = varInterest dans les arguments de la fonction 

  Ab=as.data.frame(SpNuit)[,ColA]  ### pour recuperer la colonne avec les abondances

  

  SpNuitwoNA=subset(SpNuit,!is.na(Ab))  #### pour supprimer les NA dans tout le jeu de données ou il n'y a pas d'abondance

  AbwoNA=subset(Ab,!is.na(Ab)) #### pour supprimer les NA dans les abondances

  

  

  SpA1=aggregate(AbwoNA,by=list(SpNuitwoNA$espece),FUN=mean) #### calcul les abondances moyennes par espèce (contient des 0 voir ci dessous)

  

    

    barplot(SpA1$x,names.arg=SpA1$Group.1,las=2,cex.names=0.6)  #### par sps (Group.1 est le nom de la colonne espece dans l'objet SpA1

    SpPos=subset(SpNuitwoNA,AbwoNA>0)  #### selectionner dans le tout le jeu de données sans NA les sps qui ont une abondance > 0

    AbPos=subset(Ab,Ab>0)  #### selectionner dans les abondances (avec NA) les sps avec une abondance > 0

    print(length(AbPos)) ### donne la taille du jeu de données sans NA et sans les abondance à 0

    if(length(AbPos)<=length(VarSimple))  #### regarde si le jeu de données est plus petit que le nombre de variables explicative choisies (VarSimple=listEffects)

    {

      print(paste(FAct,": too few positive data to fit model"))

    }else

    {

      

    SpOcc=aggregate(AbPos,by=list(SpPos$espece),FUN=length) #### donne la longueur de chaque jeu de données qd coupés par sps 

    barplot(SpOcc$x,names.arg=SpOcc$Group.1,las=2,cex.names=0.6)  ####

    

    SpAbIfP=aggregate(AbPos,by=list(SpPos$espece),FUN=mean)  ####  moyenne des abondances par sps dans le jeu de données sans 0 

    barplot(SpAbIfP$x,names.arg=SpAbIfP$Group.1,las=2,cex.names=0.6)

    

    ######### Calcul du VIF (adapté à glmmTMB, sinon il faut adapter v et nam)

    vif.mer <- function (fit) {

      ## adapted from rms::vif

      

      v <- vcov(fit)$cond

      nam <- names(fixef(fit)$cond)

      

      ## exclude intercepts

      ns <- sum(1 * (nam == "Intercept" | nam == "(Intercept)"))

      if (ns > 0) {

        v <- v[-(1:ns), -(1:ns), drop = FALSE]

        nam <- nam[-(1:ns)]

      }

      

      d <- diag(v)^0.5

      v <- diag(solve(v/(d %o% d)))

      names(v) <- nam

      v

    }

    

    

    ######## Pour correction autocorrelation spatiale

    #MaDataActiNew$FauxGroupe=rep(1,nrow(MaDataActiNew))

    #MaDataActiNew$Coord=numFactor(MaDataActiNew$X,MaDataActiNew$Y)

    

    

    SpNuit_SLPAGN=as.data.frame(SpNuitwoNA)   ######  renomme l'objet contenant toutes les données sans NA  dans les abondances

    

    OtherVariables=subset(names(SpNuit),!(names(SpNuit) %in% VarSimple)) ##### toutes les noms de colonnes pas utilisés dans le model comme variables predictrices

    

    SpNuit_Scale=subset(SpNuit_SLPAGN,select=OtherVariables)  ##### selectionne les colonnes pas utilisés dans le model comme variables predictrices

    
  
    ############## Recuperation et calcul de la moyenne et de l'ecart type  pour pouvoir faire la back transformation

    Mean=vector()

    Sdev=vector()

    VarList=vector()

    for (i in 1:length(VarSimple))

    {

      if(substr(VarSimple[i],1,5)=="poly(")   ### recherche les variables explicatives en polynome

      {

        Var=gsub("poly","",VarSimple[i]) ### supprime le terme "poly" de ces variables

        Terms=tstrsplit(Var,split=",") ### coupe la variable en polynome à la virgule ex: poly(var,2) ==> (var et 2) 

        VarTemp=substr(Terms[[1]],2,nchar(Terms[[1]])) ### sur la partie 1 de Terms avec le nom de la variable, récupère le nom de la variable (a partir du 2nd element et jusqu'au dernier)

      }else{

        VarTemp=VarSimple[i]

      }

      VarList=c(VarList,VarTemp)  ### créé un vecteur avec le nom de toutes les variables explicatives (mais pas comme declaré car sans les poly() par exemple) 

      Vinit=(SpNuit_SLPAGN)[,VarTemp]

      if(is.numeric(Vinit))

      {

        

        Vscale=scale(Vinit)  ### centre et reduit les variables 

        Mean=c(Mean,mean(Vinit)) ### calcul les moyennes des variables avant transfo et les ajoutes dans objet Mean

        Sdev=c(Sdev,sd(Vinit))  ### calcul l ecart type

      }else{

        Vscale=Vinit

        Mean=c(Mean,NA)

        Sdev=c(Sdev,NA)

      }

      SpNuit_Scale=cbind(SpNuit_Scale,Vscale) #### créé le jeu de données avec les variables rescale

if (exists(VarTemp,SpNuit_Scale)) ##### BEN MODIF pour eviter les doublons de variables mis en polynome 
{nline=match(VarTemp,names(SpNuit_Scale)) #### BEN MODIF 
SpNuit_Scale[,-nline]}##### BEN MODIF

      names(SpNuit_Scale)[ncol(SpNuit_Scale)]=VarTemp  #### renommer chaque variable nouvellement ajouté (ligne ci dessus) à SpNuit_Scale	

      if(i%%10==1){print(paste(i,Sys.time()))}  #### tous les 1er chiffre ex: 1, 11, 21, 31  R renvoie ce qu'il y a dans {}  ici un print 
#######  (e.g. use if (i %% 10 == 0) { #do something} to do something every 10th iteration)
      

    }   ################################### !!!!!!!!! UNE ERREUR DANS LA LOUPE CI DESSUS CAR MET 2 FOIS LES VARIABLES EN POLY AVEC LE MEME NOM  !!!!!!!!!

###browser()

    forBackTransform=data.frame(cbind(VarList,Mean,Sdev)) 

    fwrite(forBackTransform,paste0("./GLMs/forBackTransform/forBackTransform_"

                                   ,TagModel,".csv"))

    

    ColNumTest=unlist(lapply(SpNuit_Scale[1,],FUN=function(x) is.numeric(x)))  #### pour indiquer quelle variable est du numeric ou pas

    ColNum=subset(names(SpNuit_Scale),ColNumTest) ### recupere les noms de colonnes des variables numeric 

    SpNuit_ColNum=subset(SpNuit_Scale,select=ColNum) ### recupere les colonnes numeric

    MatCor=cor(SpNuit_ColNum) #### calcul la matrice de correlation

    corrplot(MatCor)
	
	
	VarAutoCor = subset(SpNuit_Scale,select=varAutoCor)#### BEN MODIF
	VarAutoCor = numFactor(VarAutoCor)  ####### transformation de la variable "portant" l'autocorrelation en numFactor() necessaire pour glmmTMB() BEN MODIF 

	FormulaAutoCor = paste0("+ar1(","VarAutoCor","- 1|groupingFact)")#### BEN MODIF 

if (AutoCor)#### BEN MODIF 
{
    Formula=as.formula(paste0(FormulaFix_TtSp#### BEN MODIF 

                              ,FormulaRandom,FormulaAutoCor))#### BEN MODIF 

}else{#### BEN MODIF 
    Formula=as.formula(paste0(FormulaFix_TtSp#### BEN MODIF 

                              ,FormulaRandom))#### BEN MODIF 
}#### BEN MODIF 
    

    

    if(SelSample<nrow(SpNuit_Scale))

    {

      SpNuit_Sample=SpNuit_Scale[sample.int(nrow(SpNuit_Scale),SelSample),]

    }else{

      SpNuit_Sample=SpNuit_Scale

    }

	n=nrow(SpNuit_Sample)  #### BEN MODIF
	groupingFact= factor(rep(1,n)) #### BEN MODIF


    Sys.time()

    ModSp=glmmTMB(Formula,data=SpNuit_Sample, family=familyMod)  #37 min

    Sys.time()

    beep()

    Res=residuals(ModSp)

    SpNuit_Sample$Res=Res

    

    Estimates=as.data.frame(coef(summary(ModSp))$cond)

    Estimates=cbind(term=row.names(Estimates),Estimates)

    save(ModSp,file=paste0("./GLMs/",TagModel,".glm"))

    VIFMod=c(1,vif.mer(ModSp))

    Estimates$VIF=VIFMod

    Suffix=tstrsplit(basename(FAct),split="[.]")[[1]]

    

    fwrite(Estimates,paste0("./GLMs/Summaries/",TagModel,"_",Suffix,"_Coefs.csv"),sep=";")

    

    fwrite(as.list(FormulaFix_TtSp),paste0("./GLMs/logs/",substr(Sys.time(),1,13),".log"))

    

    fwrite(SpNuit_Sample,paste0("./GLMs/",TagModel,"_",Suffix,"_Res.csv"))

  }

}

#for test

if(test)

{

  Sp_GLM_short(

    dataFile="./Pipkuh.csv"

    ,

    varInterest="abond"

    ,

    listEffects=c("annee","poly(julian,2)","sample_cat","nb_Tron_strict"

                  ,"temps_enr_strict","latitude","longitude","expansion_direct"

    )

    ,

    interactions=NA

    ,

    formulaRandom="+(1|site)"

    ,

    selSample=1e10

    ,

    tagModel="GLMalphatest_tendancesFY"

    ,

    family="nbinom2"

    ,

    asfactor="year"

  )

  

}
