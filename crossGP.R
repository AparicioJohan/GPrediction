library(tidyverse)
library(sommer)
library(asreml)
library(BGLR)
library(readr)


#geno <- "D:/OneDrive - CGIAR/2020/Documento presentacion 4/ADP_VEC_VEF/ADP_VEC_VEF_repMasked_q40_s_fi_maf1%_oh1%_i900_imputed_rrBLUP.in"
#samp <- "D:/OneDrive - CGIAR/2020/Documento presentacion 4/ADP_VEC_VEF/ADP_VEC_VEF_repMasked_q40_s_fi_maf1%_oh1%_i900_imputed_rrBLUP_samples.txt"
#phen <- "D:/OneDrive - CGIAR/2020/Documento presentacion 4/Phenotypic Analysis.csv"

#prior <- c("ASReml","RKHS","sommer")
#traits <- c("Dar16C_hiP", "Dar16C_loP")

"crossGP" <- function(geno,samp,phen,prior, niter=50, testporc=0.3, traits=NULL){
  
  # Samples 
  samp = read.delim(samp, header = F)[,1]
  
  # Genotypic
  G = read.table(geno, row.names = as.character(samp), header = F)
  G <- G[ order(row.names(G)) , ]
  
  # phenotypic
  phen <- read.csv(phen) %>% data.frame()
  genoname <<- names(phen)[1]
  phen <- arrange(phen, get(genoname))
  
  # traits
  if(is.null(traits)) traits <- names(phen)[names(phen)!=genoname]
  
  if (sum(traits%in%names(phen))!=length(traits)) {
    stop("There are missings traits")
  } else {
    message("[]==============================================================[]")
    message("[]==================== Genomic Prediction  =====================[]")
    message("[]============= ASReml - BGLR - sommer package  ================[]")
    message("[]======= Last update: 2020-03-10  Johan Aparicio ==============[]")
    message("[]==============================================================[]\n")
  }
  
  Gen <- list()
  for (i in traits) {
    LinesA <- as.character(phen[!is.na(phen[,i]),genoname])  # A 
    LinesB <- as.character(rownames(G))                      # B
    Gen[[i]] <- intersect(LinesA,LinesB)
    message(i, "\t == " ,length(intersect(LinesA,LinesB)))
  }
  
  
  # Dataframe with results --------------------------------------------------
  out_table = data.frame(prior = as.character(), trait = as.character(), 
                         randomPop = as.character(), hrk=as.numeric(), 
                         corr = as.numeric(), finishedAt = as.character())
  
  message("\nprior", "\ttrait", "\trandomPop", "\thrk", "\tcorr" , "\tfinishedAt" )
  #--------------------------------------------------------------------------
  
  modMar <- c("BRR","BayesA","BayesB","BayesC","BLasso")
  stI_list <- matrix(data=list(), nrow=length(traits), ncol=length(modMar), 
                     dimnames=list(traits, modMar))
  
  #--------------------------------------------------------------------------
  
  niter <- ifelse(testporc==0,yes = 1,no = niter )
  
  for (i in traits  ) {
    
    LinesA <- as.character(phen[!is.na(phen[,i]),genoname])  # A 
    LinesB <- as.character(rownames(G))                      # B
    GT <- G[rownames(G)%in%intersect(LinesA,LinesB),]   # Marker
    
    phen2 <- phen %>% subset(.[,1]%in%rownames(GT)) %>% droplevels() 
    phen2$var <- phen2[,i]
    phen2$level <- phen2[,genoname]
    
    #------ ainverse  ---------
    
    
    G_hat <- GT
    AHAT<-A.mat(G_hat)
    AHAT_blend<-0.98*AHAT+0.02*diag(x=1,nrow=nrow(G_hat),ncol=nrow(G_hat))
    det(AHAT_blend)
    Ainv.blend<-solve(AHAT_blend)
    Ainv.blend[1:4,1:4]
    
    source("https://raw.githubusercontent.com/AparicioJohan/GPrediction/master/full2sparse.R")
    
    rownames(Ainv.blend) <- colnames(Ainv.blend) <- NULL
    AHAT.inv.sparse<-full2sparse(Ainv.blend)  #lower diag sparse matrix (Ainv.bend)
    colnames(AHAT.inv.sparse)<-c('Row','Column','Ainverse')
    head(AHAT.inv.sparse)
    
    # Preparing Ginverse
    ahatinv<<-data.frame(AHAT.inv.sparse)
    attr(ahatinv,"rowNames")<<-as.character(phen2[,genoname]) # A vector with ID names in order
    attr(ahatinv,"colNames")<<-as.character(phen2[,genoname])
    attr(ahatinv,"INVERSE")<<-TRUE  # Very important!
    head(ahatinv)
    
    phen2[,genoname] <- as.factor(phen2[,genoname])
    
    #-----------------------
    
    y <- phen2$var
    n <- length(y)
    
    
    # asr / rkhs / som / BRR / bayA / bayB / bayC / BL
    hrk <- hrk2 <- hrk3 <- hrk3 <- hrk4 <- hrk5 <- hrk6 <- hrk7 <- hrk8 <- c()
    corrk <- corrk2 <- corrk3 <- corrk4 <- corrk5 <- corrk6 <- corrk7 <- corrk8 <- c()
    fm <- fm2 <- fm3 <- fm4 <- fm5 <- fm6 <- fm7 <- fm8 <- list()
    
    for (pop in 1:niter) {
      
      tst<-sample(1:n,size=round(n*testporc),replace=FALSE)
      
      yNA<<-y
      yNA[tst]<<-NA
      
      if (testporc==0) tst <- 1:n
      
      # ------ asreml 
      
      if("ASReml"%in%prior){
        
        fm <-asreml(fixed=yNA~1,
                    random=~vm(level,ahatinv),
                    workspace=128e06,na.action=na.method(y="include"),data=phen2, trace=F)  
        
        # hrk[pop] <- as.numeric(vpredict(fm , h2~V1/(V1+V2))[1])   
        
        Vary<-var(yNA,na.rm=TRUE)
        hrk[pop] <- round(1-summary(fm)$varcomp[2,1]/Vary,4)
        
        predGBLUPcv<-predict(fm,classify="level",sed=T)$pvals
        predGBLUPcv<-predGBLUPcv[,2]
        
        corrk[pop] <- cor(y[tst],predGBLUPcv[tst]) %>% round(.,4)
        message(prior[prior=="ASReml"],"\t",i,"\t",paste("Pop",pop,sep="_"), "\t" ,hrk[pop], "\t",corrk[pop], "\t" , paste0(Sys.time()))
        
      }
      
      # ------ RKHS 
      
      if("RKHS"%in%prior){
        
        
        D = as.matrix( dist( GT, method="euclidean")) ^2
        D = D / mean(D)
        
        h = 0.5                                # fixed = list(~coord,data=phen2,model='FIXED')
        K = exp(-h * D)                        # fixed = list(~factor(Haplotype),data=phen2,model='FIXED')
        ETA = list(list(model="RKHS",   K=K))  # fixed = list(~factor(Meso.or.Andean),data=phen2,model='FIXED')
        
        
        fm2 = BGLR( y=yNA, ETA=ETA, 
                    nIter=10000, burnIn=1000, thin=5, 
                    verbose=FALSE)
        
        # hrk2[pop] <- fm2$ETA[[1]]$varU/(fm2$ETA[[1]]$varU + fm2$varE)
        
        Vary<-var(yNA,na.rm=TRUE)
        VarE<-fm2$varE
        hrk2[pop]<-round(1-VarE/Vary, 4)
        
        
        corrk2[pop] <- cor(y[tst],fm2$yHat[tst]) %>% round(.,4)
        message(prior[prior=="RKHS"],"\t",i,"\t",paste("Pop",pop,sep="_"), "\t" ,hrk2[pop], "\t",corrk2[pop], "\t" , paste0(Sys.time()))
        
      }
      
      # ------ sommer 
      
      if("sommer"%in%prior){
        
        phen2$vSom <- yNA
        
        fm3 = mmer(vSom~1,
                   random=~vs(level,Gu=AHAT_blend),
                   rcov=~units,
                   data=phen2,verbose=FALSE )
        
        # hrk3[pop] <- as.numeric(pin(fm3, h2 ~ (V1) / ( V1+V2) )[1])
        
        Vary<-var(yNA,na.rm=TRUE)
        hrk3[pop] <- round(1-summary(fm3)$varcomp[2,1]/Vary, 4)
        
        fm3$U$`u:level`$vSom <- as.data.frame(fm3$U$`u:level`$vSom)
        rownames(fm3$U$`u:level`$vSom) <- gsub("level","",rownames(fm3$U$`u:level`$vSom))
        
        corrk3[pop] <- cor(fm3$U$`u:level`$vSom[tst,],y[tst], use="complete") %>% round(.,4) 
        
        message(prior[prior=="sommer"],"\t",i,"\t",paste("Pop",pop,sep="_"), "\t" ,hrk3[pop], "\t",corrk3[pop], "\t" , paste0(Sys.time()))
        
      }
      
      # ------ BRR 
      
      if("BRR"%in%prior){
        
        ETA = list(list(model="BRR",    X=GT))  
        
        fm4 = BGLR( y=yNA, ETA=ETA, 
                    nIter=10000, burnIn=1000, thin=5, 
                    verbose=FALSE)
        
        Vary<-var(yNA,na.rm=TRUE)
        VarE<-fm4$varE
        hrk4[pop]<-round(1-VarE/Vary, 4)
        
        
        corrk4[pop] <- cor(y[tst],fm4$yHat[tst]) %>% round(.,4)
        message(prior[prior=="BRR"],"\t",i,"\t",paste("Pop",pop,sep="_"), "\t" ,hrk4[pop], "\t",corrk4[pop], "\t" , paste0(Sys.time()))
        
      }
      
      # ------ BayesA
      
      if("BayesA"%in%prior){
        
        ETA = list(list(model="BayesA",   X=GT))  
        
        fm5 = BGLR( y=yNA, ETA=ETA, 
                    nIter=10000, burnIn=1000, thin=5, 
                    verbose=FALSE)
        
        Vary<-var(yNA,na.rm=TRUE)
        VarE<-fm5$varE
        hrk5[pop]<-round(1-VarE/Vary, 4)
        
        
        corrk5[pop] <- cor(y[tst],fm5$yHat[tst]) %>% round(.,4)
        message(prior[prior=="BayesA"],"\t",i,"\t",paste("Pop",pop,sep="_"), "\t" ,hrk5[pop], "\t",corrk5[pop], "\t" , paste0(Sys.time()))
        
      }
      
      # ------ BayesB
      
      if("BayesB"%in%prior){
        
        ETA = list(list(model="BayesB",   X=GT))  
        
        fm6 = BGLR( y=yNA, ETA=ETA, 
                    nIter=10000, burnIn=1000, thin=5, 
                    verbose=FALSE)
        
        Vary<-var(yNA,na.rm=TRUE)
        VarE<-fm6$varE
        hrk6[pop]<-round(1-VarE/Vary, 4)
        
        
        corrk6[pop] <- cor(y[tst],fm6$yHat[tst]) %>% round(.,4)
        message(prior[prior=="BayesB"],"\t",i,"\t",paste("Pop",pop,sep="_"), "\t" ,hrk6[pop], "\t",corrk6[pop], "\t" , paste0(Sys.time()))
        
      }
      
      # ------ BayesC
      
      if("BayesC"%in%prior){
        
        ETA = list(list(model="BayesC",   X=GT))  
        
        fm7 = BGLR( y=yNA, ETA=ETA, 
                    nIter=10000, burnIn=1000, thin=5, 
                    verbose=FALSE)
        
        Vary<-var(yNA,na.rm=TRUE)
        VarE<-fm7$varE
        hrk7[pop]<-round(1-VarE/Vary, 4)
        
        
        corrk7[pop] <- cor(y[tst],fm7$yHat[tst]) %>% round(.,4)
        message(prior[prior=="BayesC"],"\t",i,"\t",paste("Pop",pop,sep="_"), "\t" ,hrk7[pop], "\t",corrk7[pop], "\t" , paste0(Sys.time()))
        
      }
      
      # ------ BLasso
      
      if("BLasso"%in%prior){
        
        ETA = list(list(model="BL",   X=GT))  
        
        fm8 = BGLR( y=yNA, ETA=ETA, 
                    nIter=10000, burnIn=1000, thin=5, 
                    verbose=FALSE)
        
        Vary<-var(yNA,na.rm=TRUE)
        VarE<-fm8$varE
        hrk8[pop]<-round(1-VarE/Vary, 4)
        
        
        corrk8[pop] <- cor(y[tst],fm8$yHat[tst]) %>% round(.,4)
        message(prior[prior=="BLasso"],"\t",i,"\t",paste("Pop",pop,sep="_"), "\t" ,hrk8[pop], "\t",corrk8[pop], "\t" , paste0(Sys.time()))
        
      }
      
      # final results ---------
      
      hf <- c(hrk[pop],hrk2[pop],hrk3[pop],hrk4[pop],hrk5[pop],hrk6[pop],hrk7[pop],hrk8[pop])
      corf <- c(corrk[pop],corrk2[pop],corrk3[pop],corrk4[pop],corrk5[pop],corrk6[pop],corrk7[pop],corrk8[pop])
      tmp = data.frame(prior = prior, trait = i, randomPop = paste0("Pop",pop), hrk = hf,  corr=corf , finishedAt = paste0(Sys.time()))
      out_table <- rbind(out_table,tmp)
      
    }
    
    
    if(testporc==0&sum(prior%in%modMar)>=1){
      stI_list[[i, 1]] <- fm4
      stI_list[[i, 2]] <- fm5
      stI_list[[i, 3]] <- fm6
      stI_list[[i, 4]] <- fm7
      stI_list[[i, 5]] <- fm8
    }
    
    
  }
  
  message("\n[]============================ End =============================[]\n") 
  
  return(list(data = out_table, 
              models = stI_list ))

}

## Example ##
# out_table <- crossGP(geno,samp,phen,prior,niter = 2,testporc = 0.3)

