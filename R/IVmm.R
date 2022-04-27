
IV_MM  <- function(formulaX,formulaY,random,subject,dataX,dataY,  binary=F, verbose=T){
  Wald2 <- function (Mod, pos = NULL, contrasts = NULL, name = NULL, value = NULL,seCorrec)
  {
    if (!(class(Mod) %in% c("hlme", "lcmm", "multlcmm",
                            "Jointlcmm")))
      stop("applies to \"hlme\" or \"lcmm\" or \"multlcmm\" or \"Jointlcmm\" objects only")
    if (inherits(Mod, "hlme") | inherits(Mod, "lcmm")) {
      nea <- sum(Mod$idea)
      nef <- Mod$N[2]
      nvc <- Mod$N[3]
      nprob <- Mod$N[1]
      idiag <- ifelse(Mod$idiag == 1, TRUE, FALSE)
    }
    if (inherits(Mod, "multlcmm")) {
      nea <- sum(Mod$idea0)
      nef <- Mod$N[3]
      nvc <- Mod$N[4]
      nprob <- 0
      idiag <- ifelse(Mod$idiag == 1, TRUE, FALSE)
    }
    if (inherits(Mod, "Jointlcmm")) {
      nea <- sum(Mod$idea)
      nef <- Mod$N[4]
      nvc <- Mod$N[5]
      nprob <- sum(Mod$N[1:3])
      idiag <- ifelse(Mod$idiag == 1, TRUE, FALSE)
    }
    if (nvc > 0) {
      debut <- nprob + nef + 1
      fin <- nprob + nef + nvc
      cholesky <- Mod$cholesky
      if (inherits(Mod, "multlcmm"))
        cholesky[1] <- NA
      if (isTRUE(idiag))
        cholesky[setdiff(1:(nea * (nea + 1)/2), 1:nea * (1:nea +
                                                           1)/2)] <- NA
      Mod$best[debut:fin] <- na.omit(cholesky)
    }
    l <- length(Mod$best)
    V <- matrix(0, nrow = l, ncol = l)
    V[upper.tri(V, diag = TRUE)] <- Mod$V
    V[lower.tri(V, diag = FALSE)] <- t(V)[lower.tri(V, diag = FALSE)]
    if (is.null(pos)) {
      stop("pos must be specified")
    }
    else {
      if (!is.vector(pos))
        stop("Error : pos must be a numeric vector")
      Mat <- matrix(0, nrow = length(pos), ncol = length(pos))
      #Mat <- V[pos, pos]
      Mat <- seCorrec[pos,pos]
      Vect <- Mod$best[pos]

      if (is.null(contrasts)) {
        if (!is.null(value)) {
          if (!is.vector(value))
            stop("Error : value must be a numeric vector")
          if (length(value) != length(pos))
            stop("value must have the same length as the vector pos")
          Vect <- Mod$best[pos] - value
        }
        Wald <- t(Vect) %*% solve(Mat) %*% Vect
        ddl <- length(pos)
        p_value <- 1 - pchisq(Wald, df = ddl)
        Results <- matrix(NA, nrow = 1, ncol = 2)
        colnames(Results) <- c("Wald Test", "p_value")
        if (is.null(name)) {
          if (!is.null(value)) {
            rownames(Results) <- paste(names(Mod$best[pos]),
                                       " = ", value, collapse = " and ",
                                       sep = "")
          }
          else {
            rownames(Results) <- paste(paste(names(Mod$best[pos]),
                                             collapse = " = "), "= 0")
          }
        }
        else {
          rownames(Results) <- name
        }
        Results[, 1] <- round(Wald, 5)
        Results[, 2] <- round(p_value, 5)
      } else {
        if (length(contrasts) != length(pos)) {
          stop("contrasts must have the same length as the vector pos")
        }
        if (sum(abs(contrasts)) == 0) {
          stop("The absolute value of the sum of contratsts components must be different from 0")
        }
        Scalaire <- sum(Vect * contrasts)
        if (!is.null(value)) {
          if (!is.vector(value))
            stop("value must be a numeric vector")
          if (length(value) != 1)
            stop("value must be a vector with a unique argument")
          Scalaire <- sum(Vect * contrasts) - value
        }
        Var <- t(contrasts) %*% Mat %*% contrasts
        Wald <- Scalaire/sqrt(Var)
        p_value <- 2 * (1 - pnorm(abs(Wald)))
        Results <- matrix(NA, nrow = 1, ncol = 4)
        colnames(Results) <- c("coef", "Se",
                               "Wald Test", "p_value")
        if (is.null(name)) {
          if (is.null(value))
            value <- 0
          rownames(Results) <- paste(paste(names(Mod$best[pos]),
                                           "*", contrasts, collapse = " + "),
                                     "= ", value)
        }
        else {
          rownames(Results) <- name
        }
        Results[, 1] <- round(sum(Vect * contrasts), 5)
        Results[, 2] <- round(sqrt(Var), 5)
        Results[, 3] <- round(Wald, 5)
        Results[, 4] <- round(p_value, 5)
      }
      return(Results)
    }
  }



  library(lcmm)

  ptm <-proc.time() # je ne sais pas Ã  quoi Ã§a correspond
  if(verbose==TRUE) cat("Be patient, IV estimation is running ... \n")

  cl <- match.call() # On retourne ce qui a Ã©tÃ© appelÃ©
  #args <- as.list(match.call(hlme))[-1]

  #### Message d'erreur
  if(missing(random)) random <- ~-1
  if(class(random)!="formula") stop("The argument random must be a formula")
  if(class(formulaX)!="formula") stop("The argument formulaX must be a formula")
  if(class(formulaY)!="formula") stop("The argument formulaY must be a formula")
  if(missing(dataX)){ stop("The argument dataX should be specified and defined as a data.frame")}
  if(nrow(dataX)==0) stop("DataX should not be empty")
  if(missing(dataY)){ stop("The argument dataY should be specified and defined as a data.frame")}
  if(nrow(dataY)==0) stop("DataY should not be empty")
  if(missing(subject)){ stop("The argument subject must be specified in any model even without random-effects")}

  if(!is.numeric(dataY[[subject]])) stop("The argument subject must be numeric")


  ## TO DO RECUPERER LES DONNEES
  data_suj <- model.frame(formulaX, dataX, na.action=NULL) # Avoir un jeu de donnÃ©es uniquement avec les variables que l'on souhaite
  XforXe <- model.matrix(formulaX, data = data_suj)
  data_suj$ID <- dataX[,subject]
  expo <- colnames(data_suj)[1]

  dataY2 <- model.frame(formulaY, dataY, na.action=NULL) # Avoir un jeu de donnÃ©es uniquement avec les variables que l'on souhaite
  dataY2$ID <- dataY[,subject]


  ## TO DO FAIRE LA PREMIERE REGRESSION POUR PREDIRE
  ## TO DO FAIRE LA DEUXIEME REGRESSION (MODELE MIXTE AVEC HLME)
  ## CAS CONTINU & CAS BINAIRE
  if(binary == F){
    regLM <- lm(formulaX, data=data_suj)
    #data_suj$predLM <- XforXe%*%coef(regLM)
    data_suj$predLM <- predict(regLM)
    dataY2 <- merge(dataY2,data_suj[,c("predLM","ID")], by="ID")
    expo <- paste("\\b",expo,"\\b",sep="")
    formula <- as.formula(gsub(expo,"predLM",deparse(formulaY)))
    regMM <- hlme(formula, random = random, subject = subject, data=dataY2, verbose = F)
    #regMM <- hlme(ISA ~ (Temps+I(Temps^2))*(predLM+niv2+niv3+niv4+niv5+SEX+AGE75)+Primo, random=~1+Temps+I(Temps^2), subject="ID",data=dataY2)
    #summaryMM <- summary(regMM)
    #data.frame(dataY2$predLM,dataY$PRED)
  }else{
    ## CAS BINAIRE
    regGLM <- glm(formulaX,data=data_suj, family="binomial")
    data_suj$predGLM <- 1/(1+exp(-(XforXe%*%coef(regGLM))))
    dataY2 <- merge(dataY2,data_suj[,c("predGLM","ID")], by="ID")
    expo <- paste("\\b",expo,"\\b",sep="")
    formula <- as.formula(gsub(expo,"predGLM",deparse(formulaY)))
    regMM <- hlme(formula, random = random, subject = subject, data=dataY2, verbose=F)
    #summaryMM <- summary(regMM)
  }

  ## TODO FAIRE LA CORRECTION DE LA STANDARD ERROR :
  # Faire tourner un modÃ¨le naif  puis un deuxiÃ¨me modÃ¨le naif
  #regMM_N2 <- hlme(formulaY, random=random, subject = subject, data= dataY2, verbose=F)

  nbEA <- regMM$N[3]+1
  #coef2SLS <- regMM$best[1:(length(regMM$best)-nbEA)]
  #binit <- regMM_N2$best
  #binit[1:(length(regMM$best)-nbEA)] <- coef2SLS


  regMM$best
  #regMM_N <- hlme(formulaY, random=random, subject = subject, data= dataY2, B = binit, posfix = 1:regMM$N[2], verbose=F)

  l<- sum(regMM$idea0)
  B <- matrix(0,nrow=l,ncol=l)
  B[upper.tri(B,diag=TRUE)] <- regMM$best[(length(regMM$best)-nbEA+1):(length(regMM$best)-1)]
  B[lower.tri(B,diag=FALSE)] <- t(B)[lower.tri(B,diag=FALSE)]
  sigma2 <- (regMM$best[length(regMM$best)])**2


  Vbeta <- matrix(data=0,nrow=length(regMM$best)-nbEA,ncol=regMM$N[2])
  Vrond <- matrix(data=0,nrow=length(regMM$best)-nbEA,ncol=regMM$N[2])


  for(i in unique(dataY2$ID)){
    #Variance de Yi
    Z <- model.matrix(random, dataY2)
    Z <- cbind(Z, dataY2$ID)
    Z <- matrix(Z[which(Z[,dim(Z)[2]]==i),-dim(Z)[2]],ncol=ncol(B))
    Vi <-as.matrix(Z)%*%B%*%t(as.matrix(Z))+diag(sigma2,nrow =length(which(dataY2$ID==i)))
    Vi1 <- solve(Vi)

    #Matrice des variables explicatives

    options(na.action='na.pass')
    X <- model.matrix(formula,dataY2)
    X <- cbind(X, dataY2$ID)
    X <- matrix(X[which(X[,dim(X)[2]]==i),-dim(X)[2]], ncol=regMM$N[2])
    Var <- as.matrix(t(X))%*%Vi1%*%as.matrix(X)

    # Calcul de la covariance :
    Xtrue <- cbind("int"=rep(1,k),dataY2[which(dataY2$ID==i),3],dataY2[which(dataY2$ID==i),colnames(data_suj)[1]],dataY2[which(dataY2$ID==i),colnames(data_suj)[1]]*dataY2[which(dataY2$ID==i),3]) # A REVOIR
    colnames(Xtrue) <- c("int","Temps","Xe","tXe")
    M <- as.matrix(Xtrue)%*%regMM$best[1:regMM$N[2]]
    Ma <- dataY2$Y[which(dataY2$ID==i)] - M
    covY <- Ma %*% t(Ma)
    VarCorr <- (as.matrix(t(X))%*%Vi1%*%covY%*%Vi1%*%as.matrix(X))

    Vbeta <- Vbeta+ Var
    Vrond <- Vrond + VarCorr


  }
  VarBeta <- solve(Vbeta)
  varcovNew0 <- sqrt(diag(VarBeta))

  var_corrig <- VarBeta %*% Vrond %*% VarBeta
  se_corr <- sqrt(diag(var_corrig))

  # Affichage
  coef <- regMM$best[1:regMM$N[2]]
  se <- se_corr
  wald <- NULL
  pvalue <- NULL
  for(i in 1:regMM$N[2]){
    wald <-append(wald,Wald2(regMM,pos=i,contrast=1, seCorrec = var_corrig)[3])
    pvalue <-append(pvalue,Wald2(regMM,pos=i, contrasts = 1, seCorrec = var_corrig)[4])
  }

  newsummary <- data.frame(coef,se,wald,pvalue)
  #newsummary[,2] <- varcovNew

  cost <- proc.time() - ptm
  if (verbose == TRUE){
    cat("The program took", round(cost[3], 2), "seconds \n")
  }
  return(list("call"=cl,"Fixed"=newsummary))
  # Recalculer P-value + wald
}

