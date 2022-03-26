
### Metabric data

pen.z <- seq(0.01,0.05,by = 0.01)
pen.x <- seq(0.01,0.11,by = 0.02)
cor1 <- list()
for(i in 1:length(pen.z)){
  print(i)
  print(pen.z[i])
  c1 <- NULL
  system.time(for(j in 1:length(pen.x)){
    print(i)
    print(pen.z[i])
    
    print(j)
    print(pen.x[j])
    
    s1 = CCA(x = cdna.val.mod,z = exp.val.mod.norm,trace = T,K = 20,penaltyz = pen.z[i],penaltyx = pen.x[j])
    c1 <- cbind(c1,s1$cors)
  })
  cor1[[i]] <- c1
  
}

g <- cn <- list()
for(i in 1:14){
  g[[i]] <- which(s.alt$v[,i] !=0)
}

for(i in 1:14){
  cn[[i]] <- which(s.opt$u[,i] !=0)
}

log.reg <- list()

sig.genes <- sig.genes.fdr <- sig.genes.bonf <- list()
for(j in 1:6){
  sg <- sg1 <- sg2 <- NULL 
  log.reg <- list()
  for(i in 1:14){
    
    log.reg[[i]] <- glm(y.mat[,j] ~ exp.val.mod.norm[,g[[i]]],family = "binomial")
    s.lr <- summary(log.reg[[i]])
    pv.lr <- s.lr$coefficients[-1,4]
    pv.fdr <- p.adjust(pv.lr,method = "BH")
    sg <- c(sg,colnames(exp.val.mod.norm[,g[[i]]])[which(s.lr$coefficients[-1,4] < 2.5e-06)])
    sg1 <- c(sg1,colnames(exp.val.mod.norm[,g[[i]]])[which(pv.fdr < 0.05)])
    sg2 <- c(sg2,colnames(exp.val.mod.norm[,g[[i]]])[which(s.lr$coefficients[-1,4] < 0.05/900)])
  }
  sig.genes[[j]] <- unique(sg)
  sig.genes.fdr[[j]] <- unique(sg1)
  sig.genes.bonf[[j]] <- unique(sg2)
}

sig.cauchy.bonf <- sig.cauchy.fdr <- list()
for(j in 1:14){
  print(j)
  sel.gene.x = exp.val.mod.norm[,which(s.opt$v[,j]!=0)]
  
  pv1 <- NULL
  for(i in 1:6){
    print(i)
    logis1 <- glm(y.mat[,i] ~ sel.gene.x,family = "binomial")
    pv1 <- cbind(pv1,summary(logis1)$coefficients[-1,4])
    
  }
  rownames(pv1) = colnames(sel.gene.x)
  colnames(pv1) = colnames(y.mat)
  
  cauchy1 = qcauchy(pv1,lower.tail=F)
  r.cauchy = rowMeans(cauchy1)
  p.cauchy = pcauchy(abs(r.cauchy),lower.tail=F)*2
  p.cauchy.fdr <- p.adjust(p.cauchy,method = "BH")
  sig.cauchy.bonf[[j]] <- names(which(p.cauchy < 0.05/length(p.cauchy)))
  sig.cauchy.fdr[[j]] <- names(which(p.cauchy.fdr < 0.05))
}
writeLines(names(which(p.cauchy < 0.05)),"./All_comp3.txt")

###BRCA

sig.surv <- list()
for(cols in 1:14){
  print(cols)
  sel.gene.x = exp.val.mod.norm[,which(s.opt$v[,cols]!=0)]
  surv1 = summary(coxph(Surv(clinical.matched$OS_MONTHS,dead.status) ~ sel.gene.x))
  surv.coeff = surv1$coefficients
  rownames(surv.coeff) = colnames(sel.gene.x)
  p.surv <- surv.coeff[,5]
  p.surv.fdr <- p.adjust(p.surv,method = "BH")
  # surv.coeff[which(surv.coeff[,5] < 0.05),]
  sig.surv[[cols]] <- names(p.surv.fdr[which(p.surv.fdr < 0.05)])
}

###ColCan

pen.z <- seq(0.05,0.1,by = 0.01)
pen.x <- seq(0.01,0.21,by = 0.02)
cor1 <- list()
for(i in 1:length(pen.z)){
  print(i)
  print(pen.z[i])
  c1 <- NULL
  system.time(for(j in 1:length(pen.x)){
    print(i)
    print(pen.z[i])
    
    print(j)
    print(pen.x[j])
    
    s1 = CCA(x = cdna.val,z = exp.val.norm,trace = T,K = 30,penaltyz = pen.z[i],penaltyx = pen.x[j])
    c1 <- cbind(c1,s1$cors)
  })
  cor1[[i]] <- c1
  
}






sig.surv <- sig.surv1 <- list()
for(cols in 1:14){
  print(cols)
  sel.gene.x = exp.val.mod.norm[,which(s.alt$v[,cols]!=0)]
  surv1 = summary(coxph(Surv(clinical.matched$OS_MONTHS,dead.status) ~ sel.gene.x))
  surv.coeff = surv1$coefficients
  rownames(surv.coeff) = colnames(sel.gene.x)
  p.surv <- surv.coeff[,5]
  p.surv.fdr <- p.adjust(p.surv,method = "BH")
  #surv.coeff[which(surv.coeff[,5] < 0.05),]
  sig.surv[[cols]] <- names(p.surv.fdr[which(p.surv.fdr < 0.05)])
  sig.surv1[[cols]] <- names(p.surv[which(p.surv < 2.5e-06)])
}
unique(unlist(sig.surv))
unique(unlist(sig.surv1))












pv.ea <- ci.ea <- odr.ea <- NULL
enr.ea <- NULL
for(i in 2:18){
  print(i)
  print(colnames(funct.no.eqtl.all)[i])
  f.ea = fisher.test(as.factor(funct.no.eqtl.all[[i]]),funct.no.eqtl.all$cond)
  pv.ea <- c(pv.ea,f.ea$p.value); odr.ea <- c(odr.ea,f.ea$estimate)
  ci.ea <- rbind(ci.ea,f.ea$conf.int);
  enr.ea <- cbind(odr.ea,pv.ea,ci.ea);
}



pv.aa <- ci.aa <- odr.aa <- NULL
enr.aa <- NULL
for(i in 2:18){
  print(i)
  print(colnames(funct.aa)[i])
  f.aa = fisher.test(as.factor(funct.aa[[i]]),funct.aa$cond)
  pv.aa <- c(pv.aa,f.aa$p.value); odr.aa <- c(odr.aa,f.aa$estimate)
  ci.aa <- rbind(ci.aa,f.aa$conf.int);
  enr.aa <- cbind(odr.aa,pv.aa,ci.aa);
}

apply(cor1[[1]],2,FUN = function(x){a1 = x;a2 = c(1,a1[-length(a1)]);c1 <- which.max(a1/a2);a <- mean(x[1:c1]); return(a)})


for(i in 1:14){
  print(i)
  writeLines(colnames(cdna.val.mod.norm)[which(s.alt$u[,i]!=0)],paste0("cdna_comp",i,".txt"))
  writeLines(colnames(exp.val.mod.norm)[which(s.alt$v[,i]!=0)],paste0("gene_comp",i,".txt"))
  
}

cx <- colSums(s.opt$v!=0)
sx.all <- NULL
for(j in 1:14){
  print(j)
  sx <- NULL
  for(k in 1:14){
    ax <- NULL
    for(i in 1:cx[j]){
      y <- cdna.val.mod[,which(s.opt$u[,j] !=0)[i]]
      ax <- cbind(ax,unlist(apply(exp.val.mod.norm[,which(s.opt$v[,k] !=0)],2,FUN = function(x){cor(x,y)})))
    }
    sx <- c(sx,mean(ax^2))
  }
  sx.all <- cbind(sx.all,sx)
}

####
sx.all <- NULL
for(j in 1:14){
  print(j)
  sx <- NULL
  for(k in 1:14){
    ax <- NULL
    for(i in 1:cx[j]){
      y <- cdna.val.mod.norm[,which(s.alt$u[,j] !=0)[i]]
      ax <- cbind(ax,unlist(apply(exp.val.mod.norm[,which(s.alt$v[,k] !=0)],2,FUN = function(x){cor(x,y)})))
    }
    sx <- c(sx,mean(ax^2))
  }
  sx.all <- cbind(sx.all,sx)
}



apply(ax^2,1,max)


rep.x <- NULL
for(i in 1:14){
  print(i)
  gx <- colnames(exp.val.mod.norm)[which(s.opt$v[,i] !=0)]
  cx <- colnames(cdna.val.mod)[which(s.opt$u[,i] !=0)]
  
  sub.tab.gx <- exp.tcga.val[,which(colnames(exp.tcga.val) %in% gx)]
  sub.tab.cx <- cdna.tcga.val[,which(colnames(cdna.tcga.val) %in% cx)]
  cor.x = cor(sub.tab.gx,sub.tab.cx,use = "pairwise.complete.obs")
  rep.x <- c(rep.x,mean(cor.x^2,na.rm = T))
}

rand.gx <- rand.cx <- rand.x <- NULL
for(i in 1:1000){
  print(i)
  sub.tab.gx.1 = exp.tcga.val[,sample(1:dim(exp.tcga.val)[2],size = dim(sub.tab.gx)[2],replace = F)]
  sub.tab.cx.1 = cdna.tcga.val[,sample(1:dim(cdna.tcga.val)[2],size = dim(sub.tab.cx)[2],replace = F)]
  cor.x.cx = cor(sub.tab.gx.1,sub.tab.cx,use = "pairwise.complete.obs")
  cor.x.gx = cor(sub.tab.gx,sub.tab.cx.1,use = "pairwise.complete.obs")
  cor.x.1 = cor(sub.tab.gx.1,sub.tab.cx.1,use = "pairwise.complete.obs")
  
  rand.gx <- c(rand.gx,mean(cor.x.gx^2,na.rm = T))
  rand.cx <- c(rand.cx,mean(cor.x.cx^2,na.rm = T))
  rand.x <- c(rand.x,mean(cor.x.1^2,na.rm = T))
  
}



for(i in 1:14){
  
  gx <- colnames(exp.val.mod.norm)[which(s.opt$v[,i] !=0)]
  cx <- colnames(cdna.val.mod)[which(s.opt$u[,i] !=0)]
  
  sub.tab.gx <- exp.tcga.val[,which(colnames(exp.tcga.val) %in% gx)]
  sub.tab.cx <- cdna.tcga.val[,which(colnames(cdna.tcga.val) %in% cx)]
  
  cutoff <- ceiling();
  tx <- table(colSums(tf.encode.data[which(tf.encode$GeneSym %in% colnames(sub.tab.gx)),]))
  sum()
  
  
}




sg <- sg1 <-log.reg <-  list()
ps <- NULL
for(i in 1:14){
  
  log.reg[[i]] <- glm(y.mat[,1] ~ exp.val.mod.norm[,g[[i]]],family = "binomial")
  s.lr <- summary(log.reg[[i]])
  pv.lr <- s.lr$coefficients[-1,4]
  pv.fdr <- p.adjust(pv.lr,method = "BH")
  ps <- cbind(ps,as.vector(pv.lr))
  sg[[i]] <- colnames(exp.val.mod.norm[,g[[i]]])[which(s.lr$coefficients[-1,4] < 2.5e-06)]
  sg1[[i]] <-colnames(exp.val.mod.norm[,g[[i]]])[which(pv.fdr < 0.05)]
}



sg <- sg1 <- log.reg <-  list()
ps <- NULL
for(i in 1:14){
  
  log.reg[[i]] <- glm(y.mat[,3] ~ exp.val.mod.norm[,g[[i]]],family = "binomial")
  s.lr <- summary(log.reg[[i]])
  pv.lr <- s.lr$coefficients[-1,4]
  pv.fdr <- p.adjust(pv.lr,method = "BH")
  sg[[i]] <- colnames(exp.val.mod.norm[,g[[i]]])[which(s.lr$coefficients[-1,4] < 2.5e-06)]
  sg1[[i]] <-colnames(exp.val.mod.norm[,g[[i]]])[which(pv.fdr < 0.05)]
}

writeLines(Reduce(union,sg),"her2_bonf.txt")
writeLines(Reduce(union,sg1),"her2_fdr.txt")


pv.list <- NULL
for(i in 1:14){
  print(i)
  log.reg <- list()
  pvs <- NULL
  for(j in 1:6){
    log.reg[[i]] <- glm(y.mat[,j] ~ exp.val.mod.norm[,g[[i]]],family = "binomial")
    s.lr <- summary(log.reg[[i]])
    pv.lr <- s.lr$coefficients[-1,4]
    pvs <- cbind(pvs,pv.lr)
  }
  surv1 = summary(coxph(Surv(clinical.matched$OS_MONTHS,dead.status) ~ exp.val.mod.norm[,g[[i]]]))
  surv.coeff = surv1$coefficients
  p.surv <- surv.coeff[,5]
  pvs <- cbind(pvs,p.surv)
  rownames(pvs) <- colnames(exp.val.mod.norm)[g[[i]]]
  pv.list[[i]] <- pvs  
  
  
}

cauchy.test <- function(x){
  
  cauchy1 = apply(x,1,FUN = function(x){qcauchy(x,lower.tail=F)})
  r.cauchy = colMeans(cauchy1)
  p.cauchy = pcauchy(abs(r.cauchy),lower.tail=F)*2
  p.cauchy.fdr <- p.adjust(p.cauchy,method = "BH")
  return(list("pv" = p.cauchy,"FDR" = p.cauchy.fdr))
  
}



pv.all <- g.all <- g.bonf <- minp <- NULL
for(i in 1:14){
  
  print(i)
  cauchy1 = apply(pv.list[[i]],1,FUN = function(x){qcauchy(x,lower.tail=F)})
  r.cauchy = colMeans(cauchy1)
  p.cauchy = pcauchy(r.cauchy,lower.tail=F)
  p.cauchy.fdr <- p.adjust(p.cauchy,method = "BH")
  pv.all[[i]] <- cbind(apply(pv.list[[i]],1,min,na.rm = T),p.cauchy,p.cauchy.fdr)
  g.all[[i]] <- names(p.cauchy.fdr[which(p.cauchy.fdr < 0.05)])
  g.bonf[[i]] <- names(p.cauchy[which(p.cauchy < 2.5e-06)])
  minp[[i]] <- rownames(pv.all[[i]][which(pv.all[[i]][,1] < 2.5e-06)])
  
}
length(Reduce(union,g.all))
length(Reduce(union,g.bonf))
