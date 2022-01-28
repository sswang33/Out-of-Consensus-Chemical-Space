library(rcdk)
args <- commandArgs(trailingOnly = TRUE)
rules <- list(
  list("include"=c(791),"exclude"=c(34)),
  list("include"=c(34,286),"exclude"=c(4)),
  list("include"=c(212),"exclude"=c(34,791)),
  list("include"=c(140),"exclude"=c(34,791,212,242)),
  list("include"=c(471),"exclude"=c(34,791,212,140))
)
dat_path <- args[1]
dat <- scan(dat_path,what = c("character"),sep="\n")
header <- unlist(stringr::str_split(dat[1],"\t"))
dat <- as.data.frame(do.call("rbind",lapply(dat[-1],function(s){
  unlist(stringr::str_split(s,"\t"))
})),stringsAsFactors = F)
colnames(dat) <- header
dat$`SkinSensPred Score` <- as.numeric(dat$`SkinSensPred Score`)
mols <- try(rcdk::parse.smiles(dat$CanonicalSMILES))
if(class(mols)=="try-error"){
  print("smiles can not transfer")
}else{
  fp <- sapply(mols,function(m){
    get.fingerprint(m,type="graph")@bits
  })
  dat$`Out of Consensus Chemical Space` <- sapply(1:nrow(dat),function(i){
    if(dat$`SkinSensPred Score`[i]<=0.4){
      return(any(sapply(rules,function(rule){
        return(all(rule$include%in%fp[[i]])&&all(!(rule$exclude%in%fp[[i]])))
      })))
    }else{
      return("-")
    }
  })
}
write.table(dat,file="output.tsv",sep="\t",row.names = F)
