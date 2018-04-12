suppressWarnings(suppressPackageStartupMessages(require(optparse,quietly = T)))
suppressWarnings(suppressPackageStartupMessages(require(msPurity)))

# Get the parameter
option_list <- list(
  make_option(c("-i","--purity"),type="character"),
  make_option(c("-p","--ppm"),type="numeric"),
  make_option(c("-m","--mode"),type="character")
)
opt <- parse_args(OptionParser(option_list=option_list))

print(opt)

load(opt$purity)

grped_df <- pa@grped_df

msms <- pa@grped_ms2

puritydf <- pa@puritydf
print('filelist')
print(pa@fileList)
print('grped_df')
print(unique(grped_df$filename))
print('puritydf')
print(unique(puritydf$filename))




grped_df$fileid <- sapply(grped_df$filename, function(x) which(basename(pa@fileList)==x))

puritydf$fileid <- sapply(puritydf$filename, function(x) which(basename(pa@fileList)==x))

selfrag <- as.numeric(unique(grped_df$grpid)) 

of <- file(description = "outfile.msp", open = "w+a")

write.msp <- function(name,precmz,prectype,spectra,ofile){

    cat(paste("NAME: ", name, "\r\n", sep = ""), file = ofile)

    cat(paste("PRECURSORMZ: ", precmz , "\r\n", sep = ""), file = ofile)

    #cat(paste("PRECURSORTYPE: ", prectype, "\r\n", sep = ""), file = ofile) # No adducts? Annotation

    cat("Comment:\r\n", file = ofile)

    cat(paste("Num Peaks: ", nrow(spectra), "\r\n", sep = ""), file = ofile)

    cat(paste(paste(spectra[,1], spectra[,2], sep = "\t"), sep = "\r\n"), sep = "\r\n", file = ofile)

    cat("\r\n", file = ofile)
}

for(i in selfrag){

    j <- which(grped_df$grpid==i) 

    spec <- msms[[as.character(i)]]

    if (length(j)>1){

        grpd <- grped_df[j,] 

        if (opt$mode=="all"){ 

            for(jj in 1:length(j)){

                idj <- paste(i,jj,sep=".")

                specj <- spec[[jj]]

                grpdj <- grpd[jj,]

                name = paste(i, grpdj$sample, grpdj$pid, sep='-')

                write.msp(name,grpdj$precurMtchMZ,"",specj,of)
            }

        }else{

            prec_int <- sapply(grpd$precurMtchID, function(x) puritydf[ which(puritydf$seqNum==x & puritydf$fileid==grpd$fileid[1]), "precursorIntensity"] )

            if (opt$mode=="max"){

                idx <- which(prec_int==max(prec_int))

                grpd <- grpd[idx,]

                write.msp(i,grpd$precurMtchMZ,"",specj[[idx]],of)
            }

            if (opt$mode=="average"){

                file_weights <- prec_int/prec_int[which(prec_int==min(prec_int))] # spectra of the most intense precursor, largest weight

                merged_msms <- do.call("rbind", spec)

                colnames(merged_msms) <- c("mz","int")

                file_weights <-rep(file_weights,sapply(spec,nrow))

                merged_msms <- data.frame(cbind(merged_msms,file_weights))

                umzs <- sort(merged_msms$mz,index.return=T)$ix

                merged_msms <- merged_msms[umzs,]

                umzs <- merged_msms$mz

                mz_groups <- list() # mz windows to bind

                mz_grouped <- c() # used mzs

                ppm <- opt$ppm # ppm level to bind mzs

                for(y in 1:length(umzs)){ # spectra averaging
                    
                    z <- umzs[y]
                    
                    if(!(z %in% mz_grouped)){

                        mz_range <- z*(ppm/1e6)

                        mz_range <- c(z-mz_range,z+mz_range)

                        mz_group <- which(umzs>mz_range[1] & umzs<mz_range[2])

                        if(length(mz_group)>1){

                            zz <- umzs[max(mz_group)]
                    
                            mz_range <- zz*(ppm/1e6)
                    
                                mz_range <- c(zz-mz_range,zz+mz_range)

                                mz_group2 <- which(umzs>mz_range[1] & umzs<mz_range[2])

                                mz_group <- append(mz_group,mz_group2)
                        }

                        mz_group <- unique(mz_group)

                        mz_grouped <- append(mz_grouped,umzs[mz_group])

                        mz_group <- list(mz_group)

                        mz_groups <- append(mz_groups,mz_group)
                        }
                    }

                    averaged_spec <- t(sapply(mz_groups,function(x){
                    
                    if(length(x)==1){

                        mz <- merged_msms$mz[x]

                        nint <- merged_msms[x,]

                        nint <- nint$int*nint$file_weights

                    }else{

                        mz <- mean(merged_msms$mz[x])
                    
                        nint <- sum(sapply(x,function(y){
                          nint <- merged_msms[y,]
                          nint <- nint$int*nint$file_weights
                        }))
                    }
                        return(c(mz,nint))
                    }))
                    name = paste(i, grpd$sample, grpd$pid, sep='-')
                    write.msp(name,mean(grpd$precurMtchMZ),"",averaged_spec,of)

            }
        }
    }else{
        spec <- spec[[1]]

        grpd <- grped_df[j,]
        name = paste(i, grpd$sample, grpd$pid, sep='-')
        write.msp(name,grpd$precurMtchMZ,"",spec,of)
    }
}

