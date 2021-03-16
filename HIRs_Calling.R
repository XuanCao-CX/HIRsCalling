#!/usr/bin/Rscript --vanilla

library(optparse)
library(grid)
library(gridExtra)
library(RColorBrewer)
library(pheatmap)
library(data.table)


option_list = list(
    make_option(c("-i", "--input"), type="character", default=NULL,
    			help="Input merged square matrix file. File has header and rownames.",metavar="character"),
    make_option(c("-o", "--out_dir"), type="character", default=NULL,
                help="Output dir contain 4 files and 2 hetmap dir.\nFiles: HIRs.intra.txt, HIRs.inter.txt, Regions.intra.bdg, Regions.inter.bdg.\nHeatmap dir: heatmap_intra and heatmap_inter. Each heatmap_dir contain contact heatmap for each chr",metavar="character"),
    make_option(c("-r", "--regions"), type="character", default=NULL,
    			help="The genome regions corresponding to each row of matirx.Format: first 3 columns: chr,star,end.",metavar="character"),    
    make_option("--rm_diag_ratio", type="numeric", default=0.3,
              	help="For each intra matrix, it's the ratio of each chr bins number to remove the diagonal bins. Not used for inter chrs. Default: 0.3", metavar="character"),
    make_option("--quan_min_intra", type="numeric", default=0.9,
              	help="For each intra matrix, the minimum quantile ratio of each matrix to define the bins as signal bins to next step or noise bins. Default:0.9", metavar="character"),
    make_option("--quan_min_inter", type="numeric", default=0.9,
              	help="For each intra matrix, the minimum quantile ratio of each matrix to define the bins as signal bins to next step or noise bins. Default:0.9", metavar="character"),
    make_option("--continous_bins", type="integer", default=3,
          	  	help="For each matrix after removing the noise bins, we kept the continous bins that number >= continous_bins", metavar="character"),
    make_option("--inter_bin_apper", type="integer", default=3,
                help="For call inter HIRs of each chrs, we kept the continous bins that appaered number >= inter_bin_apper", metavar="character")
    )


opt_parser = OptionParser(usage="%prog: Call Pf High Interacted Regions(HIRs) from merged hic concact matrix.", option_list=option_list)
opt = parse_args(opt_parser)

####
read_matrix <- function(mat_file){
	### matrix file format: has header and raw.names
	### matrix row * column can be 10*10 or 10*20 or 20*10 and so on.
	mat=read.table(mat_file,header=TRUE,row.names=1,check.names=FALSE)
	mat=as.matrix(mat)
	return(mat)
}

read_binned_rgs <- function(region_file){
    ### region_file: chr, start, end, name. region_file has the same order with col_names.
    ### rgs(data.frame): chr, start, end, idx_in_matrix, name
    rgs=read.table(region_file,header=FALSE)
    #rgs[,1]=as.vector(rgs[,1])
    rgs=cbind(rgs[,1:3],bin_idx=seq(dim(rgs)[1]))
    #chrs=sort(unique(as.vector(rgs[,1])))
	return(rgs)
}

cen_rm_bins_fun <- function(cen_region, extend_dis, bin_size){
    ### region_file: chr, start, end. region_file has the same order with col_names.
    ### rgs(data.frame): chr, start, end, idx_in_matrix, name
    cen=read.table(cen_region,header=FALSE)
    cen_loci_bin=round(((cen[,3]-cen[,2])/2 + cen[,2])/10000)
    extend_bins=extend_dis/bin_size
    start_bin=cen_loci_bin -extend_bins
    end_bin=cen_loci_bin + extend_bins
    cen_rm_bins=cbind(data.frame(chr=as.character(cen[,1])),start_bin,end_bin)
  return(cen_rm_bins)
}


write_regions <- function(regions_to_write, out_path){
	write.table(regions_to_write, file= out_path, quote=FALSE,sep="\t", col.names = FALSE, row.names = FALSE)
}
###########

mat_scale <- function(mat,quan_max=0.95,quan_min=0.05){
	#### To define the max_value and min_vaule for ploting contact heatmap. For intra matrix, remove the diag 1 bin.
  	### For intra, we used: quan_max=0.9, quan_min=0.1,
  	### For inter, we used: quan_max=0.9, quna_min=0.85
	nx=dim(mat)[1]
  	ny=dim(mat)[2]
  	x=seq(nx)
  	if (nx == ny){for (i in x){mat[i,i]=0}}
	cut_value_max = quantile(mat[mat>0],quan_max)
	cut_value_min = quantile(mat[mat>0],quan_min)
	cut_value_max = as.vector(cut_value_max)
	cut_value_min = as.vector(cut_value_min)
	return(c(cut_value_max,cut_value_min))
}

plot_contact_mat<-function(mat,cut_value_max,cut_value_min,title_text="",col_text=colorRampPalette(c("white",brewer.pal(9,"Reds")))(100)){
  	if (max(mat) != 0){
		mat[mat>cut_value_max]=cut_value_max
  		mat[mat<cut_value_min]=cut_value_min
  		p=pheatmap(mat, cluster_cols=F,cluster_rows=F,
             show_colnames=F,show_rownames = F,
             border_color=NA,color = col_text,main=title_text)
  		}else{
  			mat[1,1]=1
  			p=pheatmap(mat, cluster_cols=F,cluster_rows=F,show_colnames=F,show_rownames = F,border_color=NA,color = colorRampPalette(c(brewer.pal(9,"Reds")[2]))(100),main=title_text)
  		}
  return(p[[4]])
}

remove_diag_values <- function(mat, remove_diag_bins_ratio=0.5){
	#### Remove bins ratio only for intra matirx
  n=dim(mat)[1]
  remove_diag_bins=round(n*remove_diag_bins_ratio)
  x=seq(n)
  for (i in x){
    j1=max(1,i-remove_diag_bins)
    j2=min(i+remove_diag_bins,n)
    mat[i,j1:j2]=0}
  return(mat)
}

remove_noise_mat <- function(mat,mat_min,continous_bins=4) {
	### mat_min : min value to filter
  	row_n = dim(mat)[1]
  	col_n=dim(mat)[2]
  	for (i in 1:row_n){
  		l = mat[i,]
    	idxes = which(l>mat_min)
    	idx_n = length(idxes)
    	idxes_select = c()
    	if (idx_n >= continous_bins){
      		idxes_temp=c(idxes[1])
      		for (j in 2:idx_n){
        		i1 = idxes_temp[length(idxes_temp)]
        		i2 = idxes[j]
        		if (i1+1 == i2){
        			idxes_temp = c(idxes_temp, i2)
        		}else{
          			if (length(idxes_temp) >= continous_bins){idxes_select = c(idxes_select, idxes_temp)}
          			idxes_temp = c(i2)
          		}
        		if (length(idxes_temp) >= continous_bins){
        			idxes_select = c(idxes_select, idxes_temp)
        		}
        		idxes_0 = setdiff(seq(col_n),idxes_select)
      		}
    	}else{
      		idxes_0 = seq(col_n) ######idxes for make value to zeros 
    	}
    	mat[i,idxes_0]=0
  	}
  return(mat)
}

remove_noise_oneChr <- function(mat,quan_min=0.9,continous_bins=3){
	mat_min=quantile(mat[mat>0],quan_min)
	mat_min=as.vector(mat_min)
  	mat[mat<mat_min]=mat_min

	### remover nosie three times
	times=5
	for (i in 1:times){
		#mat =as.data.frame(mat)
		mat = remove_noise_mat(mat, mat_min, continous_bins)
		mat = t(mat)
	}
	if (times%%2 == 1){mat=t(mat)}
	col_sum = as.vector(colSums(mat))
	return(list(mat, col_sum))	
}

call_HIRs_oneChr <- function(col_sum,chr,rgs){
	rgs_signal=cbind(rgs, col_sum=col_sum)
	###rg_signal:chr, star, end, bin_index, signal
	rgs_select=rgs_signal[which(rgs_signal$col_sum>0),]
	n=dim(rgs_select)[1]
	idx_s = 1
	idx_e = 1
	bin_s = rgs_select[idx_s,4]
	bin_e = bin_s
	HIRs = data.frame()
  rgs_select1 = data.frame()
	for (i in 2:n){
		bin_idx = rgs_select[i,4]
		bin_e = rgs_select[idx_e,4]
		if (bin_e +1 == bin_idx){
			idx_e = i
		}
		if (bin_e + 1 != bin_idx || i == n){
			d = idx_e - idx_s
			if (d > 0){
				chr = rgs_select[idx_s,1]
				start = rgs_select[idx_s,2]
				end = rgs_select[idx_e,3]
				value = mean(rgs_select[idx_s:idx_e,4])
				HIR = data.frame("chr"=chr,"start"=start,"end"=end,"value"=value)
				HIRs = rbind(HIRs,HIR)
        rgs_select1 = rbind(rgs_select1,rgs_select[idx_s:idx_e,])
			}
			idx_s = i
			idx_e = i
			bin_s = rgs_select[idx_s,4]
			bin_e = bin_s
		}
	}
	rgs_select1 = rgs_select1[,-4]
	return(list(HIRs,rgs_select1))
}

#####################
### The input file is merged matrix
#####################


###########
call_HIRs_oneChr_intra <- function(mat1, chr, rgs,remove_diag_bins_ratio=0.3, quan_min_call_HIRs=0.9,continous_bins=3, out_heatmap){
	### mat: merged matrix
    ### rgs: regions of one chr, format:chr, start, end, bin_index
    
    mat2=remove_diag_values(mat1,remove_diag_bins_ratio=remove_diag_bins_ratio)
    mat_rm_noise=remove_noise_oneChr(mat2,quan_min=quan_min_call_HIRs,continous_bins=continous_bins)
    mat3=mat_rm_noise[[1]]
    col_sum=mat_rm_noise[[2]]


    t_text=chr
    cut_values=mat_scale(mat1,quan_max=0.95,quan_min=0.05)
    cut_value_max=cut_values[1]
    cut_value_min=cut_values[2]

    p1=plot_contact_mat(mat1,cut_value_max,cut_value_min,title_text=t_text)
    p2=plot_contact_mat(mat2,cut_value_max,cut_value_min,title_text=t_text)
    ##p3=plot_contact_mat(mat3,cut_value_max,cut_value_min,title_text=t_text,col_text=colorRampPalette(c(brewer.pal(9,"Reds"))[2:9])(100))
    p3=plot_contact_mat(mat3,cut_value_max,cut_value_min,title_text=t_text)
    pdf(out_heatmap,width = 16,height = 5)
    grid.arrange(grobs=list(p1,p2,p3),ncol=3)
    dev.off()
    if (sum(col_sum) >0){
    	HIRs_oneChr = call_HIRs_oneChr(col_sum,chr,rgs)
    	}else{
    		HIRs_oneChr=c()
    	}
    return(list(list(p1,p2,p3),HIRs_oneChr))
}

########


call_HIRs_allChrs_intra <- function(mat_merged,chrs,rgs_all,remove_diag_bins_ratio=0.3, quan_min_call_HIRs=0.9,continous_bins=3, out_dir){
	heatmap_dir = paste(out_dir,"/heatmap_intra",sep="")
    dir.create(heatmap_dir)
    ps=list()
    HIRs=data.frame()
    bdg=data.frame()
    for (chr in chrs){
        idxes=which(rgs_all[,1]==chr)
        rgs=rgs_all[idxes,]
        x1=rgs[1,4]
        x2=rgs[dim(rgs)[1],4]
        mat1=mat_merged[x1:x2,x1:x2]
        out_heatmap = paste(heatmap_dir,"/",chr,'.pdf',sep="")
        HIRs_info=call_HIRs_oneChr_intra(mat1, chr, rgs,remove_diag_bins_ratio=remove_diag_bins_ratio, quan_min_call_HIRs=quan_min_call_HIRs,continous_bins=continous_bins, out_heatmap)
        ps=c(ps,HIRs_info[[1]])
        if (length(HIRs_info[[2]]) != 0){
            h = HIRs_info[[2]][[1]]
            r = HIRs_info[[2]][[2]]
            HIRs=rbind(HIRs,as.data.frame(h))
            bdg=rbind(bdg,as.data.frame(r))
        }
    }
    HIRs=HIRs[1:3]
    out_HIR=paste(out_dir,"/","HIRs.intra.bed",sep="")
    out_bdg=paste(out_dir,"/","Regions.intra.bdg",sep="")
    write_regions(HIRs, out_HIR)
    write_regions(bdg, out_bdg)
     
    out_heatmap_all=paste(heatmap_dir,"/",'1.all_chr.pdf',sep="")
    pdf(out_heatmap_all,width = 16,height = 70)
    grid.arrange(grobs=ps,ncol=3)
    dev.off()
} 

##################################
## Processing inter
##################################



call_HIRs_oneChr_inter <- function(mat_merged, chr, chrs,rgs_all, quan_min_call_HIRs=0.9,continous_bins=3, out_heatmap,inter_bin_apper=3){
  idxes=which(rgs_all[,1]==chr)
  rgs=rgs_all[idxes,]
  y1=rgs[1,4]
  y2=rgs[dim(rgs)[1],4]
  ps=list()
  i=1
  can_bin_num=rep(0,dim(rgs)[1])
  for (chr2 in chrs){
    if (chr2 != chr){
      idxes1=which(rgs_all[,1]==chr2)
      rgs1=rgs_all[idxes1,]
      x1 = rgs1[1,4]
      x2 = rgs1[dim(rgs1)[1],4]
      mat1 = mat_merged[x1:x2,y1:y2]
      mat_rm_noise=remove_noise_oneChr(mat1,quan_min=quan_min_call_HIRs,continous_bins=continous_bins)
      mat3=mat_rm_noise[[1]]
      col_sum=mat_rm_noise[[2]]
      can_bin_num[which(col_sum>0)]=can_bin_num[which(col_sum>0)]+1
      col_sum=as.data.frame(col_sum)
      if (i==1){
        col_sum_oneChr=col_sum
        }else{
        col_sum_oneChr=cbind(col_sum_oneChr,col_sum)}
      t_text=paste(chr,".vs.",chr2)
      cut_values=mat_scale(mat1,quan_max=0.95,quan_min=0.05)
      cut_value_max=cut_values[1]
      cut_value_min=cut_values[2]

      p1=plot_contact_mat(mat1,cut_value_max,cut_value_min,title_text=t_text)
      #p3=plot_contact_mat(mat3,cut_value_max,cut_value_min,title_text=t_text,col_text=colorRampPalette(c(brewer.pal(9,"Reds"))[2:9])(100))
      p3=plot_contact_mat(mat3,cut_value_max,cut_value_min,title_text=t_text)
      ps=c(ps,list(p1,p3))
      i=i+1
    } 
  }

  col_sum_oneChr=t(col_sum_oneChr)
  col_sum_oneChr=colMeans(col_sum_oneChr)
  idx0=which(can_bin_num<inter_bin_apper)
  col_sum_oneChr[idx0]=0


  pdf(out_heatmap,width = 10,height = 70)
  grid.arrange(grobs=ps,ncol=2)
  dev.off()
  if (sum(col_sum_oneChr) >0){
    HIRs_oneChr = call_HIRs_oneChr(col_sum_oneChr,chr,rgs)
    }else{
      HIRs_oneChr=c()
    }
  return(HIRs_oneChr)
}



call_HIRs_allChrs_inter <- function(mat_merged,chrs,rgs_all,quan_min_call_HIRs=0.9,continous_bins=3, out_dir,inter_bin_apper=10){
	heatmap_dir = paste(out_dir,"/heatmap_inter",sep="")
    dir.create(heatmap_dir)
    HIRs=data.frame()
    bdg=data.frame()
    bin_size=rgs_all[1,3]-rgs_all[1,2]
    for (chr in chrs){
    	out_heatmap=paste(heatmap_dir,"/",chr,'.pdf',sep="")
		HIRs_info=call_HIRs_oneChr_inter(mat_merged, chr, chrs,rgs_all, quan_min_call_HIRs=quan_min_call_HIRs,continous_bins=continous_bins, out_heatmap,inter_bin_apper)
        if (length(HIRs_info) != 0){
            h = HIRs_info[[1]]
            r = HIRs_info[[2]]
            HIRs=rbind(HIRs,as.data.frame(h))
            bdg=rbind(bdg,as.data.frame(r))
        }
    }
    HIRs=HIRs[1:3]
    out_HIR=paste(out_dir,"/","HIRs.inter.bed",sep="")
    out_bdg=paste(out_dir,"/","Regions.inter.bdg",sep="")
    write_regions(HIRs, out_HIR)
    write_regions(bdg, out_bdg)  
}

#######
HIRs_overlapped <- function(intra_file,inter_file,out_file){
  intra=read.table(intra_file)
  inter=read.table(inter_file)

  colnames(intra)=c("chr", "start", "end")
  colnames(inter)=c("chr", "start", "end")

  intra=as.data.table(intra)
  inter=as.data.table(inter)

  setkey(intra, chr, start, end)
  setkey(inter, chr, start, end)

  over <- foverlaps(intra, inter, nomatch = 0)
  # Extract exact regions
  over2 <- data.table(
      chr = over$chr,
      start = over[, ifelse(start > i.start, start, i.start)],
      end = over[, ifelse(end < i.end, end, i.end)])
  HIRs=as.data.frame(over2)
  write_regions(HIRs, out_file)
}


###########
mat_merged=read_matrix(opt$input)
rgs_all=read_binned_rgs(opt$regions)
dir.create(opt$out_dir)
chrs=sort(unique(as.vector(rgs_all[,1])))
print("Call HIRs intra")
call_HIRs_allChrs_intra(mat_merged,chrs,rgs_all,remove_diag_bins_ratio=opt$rm_diag_ratio, quan_min_call_HIRs=opt$quan_min_intra,continous_bins=opt$continous_bins, opt$out_dir)

print("Call HIRs inter")
call_HIRs_allChrs_inter(mat_merged,chrs,rgs_all,quan_min_call_HIRs=opt$quan_min_inter,continous_bins=opt$continous_bins, opt$out_dir,opt$inter_bin_apper)

intra_file=paste(opt$out_dir,"/","HIRs.intra.bed",sep="")
inter_file=paste(opt$out_dir,"/","HIRs.inter.bed",sep="")
out_file=paste(opt$out_dir,"/","HIRs.bed",sep="")

HIRs_overlapped(intra_file,inter_file,out_file)

print("Finished!")



