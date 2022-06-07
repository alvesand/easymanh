#eff: *Mandatory. The vector of SNPs numeric scores e.g. pvalues, % of variance, Random forest importance scores. Default: -log10(pvalue);
#chr: *Mandatory. The vector of numeric values for the chromossomes;
#pos: *Mandatory, The vector of numeric marker positions;
#col.op: Optional. The desirable colors for the plot. Ex: col.op = c("blue","gray"), col.op = c(4,5,2,6). Default: col.op = c("black","gray");   
#sign: Optional. A constant for indicating the significance threshold;
#sig_col: Optional. The color for the significance line;
#y_lab: Optional. The label for y_axis. default: -log10(pvalues);
#x_lab: Optional. The label for x_axis. default: Chromossomes;
#ylim: Optional. A numeric vector indicating the range for y_axis values. Default: ylim = c(min(eff), max(eff)); 
#alpha: Optional. A constant between 0 and 255 for indicating the transparency level of the plot points; 
#size: Optional. The size of the plot points. default: 1.2;
#custom_labels: Optional. A vector of characters indicating how chomossomes must be labeled. Ex: custom_labels = c(1,2,3,"X", "Y");

easymanh = function(eff, chr, pos,  col.op = col.op, 
                    sign = sign, sig_col = sig_col, 
                    ylab=ylab,xlab=xlab, ylim=ylim, alpha=alpha,size=size,custom_labels){
  
  if (missing(eff)==T) stop('Please inform the vector for the SNP effects or p-values ');
  if (missing(pos)==T) stop('Please inform the vector of SNP positions');
  if (missing(chr)==T) stop('Please inform the vector of Chromossomes IDs');
  
  if (is.numeric(chr)==F) stop('chr vector must be numeric')
  if (is.numeric(pos)==F) stop('pos vector must be numeric')
  if (is.numeric(eff)==F) stop('eff vector must be numeric')
  
  if (missing(ylim)==F) eff_max = ylim[2];
  if (missing (ylim)==F) eff_min = ylim[1];
  
  if (missing(ylim)==T) eff_max = max(eff);
  if (missing (ylim)==T) eff_min = min(eff);
  
  
  gwas_dat=data.frame(eff=eff, chr=chr,pos=pos)
  ind = order(gwas_dat$chr,gwas_dat$pos)
  gwas_dat=gwas_dat[ind,]
  
  if (missing (ylab) ==T) ylab = '-log10(pvalue)';
  if (is.character (ylab) ==F) stop('Please insert a valid label for axis y');
  if (missing (xlab) ==T) xlab = 'Chromossomes';
  if (missing(alpha)==T) alpha = 255;
  if (missing(size)==T) size = 1.2;
  if (any(gwas_dat$pos>1e6)) gwas_dat$pos<-gwas_dat$pos/1e6;
  if(missing(col.op)==F){
    k = length(col.op)  
    if (k == 1)  stop('Please set at least two colors for the plot!')}
  
  if (missing(col.op)==T) col.op = 'standard';
  
  chr = unique(gwas_dat$chr)
  
  if (missing(custom_labels)==F){
    if (length(custom_labels)!=length(chr)) stop('wrong size for custom labels')
  }
  gwas_dat$tmp = NA
  
  col.op=col.op
  
  makeTransparent<-function(someColor, alpha=100)
  {
    newColor<-col2rgb(someColor)
    apply(newColor, 2, function(curcoldata){rgb(red=curcoldata[1], green=curcoldata[2],
                                                blue=curcoldata[3],alpha=alpha, maxColorValue=255)})
  }
  
  pos1 = tail(gwas_dat$pos[which(gwas_dat$chr==chr[1])],n=1)
  pos2 = head(gwas_dat$pos[which(gwas_dat$chr==chr[2])],n=1)
  
  if (pos1>pos2){
    gwas_dat$tmp[which(gwas_dat$chr==chr[1])] = gwas_dat$pos[which(gwas_dat$chr==chr[1])]
    for (i in 2:length(chr)){
      end = range(gwas_dat$tmp[which(gwas_dat$chr==chr[i-1])])[2]
      gwas_dat$tmp[which(gwas_dat$chr==chr[i])] = gwas_dat$pos[which(gwas_dat$chr==chr[i])]+end
    }
  }else{
    gwas_dat$tmp = gwas_dat$pos
  }
  min_pos = min(gwas_dat$tmp)
  max_pos = max(gwas_dat$tmp)
  
  
  if(length (col.op) ==1){
    
    col.std = matrix(NA, length(chr),4)
    
    if (length(chr)%%2==0){
      col.std[,1] = rep(c(0,0.7),length(chr)/2)
      col.std[,2] = rep(c(0,0.7),length(chr)/2)
      col.std[,3] = rep(c(0,0.7),length(chr)/2)
      col.std[,4] = rep(0.7,length(chr))
    }else{
      col.std[,1] = rep(c(0,0.7),length(chr)%/%2+1)[-c(2*(length(chr)%/%2+1))]
      col.std[,2] = rep(c(0,0.7),length(chr)%/%2+1)[-c(2*(length(chr)%/%2+1))]
      col.std[,3] = rep(c(0,0.7),length(chr)%/%2+1)[-c(2*(length(chr)%/%2+1))]
      col.std[,4] = rep(0.7,length(chr))
    }
    
    plot(eff~tmp,data=gwas_dat[which(gwas_dat$chr=='1'),], axes= F,ylab = '', xlab = '', ylim = c(eff_min, eff_max), xlim = c(min_pos,max_pos), pch = 20, cex = size,col=rgb(col.std[1,1],col.std[1,2],col.std[1,3],col.std[2,4]))
    for (i in 2:length(chr)){
      points(eff~tmp,data=gwas_dat[which(gwas_dat$chr==chr[i]),], pch = 20, cex = size,col=rgb(col.std[i,1],col.std[i,2],col.std[i,3],col.std[i,4]))
    }
  }else{
    
    if (length(chr)%%2==0){
      col.op1 = rep(col.op,length(chr)/2)
    }else{
      col.op1 = rep(col.op,length(chr)%/%2+1)[-c(2*(length(chr)%/%2+1))]
    }
    
    plot(eff~tmp,data=gwas_dat[which(gwas_dat$chr=='1'),], axes= F,ylab = '', xlab = '', ylim = c(eff_min, eff_max), xlim = c(min_pos,max_pos), pch = 20, cex = size,col=makeTransparent(col.op1[1],alpha=alpha))
    for (i in 2:length(chr)){
      points(eff~tmp,data=gwas_dat[which(gwas_dat$chr==chr[i]),], pch = 20, cex = size,col=makeTransparent(col.op1[i],alpha=alpha))
    }
  }
  
  
  pos_id = NULL
  for (i in 1:length(chr)){
    pos_id[i] = median(gwas_dat$tmp[which(gwas_dat$chr==chr[i])])
  }
  
  if (missing(custom_labels)==F){  
    axis(1, at = pos_id,labels = custom_labels, tick = 0, line = -1)
  }else{
    axis(1, at = pos_id,labels = chr, tick = 0, line = -1)
  }
  axis(2, tick =0.1)
  
  title(ylab = ylab, line = 2.2, cex.lab = 1)
  title(xlab = xlab, line = 1.5, cex.lab = 1)
  
  if (missing(sign)==F){
    if (missing(sig_col)==T) sig_col = 'red'
    abline(a=sign, b = 0, col=sig_col,lty=2)
  }
  
}