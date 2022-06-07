 # easymanh: A flexible script for plotting Manhattan plots with different SNP scores

## Main parameters are:

* **eff:** *Mandatory*. The vector of SNPs numeric scores e.g. pvalues, % of variance, Random forest importance scores. Default: -log10(pvalue);  
* **chr:** *Mandatory*. The vector of numeric values for the *chromosomes*;  
* **pos:** *Mandatory*. The vector of numeric marker positions;  
* **col.op:** *Optional*. The desirable colors for the plot. Ex: col.op = c("blue","gray"), col.op = c(4,5,2,6). Default: col.op = c("black","gray");   
* **sign:** *Optional*. A constant for indicating the significance threshold;  
* **sig_col:** *Optional*. The color for the significance line;  
* **y_lab:** *Optional*. The label for y_axis. default: -log10(pvalues);  
* **x_lab:** *Optional*. The label for x_axis. default: Chromossomes;  
* **ylim:** *Optional*. A numeric vector indicating the range for y_axis values. Default: ylim = c(min(eff), max(eff));   
* **alpha:** *Optional*. A constant between 0 and 255 for indicating the transparency level of the plot points;   
* **size:** *Optional*. The size of the plot points. default: 1.2;  
* **custom_labels:** *Optional*. A vector of labels for the Chromosomes. Ex: custom_labels = c("1", "2", "3", "4", "X", "Y", "MT")

## Example:
#Not run
```R
easymanh(eff = eff_walk1$V3, pos = eff_walk1$V4, 
         chr = chr, col.op = c(rgb(0.4,0.4,0.6,1), "gray"),
         xlab = "Windows", alpha = 175, sign = 6, sig_col = "red")
#add a tittle if necessary
title("Trait A",cex = 1.2)
```
## This plot would look like:
