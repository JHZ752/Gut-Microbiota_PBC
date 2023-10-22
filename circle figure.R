
library(tidyverse)
library(ggplot2)
library(circlize)
library(ComplexHeatmap)

# 数据读取：
data <- read.csv("data.csv",header = T)
data[1:5, 1:5]


# 绘图:
library(circlize)

################### 外圈热图 ========================
# 设置颜色模式:
col_fun = list(col_1 = colorRamp2(c(0, 0.5, 1), 
                                  c("#fdbf6f", "white","#6779b2")),
               col_2 = colorRamp2(c(0,0.5, 1),
                                  c("#fdbf6f", "white","#6779b2")),
               col_3 = colorRamp2(c(0, 0.5, 1),
                                  c("#fdbf6f", "white","#6779b2")),
               col_4 = colorRamp2(c(0,0.5, 1),
                                  c("#fdbf6f", "white","#6779b2")),
               col_5 = colorRamp2(c(0, 0.5, 1),
                                  c("#fdbf6f", "white","#6779b2"))
)

circos.clear()
if (T) {
  pdf("plot.pdf", height = 7, width = 7)
  ################### 外圈热图 ========================
  circos.par("track.height" = 0.1,   # 每行的高度
             "start.degree" = 0,  # 环形开始的角度
             "gap.degree" = c(5, 5, 5, 5, 30), # sector之间的gap角度
             "track.margin" = c(0.01, 0.01))  # 每行之间的距离
  
  data_1 <- data[,-4]
  
  # 循环绘制热图：
  for (i in 1:6) {
    data_tmp <- as.matrix(as.numeric(data_1[, i+2]))
    if (i == 1) {
      rownames(data_tmp) <- data$ID
    }
    colnames(data_tmp) <- paste0(colnames(data_tmp), "-p")
    
    if (i < 6) {
      circos.heatmap(data_tmp, 
                     # 设置sector：
                     split = data$Level,
                     col = col_fun[[i]],
                     rownames.side = "outside", 
                     rownames.cex = 0.5,  # 行名字体大小
                     cluster = F,
                     cell.border = NA,  # 单元格边框
                     track.height = 0.05)
    } else {
      ################### 内圈散点图 ========================
      data_tmp <- as.matrix(as.numeric(data[, 4]))
      colnames(data_tmp) <- "IVW-OR"
      
      circos.track(ylim = range(data_tmp[,1]), 
                   panel.fun = function(x, y) {
                     # sector的名称：
                     circos.text(CELL_META$xcenter, 
                                 CELL_META$cell.ylim[2] + mm_y(1.5),
                                 CELL_META$sector.index,
                                 cex = 0.5
                     )
                     # y轴刻度：
                     circos.yaxis(labels.cex = 0.4, 
                                  at = seq(0.5, 1.5, 0.5), #自己根据OR值调整合适的刻度
                                  sector.index = "phylum")#位于data的第一列的第一组，可以用View(data)查看
                     y = data_tmp[,1][CELL_META$subset]
                     y = y[CELL_META$row_order]
                     # 虚线：
                     circos.lines(CELL_META$cell.xlim, c(1, 1), lty = "dashed", col = "black")
                     # 散点：
                     circos.points(seq_along(y) - 0.5, y, 
                                   pch = 16, cex = 0.5,
                                   col = "#ff6666")#散点的颜色
                     # 行名：
                     if(CELL_META$sector.numeric.index == 5) { # the last sector
                       cn = rev(c("IVW-p", "MR Egger-p", "WM-p", "MLE-p","WMODE-p", "IVW-OR"))
                       n = length(cn)
                       circos.text(rep(CELL_META$cell.xlim[2], n) + convert_x(8, "mm"), 
                                   c(1, seq(3.5, 9.5, 1.5)), cn, 
                                   cex = 0.7, adj = 0.5, facing = "inside")
                     }
                   }, 
                   track.margin = c(0.02, 0.02),
                   cell.padding = c(0.02, 0, 0.02, 0))
    }
  }
  circos.clear()
  dev.off()
}


getwd()


# For legends
# Create a blank plot
plot(0, 0, type = "n", xlab = "", ylab = "", xlim = c(0, 1), ylim = c(0, 1), main = "Blank Plot")

# Define the color ramps and titles for each plot
color_ramps <- list(
  list(c(0, 0.5, 1), c("#fdbf6f", "white","#6779b2"), "IVW-pvalue"),
  list(c(0, 0.5, 1), c("#fdbf6f", "white","#6779b2"), "MR Egger-pvalue"),
  list(c(0, 0.5, 1), c("#fdbf6f", "white","#6779b2"), "Weighted median-pvalue"),
  list(c(0, 0.5, 1), c("#fdbf6f", "white","#6779b2"), "MLE-pvalue"),
  list(c(0, 0.5, 1), c("#fdbf6f", "white","#6779b2"), "Weighted mode-pvalue")
)

# Loop through the color ramps and generate the plots
for (i in 1:length(color_ramps)) {
  # Generate the PDF file for each plot
  pdf(paste0("plot", i, ".pdf"), height = 7, width = 7)
  
  # Retrieve the color ramp and title for the current plot
  col_fun <- colorRamp2(color_ramps[[i]][[1]], color_ramps[[i]][[2]])
  title <- color_ramps[[i]][[3]]
  
  # Create the legend
  lgd <- Legend(col_fun = col_fun, title = title)
  
  # Draw the legend
  draw(lgd, test = "only col_fun")
  
  # Close the PDF file
  dev.off()
}

