source('~/Development/phd/reconstruct_kinase/plot/exp-vs-pred-figure-panel.R')

PositionPlot <- function(file.in, dir.out){
  
  # Trim factor: x offset
  trim.factor = c('1'=-0.003, '2'=-0.001, '3'=0, '4'=0, '5'=0.003,
                  '6'=0.004, '7'=0.005, '8'=0.007, '9'=0.008, '10'=0.009,
                  '11'=0.010, '12'=0.011, '13'=0.013, '14'=0.0155, '15'=0.016)
  
  pos = -7:7
  pic = Pdf2Picture(pdf.file = file.in)
  
  kin.name = gsub('.pdf', '',basename(file.in))
  for(i in 1:1){
    path = sprintf('%s/%s_%s.svg', dir.out, kin.name, pos[i])
    svg(path, width = 819/96, height=313/96)
#     path = sprintf('%s/%s_%s.pdf', dir.out, kin.name, pos[i])
#     pdf(path, width = 819/96, height=297/96)
    
    d = 0.055
    f = 0.17 - d
    f = f-trim.factor[i]
    
    fx = f + (i*d)
    
    grid.roundrect(gp=gpar(fill="#FF4040", alpha=0.3, lwd=0, col=NA), x=fx, y=0.545,
                   r=unit(0.16, "snpc"), width=0.057, height=0.745)
    grid.picture(pic)
    dev.off()
  }
}

PlotBgRects <- function(dir.out){
  
  # Trim factor: x offset
  trim.factor = c('1'=-0.003, '2'=-0.001, '3'=0, '4'=0, '5'=0.003,
                  '6'=0.004, '7'=0.005, '8'=0.007, '9'=0.008, '10'=0.009,
                  '11'=0.010, '12'=0.011, '13'=0.013, '14'=0.0155, '15'=0.016)
  
  pos = -7:7
  for(i in 1:15){
    path = sprintf('%s/%s.png', dir.out, pos[i])
    png(path, width = 819, height=331, bg='transparent')
    
    d = 0.055
    f = 0.17 - d
    f = f-trim.factor[i]
    
    fx = f + (i*d)
    
    grid.roundrect(gp=gpar(fill="#FF4040", alpha=0.3, lwd=0, col=NA), x=fx, y=0.545,
                   r=unit(0.16, "snpc"), width=0.057, height=0.745)
    dev.off()
  }
  
  
}


PlotBgRects <- function(dir.out){
  
  #pic = Pdf2Picture(pdf.file = '~/Desktop/x.pdf')
  # Trim factor: x offset
  trim.factor = c('1'=-0.003, '2'=-0.001, '3'=0, '4'=0, '5'=0.003,
                  '6'=0.004, '7'=0.005, '8'=0.007, '9'=0.008, '10'=0.009,
                  '11'=0.010, '12'=0.011, '13'=0.013, '14'=0.005, '15'=0.005)
  
  pos = -7:7
  for(i in 1:15){
#     path = sprintf('%s/%s.png', dir.out, pos[i])
#     png(path, width = 819, height=331, bg='transparent')
    path = sprintf('%s/%s.pdf', dir.out, pos[i])
    pdf(path, width = 819/96, height=331/96, bg='transparent')
    
    
    d = 0.055
    f = 0.17 - d
    #f = f-trim.factor[i]
    
    fx = f + (i*d)
    if(i < 8) fx = fx - 0.01
    if(i >= 8) fx = fx - 0.015
    
    grid.roundrect(gp=gpar(fill="pink", alpha=1, lwd=1, col='black'), x=fx, y=0.085,
                   r=unit(0.16, "snpc"), width=0.050, height=0.12)
    #grid.picture(pic)
    dev.off()
    
    path.svg = sprintf('%s/%s.svg', dir.out, pos[i])
    system(sprintf('pdf2svg %s %s', path, path.svg))
    unlink(path)
  }
  
  
}

PlotBgRects('~/Development/mimp/inst/extdata/html/images/highlight/')
#PositionPlot('~/Desktop/x.pdf','~/Desktop/temp/')
#'~/Desktop/x.pdf'