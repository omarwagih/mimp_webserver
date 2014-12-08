require(RWebLogo)
source('~/Development/phd/reconstruct_kinase/init.R')
ks = readRDS('~/Desktop/MSc project/rewiring/data/new/ksr_newman.rds')
setwd('~/Desktop/')

logo_folder = '~/Desktop/logos_newman_rect/'

for(kin in names(ks)){
  pdf.out = sprintf('%s/%s.pdf', logo_folder, kin)
  png.out = sprintf('%s/%s.png', logo_folder, kin)
  weblogo(ks[[kin]], open = F, file.out = pdf.out, yaxis = 4, errorbars = F, annotate = -7:7, color.scheme='chemistry3', title=kin, scale.width = F)
  
  
  system(sprintf('convert -density 300 -quality 100 %s %s', pdf.out, png.out))
  
  unlink(pdf.out)
  
  shannon = attr(PWM(ks[[kin]]), 'shan')
  
  xstart = 160
  xwidth = 67
  
  abs_pos = 1:15
  rel_pos = -7:7
  
  for(i in abs_pos){
    x0 = xstart + (xwidth * (i-1) )
    x1 = x0 + xwidth
    y0 = 50
    y1 = 412
    
    bit = (y1 - y0)/4.2
    z = shannon[i]
    if(z > 4) z = 4
    z = bit * z
    y0. = ceiling(  y1 - z )
    
    png.dist.out = sprintf('%s/%s_%s.png', logo_folder, kin, rel_pos[i])
    system(sprintf('convert -fill transparent -stroke red -strokewidth 3 -draw "rectangle %s,%s %s,%s" %s %s', x0, y0., x1, y1, 
                   png.out, png.dist.out))
  }
  
}
