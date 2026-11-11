# Base-R layout adapted to the user's SuSiE-mix simulation figures.
method_colours = c(fSuSiE="#2489FF",
  mvSuSiE_diagonal_fixed="#CC79A7",mvSuSiE_diagonal_estimated="#D81B60",
  mvSuSiE_canonical_fixed="#56B4E9",mvSuSiE_canonical_estimated="#009E73",
  mvSuSiE_mash_fixed="#E69F00",mvSuSiE_mash_estimated="#D55E00")
method_labels = c(fSuSiE="fSuSiE",
  mvSuSiE_diagonal_fixed="mvSuSiE diagonal: V fixed",mvSuSiE_diagonal_estimated="mvSuSiE diagonal: V updated",
  mvSuSiE_canonical_fixed="mvSuSiE canonical: V fixed",mvSuSiE_canonical_estimated="mvSuSiE canonical: V updated",
  mvSuSiE_mash_fixed="mvSuSiE mash: V fixed",mvSuSiE_mash_estimated="mvSuSiE mash: V updated")
method_pch = setNames(rep(16,length(method_colours)),names(method_colours))

save_figure = function(path,width,height,draw){
  grDevices::pdf(paste0(path,".pdf"),width=width,height=height,useDingbats=FALSE)
  tryCatch(draw(),finally=grDevices::dev.off())
  grDevices::png(paste0(path,".png"),width=width,height=height,units="in",res=180)
  tryCatch(draw(),finally=grDevices::dev.off())
}

panel_figure = function(cfg,title,xlab,ylab,xlim,ylim,xticks,yticks,
                        draw_panel,note="",reference=NULL,
                        legend_methods=names(method_colours)){
  scenarios = expand.grid(noise_type=cfg$noise_type,f_type=cfg$f_type,stringsAsFactors=FALSE)
  nr = nrow(scenarios)
  nc = length(cfg$pve)
  panels = matrix(seq_len(nr*nc),nrow=nr,byrow=TRUE)
  layout(cbind(panels,nr*nc+seq_len(nr)),widths=c(rep(3.5,nc),.95))
  par(oma=c(8,3.1,3.2,.2),mar=c(2,2.7,1.5,.45),
      mgp=c(1.3,.4,0),tcl=-.2,family="sans",cex=.9,xpd=FALSE)
  for (i in seq_len(nr)) for (j in seq_len(nc)){
    plot(NA,xlim=xlim,ylim=ylim,xaxs="i",yaxs="i",axes=FALSE,xlab="",ylab="")
    abline(v=xticks,h=yticks,col="#DEDEDE",lwd=.8)
    if (!is.null(reference)) abline(h=reference,lty=2,lwd=1.2,col="#333333")
    axis(1,at=xticks,col="#777777",col.axis="#333333")
    axis(2,at=yticks,labels=if (j==1) yticks else FALSE,las=1,col="#777777",col.axis="#333333")
    box(col="#777777",lwd=.6)
    if (i==1) mtext(paste0("PVE = ",100*cfg$pve[j],"%"),side=3,line=.4,font=2)
    draw_panel(scenarios$f_type[i],scenarios$noise_type[i],cfg$pve[j])
  }
  for (i in seq_len(nr)){
    par(mar=c(2,.05,1.5,0))
    plot.new(); plot.window(c(0,1),c(0,1))
    shape = if (scenarios$f_type[i]=="step") "Step effect" else "Smooth effect"
    noise = if (scenarios$noise_type[i]=="iid") "IID noise" else "GP noise"
    text(.02,.5,paste(shape,noise,sep="\n"),adj=0,cex=.9)
  }
  mtext(if (nzchar(note)) paste(title,"-",note) else title,side=3,outer=TRUE,line=1.5,font=2,cex=1.05)
  mtext(xlab,side=1,outer=TRUE,line=.8,cex=1)
  mtext(ylab,side=2,outer=TRUE,line=1.5,cex=1)
  par(fig=c(0,1,0,1),oma=c(0,0,0,0),mar=c(0,0,0,0),new=TRUE,xpd=NA)
  plot.new(); plot.window(c(0,1),c(0,1),xaxs="i",yaxs="i")
  legend(.5,.055,legend=method_labels[legend_methods],col=method_colours[legend_methods],
    pch=method_pch[legend_methods],lty=0,
    bty="n",ncol=min(2,length(legend_methods)),xjust=.5,yjust=.5,cex=.78)
}

plot_study = function(summary,curves,calibration,cfg,fdir,note=""){
  methods = names(method_colours)
  offsets = setNames(seq(-.18,.18,length.out=length(methods)),methods)
  width = 2.2+4.2*length(cfg$pve)
  height = 2.8+1.65*length(cfg$f_type)*length(cfg$noise_type)
  xticks = cfg$L
  for (metric in c("coverage","cs_power","cs_available","mean_cs_size")){
    if (metric=="coverage"){
      title = paste0("Observed ",100*cfg$coverage,"% credible-set coverage"); ylab = "Observed coverage (%)"
    } else if (metric=="cs_power"){
      title = "Causal SNP capture by credible sets"; ylab = "Power (%)"
    } else if (metric=="cs_available"){
      title = "Simulations reporting a credible set"; ylab = "Simulations with at least one CS (%)"
    } else {
      title = "Mean credible-set size"; ylab = "SNPs per credible set"
    }
    percent = metric!="mean_cs_size"
    column = if (percent) paste0(metric,"_mean") else metric
    factor = if (percent) 100 else 1
    finite = summary[[column]][is.finite(summary[[column]])]
    ylim = if (percent) c(-3,103) else c(0,max(c(1,finite))*1.12)
    yticks = if (percent) seq(0,100,25) else pretty(ylim,4)
    yticks = yticks[yticks>=ylim[1] & yticks<=ylim[2]]
    draw = function(){
      panel_figure(cfg,title,"Number of causal SNPs",ylab,
        range(xticks)+c(-.35,.35),ylim,xticks,yticks,note=note,
        reference=if (metric=="coverage") 100*cfg$coverage else NULL,
        draw_panel=function(shape,noise,pve){
          d = summary[summary$f_type==shape & summary$noise_type==noise & summary$pve==pve,,drop=FALSE]
          for (method in methods){
            s = d[d$method==method,,drop=FALSE]
            x = s$L+offsets[method]
            if (percent && isTRUE(cfg$show_intervals)){
              lower = factor*s[[paste0(metric,"_lower")]]
              upper = factor*s[[paste0(metric,"_upper")]]
              ok = is.finite(lower) & is.finite(upper) & upper>lower
              arrows(x[ok],lower[ok],x[ok],upper[ok],angle=90,code=3,length=.035,
                col=grDevices::adjustcolor(method_colours[method],alpha.f=.6),lwd=1)
            }
            points(x,factor*s[[column]],pch=method_pch[method],col=method_colours[method],cex=.85)
          }
        })
    }
    save_figure(file.path(fdir,paste0(metric,"_by_L")),width,height,draw)
  }
  if (nrow(curves)) for (L in cfg$L){
    draw = function(){
      panel_figure(cfg,paste0("PIP power and FDR (",L," causal SNP",if (L>1) "s" else "",")"),
        "SNP false discovery rate (%)","SNP power (%)",c(-1,100*cfg$fdr_plot_max),c(-3,103),
        pretty(c(0,100*cfg$fdr_plot_max),5),seq(0,100,25),note=note,
        draw_panel=function(shape,noise,pve){
          d = curves[curves$f_type==shape & curves$noise_type==noise & curves$pve==pve & curves$L==L,,drop=FALSE]
          for (method in methods){
            s = d[d$method==method,,drop=FALSE]
            s = s[order(s$threshold,decreasing=TRUE),,drop=FALSE]
            points(100*s$fdr,100*s$power,col=method_colours[method],pch=16,cex=.5)
            mark = abs(s$threshold-cfg$pip_threshold)<1e-10
            points(100*s$fdr[mark],100*s$power[mark],col=method_colours[method],pch=16,cex=1)
          }
        })
    }
    save_figure(file.path(fdir,paste0("PIP_power_FDR_L",L)),width,height,draw)
  }
  if (nrow(calibration)) for (method in methods){
    draw = function(){
      panel_figure(cfg,paste0("PIP calibration: ",method_labels[method]),
        "Mean PIP in bin (%)","Observed causal SNP frequency (%)",c(-3,103),c(-3,103),
        seq(0,100,25),seq(0,100,25),note=note,legend_methods=method,
        draw_panel=function(shape,noise,pve){
          d = calibration[calibration$method==method & calibration$f_type==shape &
            calibration$noise_type==noise & calibration$pve==pve,,drop=FALSE]
          abline(0,1,lty=2,col="#333333")
          points(100*d$mean_pip,100*d$observed,pch=method_pch[method],col=method_colours[method],cex=.85)
        })
    }
    save_figure(file.path(fdir,paste0("PIP_calibration_",method)),width,height,draw)
  }
  invisible(fdir)
}
