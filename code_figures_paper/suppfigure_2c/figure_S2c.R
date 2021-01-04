library(ggplot2)
library(ggpmisc)
library(scales)

source("pcawg.colour.palette.R")

scale_override <- function(which, scale) {
  if(!is.numeric(which) || (length(which) != 1) || (which %% 1 != 0)) {
    stop("which must be an integer of length 1")
  }
  
  if(is.null(scale$aesthetics) || !any(c("x", "y") %in% scale$aesthetics)) {
    stop("scale must be an x or y position scale")
  }
  
  structure(list(which = which, scale = scale), class = "scale_override")
}

CustomFacetWrap <- ggproto(
  "CustomFacetWrap", FacetWrap,
  init_scales = function(self, layout, x_scale = NULL, y_scale = NULL, params) {
    scales <- ggproto_parent(FacetWrap, self)$init_scales(layout, x_scale, y_scale, params)
    
    if(is.null(params$scale_overrides)) return(scales)
    
    max_scale_x <- length(scales$x)
    max_scale_y <- length(scales$y)
    
    for(scale_override in params$scale_overrides) {
      which <- scale_override$which
      scale <- scale_override$scale
      
      if("x" %in% scale$aesthetics) {
        if(!is.null(scales$x)) {
          if(which < 0 || which > max_scale_x) stop("Invalid index of x scale: ", which)
          scales$x[[which]] <- scale$clone()
        }
      } else if("y" %in% scale$aesthetics) {
        if(!is.null(scales$y)) {
          if(which < 0 || which > max_scale_y) stop("Invalid index of y scale: ", which)
          scales$y[[which]] <- scale$clone()
        }
      } else {
        stop("Invalid scale")
      }
    }
      scales
  }
)

facet_wrap_custom <- function(..., scale_overrides = NULL) {
  facet_super <- facet_wrap(...)
    if(inherits(scale_overrides, "scale_override")) {
    scale_overrides <- list(scale_overrides)
  } else if(!is.list(scale_overrides) || 
            !all(vapply(scale_overrides, inherits, "scale_override", FUN.VALUE = logical(1)))) {
    stop("scale_overrides must be a scale_override object or a list of scale_override objects")
  }
  
  facet_super$params$scale_overrides <- scale_overrides
  
  ggproto(NULL, CustomFacetWrap,
          shrink = facet_super$shrink,
          params = facet_super$params
  )
}

df<-read.table("input.figureS3c.txt", sep="\t", stringsAsFactors = F, header=T)

organ.colours<-pcawg.colour.palette(x=c(unique(df$colour)," ","All"),scheme="tumour.subtype")
names(organ.colours)<-c(unique(df$colour)," ", "All")

df$num_muts<-(df$num_muts/1000)


p1<-ggplot(df, aes(x=num_muts, y=(1-frac_clonal_postWCC)))+
  geom_point(aes(fill=colour, shape = pch),colour="black", stroke=0.05)+
  scale_y_continuous(name="Fraction subclonal SNVs", limit=c(0,1),oob=squish)+
  scale_x_log10(expression(paste("Mutation burden (10"^"3",")")),breaks = log_breaks())+
  scale_fill_manual(values=organ.colours,name="")+
  scale_shape_identity()+
  geom_smooth(method='lm')+
  stat_poly_eq(formula = y~x, 
               aes(label = paste(..rr.label.., sep = "~~~")), 
               parse = TRUE) +
  guides(fill = FALSE, shape = FALSE)+
  theme_bw()+
  theme(panel.grid.minor = element_blank())+
  facet_wrap_custom(~organ, scales = "free_x", nrow = 4, scale_overrides = list(
    scale_override(4, scale_x_continuous(breaks = c(log10(2),log10(5),log10(10)), labels=c(2,5,10))),
    scale_override(5, scale_x_continuous(breaks = c(log10(3),log10(10),log10(30),log10(100)), labels=c(3,10,30,100))),
    scale_override(7, scale_x_continuous(breaks = c(log10(0.5),log10(1), log10(2.0), log10(5.0)), labels=c(0.5,1,2,5))),
    scale_override(8, scale_x_continuous(breaks = c(log10(0.1),log10(0.2),log10(0.5)), labels=c(0.1,0.2,0.5))),
    scale_override(14, scale_x_continuous(breaks = c(log10(2),log10(5),log10(10)), labels=c(2,5,10))),
    scale_override(19, scale_x_continuous(breaks = c(log10(1), log10(2.0), log10(5.0)), labels=c(1,2,5))),
    scale_override(22, scale_x_continuous(breaks = c(log10(0.3),log10(1),log10(3),log10(10), log10(30)), labels=c(0.3,1,3,10,30))),
    scale_override(24, scale_x_continuous(breaks = c(log10(2),log10(5),log10(10), log10(20)), labels=c(2,5,10,20))),
    scale_override(26, scale_x_continuous(breaks = c(log10(2),log10(5),log10(10)), labels=c(2,5,10))),
    scale_override(27, scale_x_continuous(breaks = c(log10(2),log10(5),log10(10)), labels=c(2,5,10))),
    scale_override(29, scale_x_continuous(breaks = c(log10(0.5),log10(1), log10(2.0)), labels=c(0.5,1,2)))
  ))

ggsave(plot=p1,height=7,width=14, filename="SuppFig_2c.pdf", useDingbats=FALSE)
