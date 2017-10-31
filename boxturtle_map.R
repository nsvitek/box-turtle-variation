library(maps)
library(ggrepel)
library(ggmap)
library(scales) #needed to make reverse log-10 y axis, as solved by stackexchange Brian Diggs
reverselog_trans <- function(base = exp(1)) {
  trans <- function(x) -log(x, base)
  inv <- function(x) base^(-x)
  trans_new(paste0("reverselog-", format(base)), trans, inv, 
            log_breaks(base = base), 
            domain = c(1e-100, Inf))
}

# make map --------

localities_fossil<-fos_metadata %>% select(.,latitude,longitude,site) %>% unique
localities_ssp<-ssp_metadata %>% select(.,latitude,longitude,sex,ssp,subspecies) %>% unique
localities_sex<-dim_metadata %>% select(.,latitude,longitude,sex,ssp,id) %>% unique

usa<-map_data('state')
ylims<-range(localities_ssp$latitude)+c(-0.8,+0.8)
xlims<-range(localities_ssp$longitude)+c(-0.8,+0.8)

#map1 or map contains legends, map2 does not.
cairo_pdf("map2.pdf",width = 5.4, height = 4.4,family ="Arial",bg="transparent")
ggplot(data=localities_ssp,aes(x=longitude,y=latitude))+
  geom_polygon(data=usa,aes(x=long,y=lat,group=group),fill=NA,color="gray") +
  coord_cartesian(xlim=xlims,ylim=ylims)+
  theme_classic() +
  theme(axis.line=element_blank(),axis.text.x=element_blank(),
      axis.text.y=element_blank(),axis.ticks=element_blank(),
      axis.title.x=element_blank(),axis.title.y=element_blank(),
      legend.position=c(0.9,0.3),text=element_text(size=10),
      panel.background = element_rect(fill = "transparent",colour = NA), # or theme_blank()
      panel.grid.minor = element_blank(), 
      panel.grid.major = element_blank(),
      plot.background = element_rect(fill = "transparent",colour = NA))+
  geom_point(aes(color=ssp,shape=sex),size=2)+
  geom_point(data=localities_sex,aes(color=ssp,shape=sex),size=2)+
  geom_point(data=localities_fossil,color="black",size=2,shape=3)+
  # geom_label_repel(data=localities_fossil,aes(x=longitude,y=latitude,label=site),size=3,
  #                  box.padding=unit(0.45,"lines"),label.padding=unit(0.05,"lines"))+
  # scale_shape_manual(values=c(15,16),name="Sex",
  #                    labels=c("female","male"))  +
  # scale_color_manual(values= ssp_col,name="Subspecies",
  #                   labels=c("T. c. bauri","T. c. carolina","T. c. major","T. c. triunguis"))
  scale_color_manual(values= ssp_col,guide=FALSE) +
  scale_shape_manual(values=c(15,16),guide=FALSE)
dev.off()
embed_fonts("map2.pdf")



# make stratigraphic chart -----
localities_strat<-fos_metadata %>% select(.,latitude,longitude,site) %>% unique
localities_strat$site<-as.character(localities_strat$site)
localities_strat$min<-0
localities_strat$max<-0
localities_strat$group<-NA

localities_strat[which(localities_strat$site=="friesenhahn cave"),c(3:6)]<-c("A. Friesenhahn Cave",17.8,19.6,"3 neither")
localities_strat[which(localities_strat$site=="ingleside"),c(3:6)]<-c("B. Ingleside",120,130,"3 maybe neither")
localities_strat[which(localities_strat$site=="haile 8a"),c(3:6)]<-c("C. Haile VIIIA",9.5,160,"3 neither")
localities_strat[which(localities_strat$site=="devils den"),c(3:6)]<-c("D. Devil's Den",9.5,12.7,"2 maybe peninsular")
localities_strat[which(localities_strat$site=="grove's orange midden"),c(3:6)]<-c("E. Grove's Orange Midden",3.8,6.2,"2 peninsular")
localities_strat[which(localities_strat$site=="reddick 1b"),c(3:6)]<-c("F. Reddick IB",9.5,160,"3 neither")
localities_strat[which(localities_strat$site=="fort center"),c(3:6)]<-c("G. Fort Center",0.45,2.8,"2 peninsular")
localities_strat[which(localities_strat$site=="melbourne"),c(3:6)]<-c("H. Melbourne",9.5,12.7,"1 maybe mainland")
localities_strat[which(localities_strat$site=="camelot"),c(3:6)]<-c("J. Camelot",160,600,"1 not mainland")
localities_strat[which(localities_strat$site=="ardis"),c(3:6)]<-c("K. Ardis",18.5,22,"3 maybe neither")
localities_strat[which(localities_strat$site=="vero"),c(3:6)]<-c("I. Vero",9.5,12.7,"2 maybe peninsular")

localities_strat$min<-as.numeric(localities_strat$min)
localities_strat$max<-as.numeric(localities_strat$max)
localities_strat<-arrange(localities_strat,site)

cairo_pdf("strat.pdf",width = 5.4, height = 3,family ="Arial")
ggplot(data=localities_strat,aes(x=site,y=min,color=group))+
  geom_linerange(aes(ymin=min,ymax=max),size=10) + 
  scale_y_continuous(trans=reverselog_trans(10)) +
  xlab("") + 
  ylab("Years Before Present (ka)") +
  theme_minimal() +
  theme(text=element_text(size=10),
        axis.text.x=element_text(angle=-45,hjust=1,color="black")) +
  scale_x_discrete(position = "top") +
  scale_color_manual(values= fos_col,guide=FALSE)
dev.off()
embed_fonts("strat.pdf")
