library(readxl)
library(tidyverse)
library(janitor)
library(lubridate)
library(survival)
library(survminer)
library(tools)
library(cowplot)
library(car)

dates<-bind_cols(seq(1,9,1),c('07/22/21','07/29/21','08/05/21','08/12/21','08/26/21','09/09/21','09/23/21','10/07/21','11/04/21'))%>%
     rename(timepoint=1,date=2)

############################## SURVIVORSHIP ######################################## #####
meta<-read_xlsx("./data/Urchin Treatment Meta Data.xlsx")%>%rename(rack=1,treatment=2)

rack606<-read_xlsx("./data/UrchinSurvivorship.xlsx",sheet=1)%>%clean_names()%>%
     select(contains('x'),contains('ind'),contains('agg'))%>%
     drop_na()%>%
     rename(plug=1)%>%select(-contains('x'))%>%
     filter(plug!="TOTALS:")%>%
     gather(type,count,-plug)%>%
     separate(type,into=c("type","garbo"),sep="_")%>%select(-garbo)%>%
     arrange(plug)%>%
     group_by(plug,type)%>%
     mutate(timepoint = 1:n())%>%
     mutate(rack='rack_606')

rack62<-read_xlsx("./data/UrchinSurvivorship.xlsx",sheet=2)%>%clean_names()%>%
     select(contains('x'),contains('ind'),contains('agg'))%>%
     drop_na()%>%
     rename(plug=1)%>%select(-contains('x'))%>%
     filter(plug!="TOTALS:")%>%
     gather(type,count,-plug)%>%
     separate(type,into=c("type","garbo"),sep="_")%>%select(-garbo)%>%
     arrange(plug)%>%
     group_by(plug,type)%>%
     mutate(timepoint = 1:n())%>%
     mutate(rack='rack_62')

rack682<-read_xlsx("./data/UrchinSurvivorship.xlsx",sheet=3)%>%clean_names()%>%
     select(contains('x'),contains('ind'),contains('agg'))%>%
     drop_na()%>%
     rename(plug=1)%>%select(-contains('x'))%>%
     filter(plug!="TOTALS:")%>%
     gather(type,count,-plug)%>%
     separate(type,into=c("type","garbo"),sep="_")%>%select(-garbo)%>%
     arrange(plug)%>%
     group_by(plug,type)%>%
     mutate(timepoint = 1:n())%>%
     mutate(rack='rack_682')

rack165<-read_xlsx("./data/UrchinSurvivorship.xlsx",sheet=4)%>%clean_names()%>%
     select(contains('x'),contains('ind'),contains('agg'))%>%
     drop_na()%>%
     rename(plug=1)%>%select(-contains('x'))%>%
     filter(plug!="TOTALS:")%>%
     gather(type,count,-plug)%>%
     separate(type,into=c("type","garbo"),sep="_")%>%select(-garbo)%>%
     arrange(plug)%>%
     group_by(plug,type)%>%
     mutate(timepoint = 1:n())%>%
     mutate(rack='rack_165')

rack22<-read_xlsx("./data/UrchinSurvivorship.xlsx",sheet=5)%>%clean_names()%>%
     select(contains('x'),contains('ind'),contains('agg'))%>%
     drop_na()%>%
     rename(plug=1)%>%select(-contains('x'))%>%
     filter(plug!="TOTALS:")%>%
     gather(type,count,-plug)%>%
     separate(type,into=c("type","garbo"),sep="_")%>%select(-garbo)%>%
     arrange(plug)%>%
     group_by(plug,type)%>%
     mutate(timepoint = 1:n())%>%
     mutate(rack='rack_22')

rack478<-read_xlsx("./data/UrchinSurvivorship.xlsx",sheet=6)%>%clean_names()%>%
     select(contains('x'),contains('ind'),contains('agg'))%>%
     drop_na()%>%
     rename(plug=1)%>%select(-contains('x'))%>%
     filter(plug!="TOTALS:")%>%
     gather(type,count,-plug)%>%
     separate(type,into=c("type","garbo"),sep="_")%>%select(-garbo)%>%
     arrange(plug)%>%
     group_by(plug,type)%>%
     mutate(timepoint = 1:n())%>%
     mutate(rack='rack_478')

all_data<-bind_rows(rack165,rack478,rack606,rack62,rack682,rack22)%>%
     mutate(count=as.numeric(count),rack=as.factor(rack),type=as.factor(type))%>%
     left_join(.,dates,by="timepoint")%>%
     mutate(date=mdy(date))%>%
    left_join(.,meta,by="rack")%>%
    arrange(plug,timepoint)%>%
    mutate(start=case_when(timepoint==1~count))%>%
    fill(start,.direction = "down")

final_surv<-all_data%>%
    filter(timepoint==9)%>%
    mutate(surv=count/start)%>%
    mutate(surv=case_when(surv>1~1,TRUE~as.numeric(surv)))%>%
    select(type,plug,treatment,surv)%>%
    mutate(type=case_when(type=="aggregate"~"Aggregate",
                        type=="individual"~"Individual"))

leveneTest(surv~treatment,data=final_surv%>%filter(type=="Individual"))
wilcox.test(surv~treatment,data=final_surv%>%filter(type=="Individual"))
leveneTest(surv~treatment,data=final_surv%>%filter(type=="Aggregate"))
t.test(surv~treatment,data=final_surv%>%filter(type=="Aggregate"))

final_surv%>%group_by(type,treatment)%>%summarise(mean=mean(surv,na.rm=TRUE),sd=sd(surv,na.rm=TRUE))
all_data%>%filter(timepoint==1)%>%group_by(type)%>%summarise(sum=sum(count))
nrow(all_data%>%filter(timepoint==1)%>%group_by(plug)%>%distinct())

endpoint<-ggplot(final_surv)+
  geom_boxplot(aes(treatment,surv,fill=treatment))+
  geom_point(aes(treatment,surv),size=0.25)+
  facet_wrap(~type)+
  theme_classic(base_size=8)+
  xlab("Treatment")+ylab("Endpoint Survivorship")+
  scale_fill_manual(values=c("gray","#6F3C4B"),name="Treatment")+
  theme(legend.position="none")+
  scale_y_continuous(breaks=seq(0,1,0.2));endpoint

my_tag<-c("p<0.001","p=0.006")
tag_facet2 <- function(p, open = "(", close = ")", tag_pool = letters, x = -Inf, y = Inf, 
                       hjust = -0.5, vjust = 1.5, fontface = 2, family = "", ...) {
  
  gb <- ggplot_build(p)
  lay <- gb$layout$layout
  tags <- cbind(lay, label = paste0(open, tag_pool[lay$PANEL], close), x = x, y = y)
  p + geom_text(data = tags, aes_string(x = "x", y = "y", label = "label"), ..., hjust = hjust, 
                vjust = vjust, fontface = fontface, family = family, inherit.aes = FALSE)
}

endpoint<-tag_facet2(endpoint, 
              x = 1, y = 1.1, 
              vjust = 0, hjust = 0,
              open = "", close = "",
              fontface = 'italic',
              size = 2,
              tag_pool = my_tag);endpoint


survdata<-all_data%>%
    mutate(start=case_when(timepoint==1~count))%>%
    fill(start,.direction = "down")%>%
  mutate(dead=start-count)%>%
  mutate(dead=case_when(dead<0~0,TRUE~as.numeric(dead)))%>%
  rename(alive=count)%>%
  select(plug,type,rack,treatment,timepoint,date,alive,dead)%>%
  gather(condition,count,-plug,-type,-rack,-timepoint,-date,-treatment)%>%
  mutate(condition=case_when(condition=="alive"~0,TRUE~1))%>%
  uncount(count)%>%
  mutate(start=ymd('2021-07-22'))%>%
  mutate(duration=interval(start,date)/days(1))

fit_ind=survfit(Surv(duration,condition) ~treatment,data=survdata%>%filter(type=="individual"),type="kaplan-meier")
surv_pvalue(fit_ind,method="n") #early detection
surv_pvalue(fit_ind)

ind<-(ggsurvplot(fit_ind,size=0.5,conf.int=TRUE)$plot)+
  theme_classic(base_size=8)+
  scale_color_manual(values=c("gray","#6F3C4B"),labels=c("Control","Urchin"),name="Treatment")+
  scale_fill_manual(values=c("gray","#6F3C4B"),guide="none")+
  theme(legend.position=c(0.3,0.4),
        legend.key.size=unit(0.1,"cm"))+
  xlab("Days")+
  annotate("text",x=0,y=0.05,label="Individual\n~treatment p<0.001",size=2,fontface="italic",hjust=0);ind

fit_agg=survfit(Surv(duration,condition) ~treatment,data=survdata%>%filter(type=="aggregate"))
surv_pvalue(fit_agg,method="n") #early detection
surv_pvalue(fit_agg)

agg<-(ggsurvplot(fit_agg,size=0.5,conf.int=TRUE)$plot)+
  theme_classic(base_size=8)+
  scale_color_manual(values=c("gray","#6F3C4B"),guide="none")+
  scale_fill_manual(values=c("gray","#6F3C4B"),guide="none")+
  theme(legend.position=c(0.2,0.5),
        legend.key.size=unit(0.3,"cm"))+
  xlab("Days")+
  annotate("text",x=0,y=0.05,label="Aggregate\n~treatment p<0.001",size=2,fontface="italic",hjust=0);agg

quartz(w=5.2,h=1.7)
plot_grid(ind,agg,endpoint,nrow=1,rel_widths=c(1,1,1.3),labels=c("A","B","C"),label_size=8)

############################## GROWTH ############################################## #####
meta<-read_xlsx("./data/Urchin Treatment Meta Data.xlsx")%>%rename(rack=1,treatment=2)%>%
  separate(rack,into=c("a","b"))%>%unite(rack,a,b,sep="")

list<-as.vector(list.files(path="./data/",pattern="rack*"))
setwd("./data/")
for (i in list){
     assign(paste0(file_path_sans_ext(i,compression=TRUE)),read_csv(paste0("./",i))%>%
              mutate(file=paste0(file_path_sans_ext(i,compression=TRUE)))%>%
              separate(file,into=c("rack","date","trash"))%>%
              clean_names()%>%select(-trash,-chunk,-shape_type)%>%
              rename(!!paste0(.$date[1]) := "area_mm_2")%>%
              select(-date)%>%clean_names())
}
setwd("../")
out<-bind_rows(purrr::reduce(mget(ls(pattern = "rack165+")), full_join, by = c("layer","shape","rack")),
          purrr::reduce(mget(ls(pattern = "rack22+")), full_join, by = c("layer","shape","rack")),
          purrr::reduce(mget(ls(pattern = "rack478+")), full_join, by = c("layer","shape","rack")),
          purrr::reduce(mget(ls(pattern = "rack606+")), full_join, by = c("layer","shape","rack")),
          purrr::reduce(mget(ls(pattern = "rack62+")), full_join, by = c("layer","shape","rack")),
          purrr::reduce(mget(ls(pattern = "rack682+")), full_join, by = c("layer","shape","rack")))%>%
  left_join(.,meta,by="rack")%>%
  mutate(type=case_when(shape<=2~"individual",TRUE~"aggregate"))%>%
  select(rack,type,layer,treatment,everything(),-shape)%>%
  rename(plug=3,start=5,end=6)%>%
  mutate(growth=(end-start)/start)%>%
  mutate(type=case_when(type=="aggregate"~"Aggregate",
                        type=="individual"~"Individual"))

out%>%group_by(type)%>%tally()
leveneTest(growth~treatment,data=out%>%filter(type=="Individual")) #failed equal variances
wilcox.test(growth~treatment,data=out%>%filter(type=="Individual"))
leveneTest(growth~treatment,data=out%>%filter(type=="Aggregate")) #failed equal variances
wilcox.test(growth~treatment,data=out%>%filter(type=="Aggregate"))

my_tag<-c("p<0.001","p<0.001")

growth<-ggplot(out)+
  geom_hline(yintercept=0,linetype="dashed")+
  geom_boxplot(aes(treatment,growth,fill=treatment),outlier.colour = NA)+
  geom_point(aes(treatment,growth))+
  facet_wrap(~type)+
  theme_classic(base_size=8)+
  scale_fill_manual(values=c("gray","#6F3C4B"),guide="none")+
  ylab("Growth (percent change)")+
  xlab("Treatment");growth
  
tag_facet2 <- function(p, open = "(", close = ")", tag_pool = letters, x = -Inf, y = Inf, 
                       hjust = -0.5, vjust = 1.5, fontface = 2, family = "", ...) {
  
  gb <- ggplot_build(p)
  lay <- gb$layout$layout
  tags <- cbind(lay, label = paste0(open, tag_pool[lay$PANEL], close), x = x, y = y)
  p + geom_text(data = tags, aes_string(x = "x", y = "y", label = "label"), ..., hjust = hjust, 
                vjust = vjust, fontface = fontface, family = family, inherit.aes = FALSE)
}

growth<-tag_facet2(growth, 
                     x = 1, y = 2, 
                     vjust = 0, hjust = 0,
                     open = "", close = "",
                     fontface = 'italic',
                     size = 2,
                     tag_pool = my_tag);growth

quartz(w=3.5,h=2)



