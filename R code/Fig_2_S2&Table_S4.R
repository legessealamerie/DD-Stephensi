### library
library(ggplot2)
library(ggfortify)
library(tidyverse)
library(ggprism)
library(scales)
library(ggpubr)
library(dplyr)
library(palmerpenguins)

theme_set(theme_prism()+ theme(legend.position = "bottom")) # theme for ggplot

urb=read.csv("data/indiv_data_all_04Jun.csv",header=T,na.strings ="")
names(urb)
urb$rdtpos=as.factor(urb$rdtpos)
urb$micpos=as.factor(urb$micpos)
urb$casecat=as.factor(urb$casecat)
urb$casexp=as.factor(urb$casexp)
urb$site=as.factor(urb$site)
urb=urb%>%filter(pfqpos!="Negative")
urb$pfqden=urb$pfqden/1000
### summary table 
table(urb$rdtpos)
table(urb$micpos)
table(urb$casecat)
table(urb$site)

cols = c( "navy", "darkturquoise")
#cols = c( "tomato","tomato4")
#cols = c( "skyblue1","red1")
summary(urb$micpos)
summary(urb$rdtpos)
### distribution of pf18S based on  RDT
means <- urb %>% group_by(rdtpos) %>% 
  summarise(mean_pf18s = 10^mean(log10(pfqden), na.rm=TRUE))
p1 <- ggplot() +
  geom_density(data = urb, aes(x = pfqden, fill = rdtpos, color = rdtpos),
               alpha = 0.2, na.rm=T ,linewidth = 0.8) +
  geom_vline(data = means, aes(xintercept = mean_pf18s, color = rdtpos),
             linewidth = 1.3) +
  scale_x_log10(label =comma_format(big.mark = ""))+
  scale_fill_manual("  ", values = cols) +
  scale_color_manual(" ", values = cols) +
  labs(x = "Parasite per μL", y = "Density") +
 # ggtitle("Distribution of Pf18S based on RDT status") +
  ggtitle("") +guides(fill = FALSE)

# distribution of pf18S based on microscopy
# filtering excluding  micpos not available variables
mic=urb%>%filter(micpos!="NA")
means <- mic %>% group_by(micpos) %>% 
  summarise(mean_pf18s = 10^mean(log10(pfqden), na.rm=TRUE))
p2 <- ggplot() +
  geom_density(data = mic, aes(x = pfqden, fill =micpos, color = micpos),
               alpha = 0.2 ,linewidth = 0.8) +
  geom_vline(data = means, aes(xintercept = mean_pf18s, color = micpos),
             linewidth = 1.3) +
  scale_x_log10(label =comma_format(big.mark = ""))+
  scale_fill_manual(" ", values = cols) +
  scale_color_manual(" ", values = cols) +
  labs(x = "Parasite per μL", y = "Density") +
  #ggtitle("Distribution of Pf18S based on Microscopy status") +
  ggtitle(" ") +guides(fill = FALSE)

ggarrange(p1, p2, common.legend = TRUE, legend = "bottom")
ggsave("plots/RDT_Mic_den.tiff", dpi=300, height=5, width=12.3, units="in", compression = "lzw")

## t test 
t.test(pfqden ~ micpos, data = mic)

# distribution of pf18S based on study site

means <- urb %>% group_by(site) %>% 
  summarise(mean_pf18s = 10^mean(log10(pfqden), na.rm=TRUE))
summary(urb$site)
p3 <- ggplot() +
  geom_density(data = urb, aes(x = pfqden, fill = site, color = site),
               alpha = 0.2, na.rm=TRUE ,linewidth = 0.8) +
  geom_vline(data = means, aes(xintercept = mean_pf18s, color = site),
             linewidth = 1.3) +
  scale_x_log10(label =comma_format(big.mark = ""))+
  scale_fill_manual("Study site ", values = cols) +
  scale_color_manual("Study site  ", values = cols) +
  labs(x = "Parasite per μL", y = "Density") +scale_y_continuous(limits = c(0, 0.5), breaks = seq(0, 0.5, by =0.1))+
  ggtitle("")+guides(fill = FALSE)
ggsave("plots/site_new_dens_a.tiff", dpi=300, height=5, width=7, units="in", compression = "lzw")


# distribution of pf18S based on  case control category 

means <- urb %>% group_by(casecat) %>% 
  summarise(mean_pf18s = 10^mean(log10(pfqden), na.rm=TRUE))
table(urb$casecat)

cols = c( "navy", "darkturquoise","goldenrod1","olivedrab4")
#cols = c( "goldenrod1","olivedrab4")
summary(urb$site)
p4 <- ggplot() +
  geom_density(data = urb, aes(x = pfqden, fill = casecat, color = casecat),
               alpha = 0.2, na.rm=TRUE ,linewidth = 0.8) +
  geom_vline(data = means, aes(xintercept = mean_pf18s, color = casecat),
             linewidth = 1.3) +
  scale_x_log10(label =comma_format(big.mark = ""))+
  scale_fill_manual("casecat ", values = cols) +
  scale_color_manual("casecat  ", values = cols) +
  labs(x = "Parasite per μL", y = "Density") +theme_prism()+guides(fill = FALSE)
ggsave("plots/cascat_new_dens.tiff", dpi=300, height=5, width=9, units="in", compression = "lzw")

####


library(stats)
library(agricolae)
#### median comparison  for university
univ$logpf=log10(univ$pfqden+0.1)
city$logpf=log10(city$pfqden+0.1)
out<-with(univ,Median.test(logpf,casecat,console=FALSE))
summary(univ$logpf)
z<-bar.err(out$medians,variation = "range",ylim=c(0,5),
           space=2,border=4,col=3,density=10,angle=45,ylab="Pf18S copies per μL (Log)",main="University")
# median 
out<-with(univ,Median.test(pfqden,casecat,console=FALSE,group=FALSE))
print(out$comparison)
###City
urb$logpfde=log10(urb$pfqden+0.1)
out<-with(urb,Median.test(logpfde,casecat,console=FALSE))
summary(urb$logpfde)
z<-bar.err(out$medians,variation = "range",ylim=c(0,5),
           space=2,border=4,col=3,density=10,angle=45,ylab="Pf18S copies per μL (Log)",main="All")
# median 
out<-with(urb,Median.test(pfqden,casecat,console=FALSE,group=FALSE))
print(out$comparison)
####


urb=urb%>%filter(pfqpos!="Negative")
urb$pfqden=urb$pfqden/1000
means <- urb %>% group_by(casecat) %>% 
  summarise(mean_pf18s = 10^mean(log10(pfqden), na.rm=TRUE))
p2 <- ggplot() +
  geom_density(data = urb, aes(x = pfqden, fill =casecat, color = casecat),
               alpha = 0.2 ,linewidth = 0.8) +
  geom_vline(data = means, aes(xintercept = mean_pf18s, color = casecat),
             linewidth = 1.3) +
  scale_x_log10(label =comma_format(big.mark = ""))+
  scale_fill_manual("casecat ", values = cols) +
  scale_color_manual("casecat", values = cols) +
  labs(x = "Pf18S copies per μL", y = "Density") +
  ggtitle("City") ++guides(fill = FALSE)+
  theme_prism()

ggarrange(p1, p2, common.legend = TRUE, legend = "bottom")
ggsave("plots/site_positive_only.tiff", dpi=300, height=5, width=12.3, units="in", compression = "lzw")

#### mean and confidence interval  for RDT , Micrscopy ,site and Case control  Category
means <- urb %>% group_by(rdtpos) %>% 
  summarise(mean_pfqpcr18s = 10^mean(log10(pfqden), na.rm=TRUE),
            pfqpcr18s_lower = 10^t.test(log10(pfqden))$conf.int[1],
            pfqpcr18s_upper = 10^t.test(log10(pfqden))$conf.int[2])

means <- mic %>% group_by(micpos) %>% 
  summarise(mean_pfqpcr18s = 10^mean(log10(pfqden), na.rm=TRUE),
            pfqpcr18s_lower = 10^t.test(log10(pfqden))$conf.int[1],
            pfqpcr18s_upper = 10^t.test(log10(pfqden))$conf.int[2])

means <- urb %>% group_by(casecat) %>% 
  summarise(mean_pfqpcr18s = 10^mean(log10(pfqden), na.rm=TRUE),
            pfqpcr18s_lower = 10^t.test(log10(pfqden))$conf.int[1],
            pfqpcr18s_upper = 10^t.test(log10(pfqden))$conf.int[2])

means <- urb %>% group_by(site) %>% 
  summarise(mean_pfqpcr18s = 10^mean(log10(pfqden), na.rm=TRUE),
            pfqpcr18s_lower = 10^t.test(log10(pfqden))$conf.int[1],
            pfqpcr18s_upper = 10^t.test(log10(pfqden))$conf.int[2])
####

## t test  
t.test(pfqden ~ micpos, data = mic)
t.test(pfqden ~ rdtpos, data = urb)
t.test(pfqden ~ site, data = urb)
table (urb$pfqpos,urb$rdtpos)
table (urb$pfqpos,urb$micpos)
