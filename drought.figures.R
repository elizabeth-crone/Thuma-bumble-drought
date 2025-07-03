setwd("C:")
dat = read.csv("palmer index coastal MA monthly.csv")
head(dat)
summary(dat)
dim(dat)

dat = dat[dat$year > 2017,]

dry.val = round(min(dat[,2:13], na.rm = T),1)-0.1
wet.val = round(max(dat[,2:13], na.rm = T),1)+0.1
PDI.vals = seq(dry.val, wet.val, 0.1)

col1 = "goldenrod1"
col2 = "blue"

pal <- colorRampPalette(c(col1, col2))
nvals = length(PDI.vals)
pal2 = pal(nvals)

#pdf("Figure 1.pdf", height = 9, width = 5)
png("Figure 1.png", units="in", width=5, height=9, res=300)
par(mfrow = c(7,1), mar = c(3,5,0,1)+0.1, oma = c(1,1,1,1)) # default mar = c(5, 4, 4, 2) + 0.1
for(j in 2018:2024){
use.year = j
use.dat = dat[dat$year == use.year,2:13]
my.mins = array()
my.cols = array()
for(i in 1:length(use.dat)){
  my.mins[i] = max(which(PDI.vals < as.numeric(round(use.dat,1)[i])))
  my.cols[i] = pal2[my.mins[i]]
}

avg = mean(t(use.dat)[6:8])
col.avg = pal2[max(which(PDI.vals < round(avg,1)))]

plot(4:9, use.dat[4:9], type = "l", lwd = 2, ylim = c(dry.val, wet.val), xlab = "", ylab = "", xaxt = "n", cex.axis = 1.25)
points(4:9, use.dat[4:9], col = ifelse(use.dat[4:9]<=0, "goldenrod1","royalblue3"), pch = 19, cex = 2)
axis(side = 1, at = 4:9, labels = c("Apr", "May", "Jun", "Jul", "Aug", "Sep"), cex.axis = 1.25)
if(j == 2024) mtext(side = 1, line = 2.5, "Month", cex = 0.9)
if(j == 2021) mtext(side = 2, line = 2.5, "Palmer Drought Index", cex = 0.9)
abline(h = 0, lwd = 2, lty = "dotted")
text(x = 8.75, y = 3.25, use.year)
points(6:8, rep(avg,3), type = "l", lwd = 3, col = ifelse(avg<=0, "goldenrod1","royalblue3"))
}
dev.off()


library(ggplot2)
library(dplyr)
library(lubridate)
dat = read.csv("pmdi monthly coastal MA 1960.2024.csv")
head(dat)

jun.oct <- dat%>% 
  mutate(pmdi = case_when(month>10 ~ 0, TRUE ~ pmdi)) %>% #adding zeroes for spacing in figure
  filter(month==6|month==7|month==8|month==9|month==10|month==11|month==12) %>% 
  mutate(color = if_else(pmdi<0,"dry","wet"),
         year.month = make_date(year=year,month=month),
         month.name=month(year.month, label=TRUE)) %>% 
  filter(year<2025)

#pdf("Figure 2.pdf", height = 6, width = 9)
png("Figure 2.png", units="in", width=9, height=6, res=300)
ggplot(jun.oct, aes(x=year,y=pmdi,alpha=month.name,fill=color))+
  geom_bar(stat="identity",position="dodge")+
  scale_fill_manual(name = 'PMDI Indication', values = (c("darkorange3",'royalblue3')), labels=c("Drought Month", "Non-Drought Month"))+
  scale_alpha_manual(name="",values = c(0.6,0.7,0.8,0.9, 1,0,0),labels=c("June","July","August","September","October","",""))+
  geom_vline(xintercept = 2017.5, color="goldenrod3", size = 0.7)+
  geom_hline(yintercept=-3, color = "red", size=0.5)+
  geom_text(aes(x=2022.5, label="Study \nPeriod", y=5), colour="black", angle=0, size=4, fontface="italic", show.legend=F)+
  labs(x="Year",y="Palmer Modified Drought Index")+
  theme(panel.background = element_rect(fill="white"),
        axis.line.x = element_line(color="black", size = 0.5),
        axis.line.y = element_line(color="black", size = 0.5),
        panel.border = element_rect(colour = "grey", fill=NA, size=1),
        axis.text.x = element_text(angle = 45, hjust=1, size=14),
        axis.text.y = element_text(size=14),
        axis.title.x = element_text(size=16, colour = "black"),    
        axis.title.y = element_text(size=16, colour = "black"),
        legend.title = element_text(size=16, colour = "black"),
        legend.text = element_text(size=14))+
  scale_x_continuous(breaks = scales::pretty_breaks(n = 24))
dev.off()