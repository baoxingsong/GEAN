library(ggplot2)
data = read.table("M017_quota_longestpath_quota_true.table")
data$V5 = as.factor(data$V5)
png(file="M017_quota_longestpath_quota_true.png",  width=14000, height=4000)
ggplot(data=data, aes(x=V2, y=V4)) + geom_point(aes(group=V5, color=V5)) + facet_grid(V3~V1, scales="free", space="free")
dev.off()


data = read.table("M017_quota_longestpath_quota_false.table")
png(file="M017_quota_longestpath_quota_false.png",  width=14000, height=4000)
ggplot(data=data, aes(x=V2, y=V4)) + geom_point() + facet_grid(V3~V1, scales="free", space="free")
dev.off()


data = read.table("M017_dagChainer_true.table")
data$V5 = as.factor(data$V5)
png(file="M017_dagChainer_true.png",  width=14000, height=4000)
ggplot(data=data, aes(x=V2, y=V4)) + geom_point(aes(group=V5, color=V5)) + facet_grid(V3~V1, scales="free", space="free")
dev.off()


data = read.table("A1013_quota_dagChainer_true.table")
data$V5 = as.factor(data$V5)
data=data[order(data$V1, data$V2),]
data$V3 = factor(data$V3, levels=unique(data$V3))
png(file="A1013_quota_dagChainer_true.png",  width=2000, height=2000)
ggplot(data=data, aes(x=V2, y=V4)) + geom_point(aes(group=V5, color=V5)) + facet_grid(V3~V1, scales="free", space="free")+theme(legend.position = "none")+
    theme(panel.spacing = unit(0, "cm"), plot.margin = unit(c(0,0,0,0), "cm"), axis.text = element_blank(), axis.ticks = element_blank(), legend.title = element_blank())
dev.off()



data = read.table("A1013_quota_longestpath_quota_true.table")
data$V5 = as.factor(data$V5)
data=data[order(data$V1, data$V2),]
data$V3 = factor(data$V3, levels=unique(data$V3))
png(file="A1013_quota_longestpath_quota_true.png",  width=2000, height=2000)
ggplot(data=data, aes(x=V2, y=V4)) + geom_point(aes(group=V5, color=V5)) + facet_grid(V3~V1, scales="free", space="free")+theme(legend.position = "none")+
    theme(panel.spacing = unit(0, "cm"), plot.margin = unit(c(0,0,0,0), "cm"), axis.text = element_blank(), axis.ticks = element_blank(), legend.title = element_blank())
dev.off()


data = read.table("sorghum_A1013_quota_longestpath_quota_true.table")
data$V5 = as.factor(data$V5)
data=data[order(data$V1, data$V2),]
data$V3 = factor(data$V3, levels=unique(data$V3))
png(file="sorghum_A1013_quota_longestpath_quota_true.png",  width=2000, height=2000)
ggplot(data=data, aes(x=V2, y=V4)) + geom_point(aes(group=V5, color=V5)) + facet_grid(V3~V1, scales="free", space="free")+theme(legend.position = "none")+
    theme(panel.spacing = unit(0, "cm"), plot.margin = unit(c(0,0,0,0), "cm"), axis.text = element_blank(), axis.ticks = element_blank(), legend.title = element_blank())
dev.off()

