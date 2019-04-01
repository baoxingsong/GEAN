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



