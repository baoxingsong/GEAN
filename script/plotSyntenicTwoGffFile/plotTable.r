library(ggplot2)
data = read.table("table")
ggplot(data=data, aes(x=V2, y=V4)) + geom_point() + facet_grid(V3~V1, scales="free", space="free")
