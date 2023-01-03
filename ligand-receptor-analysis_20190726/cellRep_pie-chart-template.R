# 
# 
# df <- data.frame(
#   group = c("Male", "Female", "Child"),
#   value = c(25, 25, 50)
# )
# head(df)
# 
# library(ggplot2)
# # Barplot
# bp<- ggplot(df, aes(x="", y=value, fill=group))+
#   geom_bar(width = 1, stat = "identity")
# bp
# 
# pie <- bp + coord_polar("y", start=0)
# pie
# 
# pie + scale_fill_brewer(palette="Dark2")


my_piePlot <- function(df) {
  library(ggplot2)
  # Barplot
  bp<- ggplot(df, aes(x="", y=perc, fill=gross_labels))+
    geom_bar(width = 1, stat = "identity", size=0.05, color="#b8bbc2")
  bp
  
  pie <- bp + coord_polar("y", start=0)
  pie <- pie + theme_bw() + theme(panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank(),
                       panel.border = element_blank(),
                       panel.background = element_blank(),
                       legend.title=element_blank(), axis.title = element_blank())
  return(pie)
}
