getwd()

rank_data = read.csv("rank.csv", header=T)

num_data = read.csv("scaled.csv", header=T)

str(num_data)

for (i in 1:nrow(num_data)) {
  num_data$unit_transformed[i] <- num_data$scale_unit[i] * 100
}

barplot(num_data$unit_transformed, main="Scaled # of Peptides Per Gene",
        xlab="Gene Names", ylab="High Scoring Peps Per Unit (Mult. x100%), col="lightblue2", names=num_data$name)


png(file="graph-it.png", width=600, height=600)

barplot(num_data$unit_transformed, main="Scaled # of Peptides Per Gene",
        xlab="Gene Names", ylab="High Scoring Peps Per Unit (Mult. x100%)", col="lightblue2", names=num_data$name)
dev.off()
