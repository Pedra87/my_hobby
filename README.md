install.packages("tidyverse")
data("iris")

#####............................ggplot2
library(ggplot2)

ggplot(data = iris, aes(x=iris$Sepal.Length, y=iris$Petal.Length, col=iris$Species))+ geom_point()
data
data = iris
df = data
2+2
# Library
library(ggplot2)
install.packages("hrbrthemes")
library(hrbrthemes)

# Dummy data
x <- LETTERS[1:20]
y <- paste0("var", seq(1,20))
data <- expand.grid(X=x, Y=y)
data$Z <- runif(400, 0, 5)

# Give extreme colors:
ggplot(data, aes(X, Y, fill= Z)) + 
  geom_tile() +
  scale_fill_gradient(low="white", high="blue") +
  theme_classic()

# Color Brewer palette
ggplot(data, aes(X, Y, fill= Z)) + 
  geom_tile() +
  scale_fill_distiller(palette = "RdPu") +
  theme_get()

# Color Brewer palette
ggplot(data, aes(X, Y, fill= Z)) + 
  geom_tile() +
  theme_bw()

############.own heat maps
library(readxl)
data <- read_excel("data/csrA-WT_sorted.xlsx")
 
# Color Brewer palette
ggplot(data, aes(X, Y, fill= Z)) + 
  geom_tile() +
  theme_bw()
