library(ggplot2)
library(directlabels)
library(plyr)
require(gridExtra)
require(scales)
library(RColorBrewer)
library(reshape2)
library(data.table)
library(dplyr)
library(tidyr)
library(zoo)
​
###change the working directory (2nd line) and saved file name (last line) before running###
​
#clear variables
rm(list=ls(all=TRUE))
​
#set working directory
setwd("~/Documents/Eichler Lab/Weekly plans/T2Tv1/Hets")
​
#load NucFreq bed file
df = fread(file.choose(), stringsAsFactors = FALSE, fill=TRUE, quote="", header=FALSE, skip=2)
cols =  c("chr", "start", "end", "first", "second")
colnames(df) <- cols
​
#determine the ratio of the first and second most common bases
df$het_ratio = round(df$second/(df$first+df$second)*100, 1)
​
#filter if the het ratio is >= 10%
df1 = df %>% 
  group_by(chr) %>% 
  filter(het_ratio >= 10)
​
#calculate the distance (in bp) between consecutive positions
df2 = df1 %>%
  group_by(chr) %>%
  mutate(distance = start - lag(start, default = start[1]))
​
#shift the distance column up one row
shift <- function(x, n){
  c(x[-(seq(n))], rep(NA, n))
}
​
df2$distance <- shift(df2$distance, 1)
​
#filter rows with a distance <=500 bp between positions (i.e. the het must have another base change within 500 bp)
df3 = df2 %>% 
  group_by(chr) %>% 
  filter(distance <= 500)
df3 = df3 %>%
  group_by(chr) %>%
  mutate(distance2 = start - lag(start, default = start[1]))
df3$distance2 <- shift(df3$distance2, 1)
​
#duplicate top row and change its distance to 501 bp (to get rows in register)
df4 = df3 %>% 
  group_by(chr) %>% 
  filter(row_number() <= 1) %>% 
  bind_rows(df3)
df5 = df4 %>%
  arrange(start, .by_group = TRUE) %>%
  mutate(distance2 = replace(distance2, row_number() == 1, 501))
​
#shift up the end column to get the range of the hets on one row
df5$end <- shift(df5$end, 1)
​
#filter only if there are 5 consecutive rows of distance <=500 bp (i.e. the het must have 5 base changes within 500 bp)
r <- with(with(df5, rle(distance2<=500)),rep(lengths,lengths))
df5$het <- with(df5,distance2<=500) & (r>=4)
​
#filter the row if it contains a het
df6 <- filter(df5, het == "TRUE")
​
###determine the max and min coordinates of a region with a distance <= 500 bp
​
#get the min and max coordinates
het_df = df6 %>%
  group_by(chr) %>%
  mutate(distance3 = start - lag(start, default = start[1]))
het_df$distance3 <- shift(het_df$distance3, 1)
het_df2 = het_df %>% 
  group_by(chr) %>% 
  filter((distance3 >= 500) | lag(distance3 >= 500) | (row_number() >= (n())) | (row_number() == 1) )
​
#shift the end column to have the min and max coordinates on one row
het_df2$end <- shift(het_df2$end, 1)
​
#take every other row (1, 3, etc.)
het_df3 = het_df2 %>%
  group_by(chr) %>% 
  filter(row_number() %% 2 == 1)
​
#take the differences of the min and max coordinates
het_df3 = het_df3 %>%
  group_by(chr) %>%
  mutate(het_length = end - start)
​
#remove those with negative lengths
het_df_filtered = het_df3 %>%
  select(chr, start, end, het_ratio, het_length) %>% 
  filter(het_length > 0)
​
#print the table
write.table(het_df_filtered, "all.hets.winnowmap.greaterThan10_2.tbl", row.names = F, quote = F, sep="\t")

