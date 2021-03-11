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
#calculate the distance (in bp) between consecutive positions
df = df %>%
  group_by(chr) %>%
  mutate(distance = start - lag(start, default = start[1]))
​
#shift the distance column up one row
shift <- function(x, n){
  c(x[-(seq(n))], rep(NA, n))
}
​
df$distance <- shift(df$distance, 1)
​
#filter those with a distance <=500 bp (i.e. the het must have another base change within 500 bp)
df = df %>% 
  group_by(chr) %>% 
  filter(distance <= 500)
df = df %>%
  group_by(chr) %>%
  mutate(distance2 = start - lag(start, default = start[1]))
​
#filter only if there are 5 consecutive rows of distance <=500 bp (i.e. the het must have 5 base changes within 500 bp)
r <- with(with(df, rle(distance2<=500)),rep(lengths,lengths))
df$het <- with(df,distance2<=500) & (r>5)
df$het <- shift(df$het, 1)
​
#filter the row if it contains a het
het_df <- filter(df, het == "TRUE")
​
###determine the max and min coordinates of a region with a distance <= 500 bp
​
#get the min and max coordinates
het_df = het_df %>% 
  group_by(chr) %>% 
  filter((distance2 >= 500) | lead(distance2 >= 500) | (row_number() >= (n())) | (row_number() == 1) )
​
#shift the end column to have the min and max coordinates on one row
het_df$end <- shift(het_df$end, 1)
​
#take every other row (1, 3, etc.)
het_df = het_df %>%
  group_by(chr) %>% 
  filter(row_number() %% 2 == 1)
​
#take the differences of the min and max coordinates
het_df = het_df %>%
  group_by(chr) %>%
  mutate(het_length = end - start)
​
#remove those with negative lengths
het_df = het_df %>%
  select(chr, start, end, het_ratio, het_length) %>% 
  filter(het_length > 0)
​
#filter if the het ratio is >= 10%
het_df_filtered = het_df %>% 
  group_by(chr) %>% 
  filter(het_ratio >= 10)
​
#print the table
write.table(het_df_filtered, "all.hets.winnowmap.greaterThan10.tbl", row.names = F, quote = F, sep="\t")
