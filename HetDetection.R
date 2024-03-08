#This script filters the output by nucfreq and reports "regions where the second most common base was present in at least 10% of reads in at least 5 positions within a 500 bp region"

require(tidyr)
require(data.table)
require(dplyr)

# Check if any command line arguments are provided
if(length(commandArgs(trailingOnly = TRUE)) == 0) {
  stop("No filename provided. Please provide a filename as a command-line argument.")
}

# Get the filename from the command line
filename <- commandArgs(trailingOnly = TRUE)[1]


#load NucFreq bed file
df = read.table(filename, stringsAsFactors = FALSE, quote="", header=TRUE)
cols =  c("chr", "start", "end", "first", "second")
colnames(df) <- cols

#determine the ratio of the first and second most common bases
df$het_ratio = round(df$second/(df$first+df$second)*100, 1)

#filter if the het ratio is >= 10%
df = df %>% 
  group_by(chr) %>% 
  filter(het_ratio >= 10)

#calculate the distances betwen variants
df <- df %>%
  arrange(chr, start) %>%
  group_by(chr) %>%
  mutate(
    lag_distance = lead(start) - end
  )

#is distance smaller than 500 bp? if yes, set to true
df$closby<-df$lag_distance<=500 

#the last rows are going to have NA, since there are more variants than distances
#let's fill it with the value of the previous row
df <- df %>%
  fill(closby, .direction = "down")

#count how many times you get consecutive TRUE or FALSE (looking for the consecutive small distances)
#this is achieved using run-length encoding, and then expanding it so that the number of rows matches
df$rle<-rep(rle(df$closby)$lengths,rle(df$closby)$lengths) 

#assign a unique group id for each cluster, since the groups should not get mixed up
df$group<-rleid(paste0(df$closby,df$rle))

#only keep rows that belong to clusters
df<-df[!df$closby==FALSE,]

#only keep rows that cluster at least 5 variants
df<-df[df$rle>=5,] 

het_ratio_per_group <- df %>%
  group_by(group) %>%
  summarise(sum_second = sum(second),
            sum_first = sum(first)) %>%
  mutate(het_ratio = round(sum_second/(sum_first+sum_second)*100, 1))


df_final <- df %>%
  group_by(chr,group) %>%
  filter(row_number() == 1 | row_number() == n()) %>%
  summarise(start = first(start),
            end = last(end))

#calculate the lengths of het regions
df_final$het_length<-df_final$end-df_final$start 

#merge with the het frequency information
df_table<-merge(df_final,het_ratio_per_group,by="group")

#only keep start and end columns
df_table<-df_table[,c("chr","start","end","het_ratio","het_length")]
df_table<-na.omit(df_table)
write.table(df_table, paste0(filename,"all.hets.filtering.tbl"), row.names = F, quote = F, sep="\t")


