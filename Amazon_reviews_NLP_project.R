#Amazon Musical Instruments Review Analysis 

#set the working directory
setwd("C:/Users/shaim/Downloads")

#used the read_csv function to read the data 
library(readr)
amazon <- read_csv('Musical_instruments_reviews.csv')

# obtain access to the placed custom functions in 
#corpus_functions_2.R
source("corpus_functions_2.R")

#create a new variable that contains only the reviews and the reviewerID
reviews_df <- amazon[, c("reviewerID", "reviewText")]

#check how many reviews are there in total 

length(reviews_df$reviewText) 

#check missing rows

sum(is.na(reviews_df$reviewText))

#remove rows with missing values

reviews_df <- reviews_df[!is.na(reviews_df$reviewText), ]

#tokenize the reviews so they can be lemmatized 

tokens <- tokenize(reviews_df$reviewText)

library(textstem)

lemma_tokens <- lemmatize_words(tokens)

#shows the number of unique tokens 
length(unique(tokens))

#shows the number of unique lemmatized tokens 
length(unique(lemma_tokens))

#shows the frequency of each lemmatized unique token
sort(table(lemma_tokens), decreasing = TRUE)

install.packages("textstem")

library(textstem)

#choose chunk size and apply it to each review using loops
# the chunked reviews will be stored in a new data frame called
#chunked_reviews_df

chunk_size <- 10  # Words per chunk
chunked_reviews_df <- data.frame(reviewerID = character(), review_chunk = character(), stringsAsFactors = FALSE)

# Loop over each review
for(i in 1:nrow(reviews_df)) {
  review <- reviews_df$reviewText[i]
  reviewerID <- reviews_df$reviewerID[i]
  
  # Tokenize the review (split into words)
  word_v <- tokenize(review)
  
  # Lemmatize the words
  lemmatized_words <- lemmatize_words(word_v)
  
  # If the review is longer than chunk_size, we chunk it
  if(length(lemmatized_words) > chunk_size) {
    # Split the review into chunks based on the chunk_size
    x <- seq_along(lemmatized_words)
    chunks_l <- split(lemmatized_words, ceiling(x / chunk_size))
    
    # Merge small chunks with the previous chunk (if last chunk is smaller than half of chunk_size)
    if(length(chunks_l[[length(chunks_l)]]) <= chunk_size / 2) {
      chunks_l[[length(chunks_l) - 1]] <- c(chunks_l[[length(chunks_l) - 1]], chunks_l[[length(chunks_l)]])
      chunks_l[[length(chunks_l)]] <- NULL  # Remove the small chunk
    }
    
    # Combine chunks back into strings
    chunk_strings_l <- lapply(chunks_l, paste, collapse = " ")
    
    # Append each chunk with the reviewerID
    for(chunk in chunk_strings_l) {
      chunked_reviews_df <- rbind(chunked_reviews_df, data.frame(reviewerID = reviewerID, review_chunk = chunk))
    }
  } else {
    # If the review is smaller than chunk_size, no chunking is done, just keep the full review
    chunked_reviews_df <- rbind(chunked_reviews_df, data.frame(reviewerID = reviewerID, review_chunk = paste(lemmatized_words, collapse = " ")))
  }
}

# View the first few rows of the chunked reviews
head(chunked_reviews_df)

# install the R Mallet Package
if(!require(mallet)) install.packages("mallet")
library(mallet)

mallet_instances <- mallet.import(
  chunked_reviews_df$reviewerID,
  chunked_reviews_df$review_chunk,
  "stoplist.csv",
  FALSE,
  token.regexp = "[\\p{L}']+" # Keep letters
  )

# choose number of topics 
topic_model <- MalletLDA(num.topics = 3)

# load the documents
topic_model$loadDocuments(mallet_instances)

# the corpus vocabulary is available as a character vector
vocabulary <- topic_model$getVocabulary()
class(vocabulary)
length(vocabulary)
head(vocabulary)
vocabulary[1:50]

# The word frequencies are also available
#  as a data frame with a row for each unique
#  word type in the corpus
word_freqs <- mallet.word.freqs(topic_model)
names(word_freqs) 
head(word_freqs)

#adjust the optimizier which is an optional step
topic_model$setAlphaOptimization(40, 80)

#set the iteration to 300 
topic_model$train(300)

topic_words_m <- mallet.topic.words(
  topic_model,
  smoothed = TRUE,
  normalized = TRUE
)

# Since normalization is set to TRUE, the values in
#  each topic (row) are converted to proportions
#  that sum to one which can be seen using rowSums
rowSums(topic_words_m)

# if normalization is FALSE, the contents are
#  counts of occurrences (integers) of that word type
topic_words_m[1:3, 1:3]

# use the vocabulary to name the columns
vocabulary <- topic_model$getVocabulary()
colnames(topic_words_m) <- vocabulary
topic_words_m[1:3, 1:3]

# investigate community violating words
keywords <- c("attack", "abuse")
topic_words_m[, keywords]

imp_row <- which(
  rowSums(topic_words_m[, keywords]) == 
    max(rowSums(topic_words_m[, keywords]))
)
imp_row


important_keywords <- topic_words_m[imp_row, keywords]

# Create a bar plot for these keywords
barplot(important_keywords,
        names.arg = keywords,  # Use the keywords as names
        main = "Keyword Occurrences in Important Topic",
        xlab = "Keywords",
        ylab = "Frequency",
        col = "lightblue")


##create word cloud and barplot for each topic and show 100 words
n_topics <- 3 # Define the number of topics
num_top_words <- 100  # Number of top words to show for each topic

if(!require(wordcloud)) install.packages("wordcloud")
library(wordcloud)

# Loop through each topic to get the top words
for (i in 1:n_topics) {
  # Get the top words for topic i
  top_words <- mallet.top.words(topic_model, topic_words_m[i, ], num_top_words)
  
  # Create a word cloud for the current topic
  wordcloud(
    words = top_words$term, 
    freq = top_words$weight, 
    scale = c(3, 0.5), 
    rot.per = 0.3, 
    random.order = FALSE, 
    main = paste("Top Words for Topic", i)
  )
  
  # create a barplot for the current topic
  barplot(
    top_words$weight, 
    names.arg = top_words$term, 
    las = 2, 
    main = paste("Top Words for Topic", i), 
    col = "steelblue"
  )
}

# create dispersion plot of the term "guitar" 

n_time_v <- seq(from = 1, to = length(lemma_tokens))
guitar_v <- which(lemma_tokens == "guitar")
g_count_v <- rep(NA, times = length(n_time_v))
g_count_v[guitar_v] <- 1

plot(
  g_count_v,
  main = "Dispersion Plot of 'guitar' in reviews",
  xlab = "Review Sequence",
  ylab = "guitar",
  type = "h",
  ylim = c(0, 1), yaxt = 'n'
)

# create dispersion plot of the term "music" 

n_time_v <- seq(from = 1, to = length(lemma_tokens))
music_v <- which(lemma_tokens == "music")
m_count_v <- rep(NA, times = length(n_time_v))
m_count_v[music_v] <- 1

plot(
  m_count_v,
  main = "Dispersion Plot of 'music' in reviews",
  xlab = "Review Sequence",
  ylab = "music",
  type = "h",
  ylim = c(0, 1), yaxt = 'n'
)

###################################################################


