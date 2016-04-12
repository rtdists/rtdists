require(dplyr)

#rr98 <- read.table("./inst/extdata/rr98-data", header = FALSE)
rr98 <- read.table(system.file("extdata", "rr98-data", package = "rtdists"), header = FALSE)

colnames(rr98) <- c("id", "session", "source", "strength", "instruction", "response", "rt", "trial")

make_block <- function(x) {
  out <- vector("numeric", length(x))
  block <- 1
  out[1] <- block
  for (i in 2:length(x)) {
    if (x[i-1] > x[i]) block <- block+1
    out[i] <- block
  }
  out
}

rr98 <- rr98[rr98$trial != 0,] # remove breaks between blocks
rr98 <- as.data.frame(rr98 %>% group_by(id, session) %>% mutate(block = make_block(trial)) %>% ungroup)  # add block variable
rr98 <- rr98[rr98$response != 2,] # remove 'other' responses
rr98 <- rr98[rr98$session != 1,] # remove first session
rr98 <- rr98[!(rr98$block == 1 & rr98$trial <= 20),] # remove first 20 trials per session
rr98 <- rr98[rr98$trial != 1,] # remove breaks between blocks
rr98$outlier <- ifelse(rr98$rt < 200 | rr98$rt > 2500, TRUE, FALSE)

rr98$source <- factor(rr98$source, levels = 1:2, labels = c("dark", "light"))
rr98$instruction <- factor(rr98$instruction, 0:1, c("speed", "accuracy"))
rr98$response <- factor(rr98$response, 0:1, c("dark", "light"))
rr98$rt <- rr98$rt/1000
rr98$block <- as.integer(rr98$block)
rr98$correct <- rr98$source == rr98$response
rr98$response_num <- ifelse(rr98$response == "dark", 1L, 2L)

rr98 <- rr98[,c("id", "session", "block", "trial", "instruction", "source", "strength", "response", "response_num", "correct", "rt", "outlier")]
rownames(rr98) <- NULL

str(rr98)

use_data(rr98, compress = "xz", overwrite = TRUE)
