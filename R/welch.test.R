require("dplyr")

welch.test <- function(.data, rate.col, time.col)
{
  rate.col <- enquo(rate.col)
  time.col <- enquo(time.col)
  
  fr <- .data %>% mutate(rate.std = abs(!!rate.col)/sqrt(!!time.col), time.std = sqrt(!!time.col))

  done <- summary(lm(data = fr, formula = rate.std ~ time.std))$coefficients[8] >=0.05
  print(done)
  
  while(!done)
  {
    min.time <- min(fr %>% pull(time.std))
    fr <- fr %>% filter(time.std > min.time)
    done <- summary(lm(data = fr, formula = rate.std ~ time.std))$coefficients[8] >=0.05
    print (done)
  }

  return (fr)
}

welch.test.2 <- function(.data, rate.col, time.col)
{
  rate.col <- enquo(rate.col)
  time.col <- enquo(time.col)
  
  fr <- .data %>% mutate(rate.std = abs(!!rate.col)/sqrt(!!time.col), time.std = sqrt(!!time.col))
  print(summary(lm(data = fr, formula = rate.std ~ time.std)))$coefficients
  done <- summary(lm(data = fr, formula = rate.std ~ time.std))$coefficients[8] >=0.05
  print(done)
  
  while(!done)
  {
    min.time <- min(fr %>% pull(time.std))
    fr <- fr %>% filter(time.std > min.time)
    done <- summary(lm(data = fr, formula = rate.std ~ time.std))$coefficients[8] >=0.05
    print (done)
  }
  
  return (fr)
}

means.test <- function(.data, first.col, last.col)
{
  first.col <- enquo(first.col)
  last.col <- enquo(last.col)
  
  fr <- .data %>% rowwise() %>% mutate(mean.col = mean(c(!!first.col, !!last.col)), d = abs(!!first.col - !!last.col))

  done <- summary(lm(data = fr, formula = mean.col ~ d))$coefficients[8] >=0.05
  print(done)
  
  while(!done)
  {
    min.mean <- min(fr %>% pull(mean.col))
    fr <- fr %>% filter(mean.col > min.mean)
    l <- summary(lm(data = fr, formula = mean.col ~ d))
    done <- l$coefficients[8] >=0.01

    print (done)
  }
  
  return (fr %>% select(!c(mean.col,d)))
}