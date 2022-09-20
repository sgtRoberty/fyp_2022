#BAMM run evaluation function

# libraries required: coda, data.table, BAMMtools

# if mcmcFile or chainFile are NULL, then it will look for the most likely files

# returns effective size and Geweke z-score for the log likelihood and the number of shifts
# as well as frequency of hot/cold chain swapping, if multi-chain analysis. 

# ideal:
# effective size > 200
# Geweke z-score between -2 and 2
# chain swap frequency between 20 and 60%.

library(coda)
library(data.table)
library(BAMMtools)

bammCheck <- function(expectedNumberOfShifts = 1, burnin=0.15, 
                      mcmcFile = NULL, chainFile = NULL, 
                      plotToPDF = FALSE, title = NULL) {
  
  if (is.null(mcmcFile)) {
    files <- list.files(pattern='mcmc_out')
    if (length(files) == 1) {
      mcmcFile <- files
    } else {
      stop('Please specify an mcmc_out file. Several were detected.')
    }
  }
  if (is.null(chainFile)) {
    files <- list.files(pattern='chain_swap')
    if (length(files) == 1) {
      chainFile <- files
    } else if (length(files) == 0) {
      chainFile <- 'chain_swap'
      cat('\n\tNo chain swap file detected. Assuming single chain analysis.\n')
    } else {
      stop('Please specify an mcmc_out file. Several were detected.')
    }
  }
  
  if (class(mcmcFile) == 'character') {
    mcmc <- data.table::fread(mcmcFile, data.table = FALSE)
  }
  mcmc2 <- mcmc[floor(burnin * nrow(mcmc)):nrow(mcmc), ]
  
  # if BAMM was running with a single chain, there will be no chain swap file.
  chainswap <- data.table::fread(chainFile, data.table = FALSE)
  chainswap <- chainswap[(burnin * nrow(chainswap)):nrow(chainswap),]
  chainswap <- chainswap[which(chainswap$rank_1 == 1), ]
  success <- length(which(chainswap$swapAccepted == 1)) / nrow(chainswap)
  success <- round(success, digits = 2)

  #get prior distribution of shifts
  obsK <- seq(from=0, to = max(mcmc2[, "N_shifts"]), by = 1)
  prior <- sapply(obsK, prob.k, poissonRatePrior=1 / expectedNumberOfShifts)
  prior <- data.frame(N_shifts = obsK, prob = prior)
  
  #get posterior distribution of shifts
  posterior <- sapply(obsK, function(x) length(which(mcmc2[,'N_shifts'] == x))) / nrow(mcmc2)
  names(posterior) <- obsK
  posterior <- data.frame(N_shifts = names(posterior), prob=posterior)
  
  if (!plotToPDF) {
    dev.new(width = 10, height = 7)
  }
  par(mfrow = c(2, 2), mar = c(4, 4, 3, 3))
  plot(mcmc[, 1], mcmc[, 4], type = 'l', lwd = 0.5, xlab = 'generations', ylab = 'loglik')
  abline(v = max(mcmc[,1]) * burnin, lty = 2)
  
  if (!is.null(title)) {
    title(main = title)
  }
  
  #barplot(table(mcmc2$N_shifts), xlab='n shifts')
  barplot(prior[,2], names.arg = prior[,1], ylim = c(0, max(c(prior[,2], posterior[,2]))), border = 'black', col = 'light blue', xlab = 'n shifts')
  barplot(posterior[,2], add = TRUE, border = 'black', col = BAMMtools::transparentColor('red', 0.4), axes=FALSE)
  legend('topright', legend = c('prior','posterior'), fill = c('light blue', BAMMtools::transparentColor('red', 0.4)), bty = 'n', cex=1.5)
  
  plot(mcmc[,1], mcmc[,2], type='p', pch=20, cex=0.5, xlab='generations', ylab='n Events')
  abline(v=max(mcmc[,1]) * burnin, lty=2)
  
  plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n', xpd=NA)
  text(x=0.5, y=0.9, paste("effective size for log lik: ", round(coda::effectiveSize(mcmc2[,4]), 2), sep=''))
  text(x=0.5, y=0.7, paste("effective size for n Events: ", round(coda::effectiveSize(mcmc2[,2]), 2), sep=''))
  text(x=0.5, y=0.5, paste("Geweke z-score for log lik: ", round(coda::geweke.diag(coda::as.mcmc(mcmc2[,4]))$z, 2), sep=''))
  text(x=0.5, y=0.3, paste("Geweke z-score for n Events: ", round(coda::geweke.diag(coda::as.mcmc(mcmc2[,2]))$z, 2), sep=''))
  
  text(x=0.5, y=0.1, paste("successful swap frequency with cold chain: ", success, sep=''))

}

prob.k <- function(k, poissonRatePrior) {
  Denom <- (poissonRatePrior + 1) ^ (k + 1)
  Prob <- poissonRatePrior / Denom
  return(Prob)
}
