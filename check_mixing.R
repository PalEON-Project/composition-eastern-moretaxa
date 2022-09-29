outputDir="/var/tmp/paleon/composition-eastern-moretaxa/output"

num <- 10

keep <- 101:300
s2 <- s2keep <- NULL
for(i in 1:num) {
    load(file.path(outputDir, paste0('sigma2_eastern_0.2-', i, '.Rda')))
    s2 <- rbind(s2, sigma2store)
    s2keep <- rbind(s2keep, sigma2store[keep,])
}

par(mfrow=c(5,7))
for(i in 1:ncol(s2))
    ts.plot(s2keep[,i])

apply(s2keep, 2, coda::effectiveSize)




## check probs after burnin removal
