---
description: Uses what is already in utils, base
---

# System functions

## Calculating elapsed times

```r
start_time <- Sys.time()
end_time <- Sys.time()
time_taken <- end_time - start_time
print(paste("Time taken: ", time_taken))
```

## Inserting a progress bar

```r
## initialising
pb <- txtProgressBar(min = 0, max = n, style = 3, width = 50, char = "=")

## in a loop
for (i in 1:n){
    ## some code
    setTxtProgressBar(pb, i)
}

close(pb)
```

## Saving files

```r
```











