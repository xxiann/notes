---
description: codes for all the various plots
---

# ggplot2

## standard ggplot colours

```r
ggplotColours <- function(n = 6, h = c(0, 360) + 15){
  # n: The number of colors to generate and return. Default is 6.
  if ((diff(h) %% 360) < 1) h[2] <- h[2] - 360/n
  hcl(h = (seq(h[1], h[2], length = n)), c = 100, l = 65)
}
```



