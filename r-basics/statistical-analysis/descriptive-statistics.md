---
description: with reference to NUS ST2137
---

# descriptive statistics

## descriptive statistics

<table><thead><tr><th width="194">Term</th><th>Definition</th></tr></thead><tbody><tr><td>Parameter</td><td>numerical summary of the population - is unknown</td></tr><tr><td>Statistic</td><td>summary of a sample taken from the population, computed from the sample to make inferences about a population parameter<br>- Descriptive statistics: numerical and graphical summaries<br>- Inferential statistics: like hypothesis testing</td></tr></tbody></table>

statistics - collected from random sample from the population that is representative of the population

* **certain assumptions are made** about the underlying population distribution to make calculations easier, but they are not necessarily true
  * can be affected by outliers, which causes the sample distribution to depart from underlying distribution assumptions
* if the statistical methods are **not robust,** the conclusions derived from the samples using these statistical methods **might not be reliable**

### numerical summaries

1. #### Location

```r
length(mark)
summary(mark)
mean(mark)
median(mark)
quantile(mark)
```

2. #### Range

```r
range(mark)
var(mark)
sd(mark)
IQR(mark)
mark[order(mark)[1:5]] # The 5 smallest observations
size<-length(mark)   #sample size = 98
mark[order(mark)[(size-4):size]] #The 5 largest observations
```

* mean (and variance) is sensitive to the outlier while the median (and IQR) is not

3. #### Skewness

* skewed right or positively skewed: right tail is longer, bulk of the data is at the left
* skewed left or negatively skewed: left tail is longer, bulk of the data is at the right
* sample skewness

<figure><img src="../../.gitbook/assets/image (4).png" alt="" width="481"><figcaption></figcaption></figure>

* skewness value
  * \= 0 : perfectly symmetrical
  * < -1, > 1 : highly skewed
  * -1 < x < - 1/2, 1/2 < x < 1 : moderately skewed
  * -1/2 < x < 1/2 : approximately symmetrical
* if mean is the same/approximately the same as median, then the data are close to symmetric
* when mean > median = right/positively skewed - affected by the high values on the right tail
* when mean < median = left/negatively skewed - affected by the low values on the left tail

```r
# skewness
skew <- function(x){
    n <- length(x)
    m3 <- mean((x-mean(x))^3)
    m2 <- mean((x-mean(x))^2)
    sk = m3/m2^(3/2)*sqrt(n*(n-1))/(n-2)
    
    return(sk)
}
```

4. #### Shape of distribution | kurtosis

* numerical measures of shape of a distribution - how tall and sharp the central peak is relative to a standard bell curve
* higher values of kurtosis == higher, sharper peak
* lower values of kurtosis == lower, less distinct peak

<figure><img src="../../.gitbook/assets/image (5).png" alt="" width="428"><figcaption><p>`</p></figcaption></figure>

```r
# kurtosis`
kurt <- function(x){
    n <- length(x)
    m4 <- mean((x-mean(x))^4)
    m2 <- mean((x-mean(x))^2)
    kurt = (n-1)/((n-2)*(n-3))*((n+1)*m4/(m2^2)-3*(n-1))
    return(kurt)
}
```

5. #### Correlation

* quantifies the relationship between 2 quantitative variables
  * can be linear or non-linear
* let X and Y be variables from a set of n objects/people (ie same sample set)

<figure><img src="../../.gitbook/assets/image (6).png" alt="" width="382"><figcaption></figcaption></figure>

```r
cor(final,midterm)
```

### graphical summaries

1. #### Histogram

* uses bars to portray frequencies or relative frequencies of possible outcomes for a quantitative variable
* look for: overall pattern - are there gaps - one or more observation deviate from the rest? / whats the mode of distribution? / is the distribution symmetric or skewed? / any suspected outliers?
* gap ≠ outliers

<pre class="language-r"><code class="lang-r"># using base r
hist(mark)
hist(mark, breaks = 30) # breaks = number of bars

# freq = T - for frequency plotting; freq = F - for relative frequency plotting
hist(mark, freq=TRUE, main = paste("Histogram of mark"),
 xlab = "mark", ylab="frequency", axes = TRUE, col = 5) 
 
<strong># using ggplot
</strong>ggplot(data) + 
  geom_histogram(aes(x = mark), binwidth = 3, fill = "grey", color = "red")
</code></pre>

<figure><img src="../../.gitbook/assets/image (7).png" alt="" width="563"><figcaption></figcaption></figure>

2. #### Density plot

* plots of smoothed histograms

```r
# using base r
hist(mark, freq=FALSE, main = paste("Histogram of mark"),xlab = "mark", ylab="Probability", axes = TRUE, 
     col = "grey",nclass = 10)
x <- seq(0, 30, length.out=98)
y <-dnorm(x, mean(mark), sd(mark)) # creates a normal distribution based on mean and sd given
lines(x, y, col = "red") # creates a normal density cuve

lines(density(mark)) #density curve of the sample, doesnt observe if the sample is normal

# using ggplot
## to add a normal density curve to a histogram
p <- ggplot(data) +
    geom_histogram(aes(x = mark, y = ..density..),
                   binwidth = 3, fill = "grey", color = "black")
x <- seq(0, 30, length.out=98) # == seq(0, 30, 0.05)
y = dnorm(x, mean(mark), sd(mark)) # using the mean and sd to create a normal curve
df <- data.frame(x = x, y = y)
p + geom_line(data = df, aes(x = x, y = y), color = "red")
```

<figure><img src="../../.gitbook/assets/image (8).png" alt="" width="302"><figcaption></figcaption></figure>

3. #### Boxplot

* provides skeletal representation of a distribution - easy to identify outliers, min, max, IQR, Q1, Q2, Q3
  * max/min whiskers: Q3 +1.5IQR / Q1 - 1.5IQR
  * outliers = large outliers > max whiskers ; small outliers < min whiskers
  * extreme outliers = > Q3 + 3IQR ; < Q1 - 3IQR
*   is well suited for showing distributions for multiple variables - fast comparison across groups of the same variables

    <figure><img src="../../.gitbook/assets/image (10).png" alt="" width="237"><figcaption></figcaption></figure>

```r
# using based r
boxplot(mark, xlab = "Midterm mark")

# using ggplot
ggplot(data) + geom_boxplot(aes(y = mark))
```

4. #### QQ plot

*   to see if the data follow (approximately) a normal distribution or not

    <figure><img src="../../.gitbook/assets/image (11).png" alt="" width="299"><figcaption></figcaption></figure>
* plotting method: standardises the sample quantiles into \~N(0,1) and plot it against the theoretical quantiles of a N(0,1) distribution \[y = x]
  * if the sample quantile coincides = can say that there is evidence that the data came from a normal distribution
* when X: sample quantiles ; Y: theoretical quantiles
  * right tail < below straight line = observed > expected : longer right tail than normal
  * right tail > above straight line = observed < expected : shorter right tail than normal
  * left tail < below straight line = observed > expected : shorter left tail than normal
  * left tail > above straight line = observed < expected : longer left tail than normal

```r
# pch - fills up the dots plotted
# datax = TRUE - sample quantile
qqnorm(mark,datax = TRUE, ylab = "Sample Quantiles", xlab = "Theorical Quantiles", 
         main = "QQ Plot", pch = 20)
qqline(mark,datax = TRUE, col = "red")

qqnorm(mark, ylab = "Sample Quantiles", xlab = "Theorical Quantiles", 
         main = "QQ Plot", pch = 20)
qqline(mark, col = "red")

# using ggplot
ggplot(data) + 
  geom_qq(aes(sample = mark)) + 
  geom_qq_line(aes(sample = mark),color = "red")+ 
  coord_flip() ## flip x and y axis
```

<figure><img src="../../.gitbook/assets/image (9).png" alt="" width="563"><figcaption></figcaption></figure>

5. #### Scatterplot

* help visualise the assocation between 2 quantitative variables
* questions to ask: are there any relationship between 2 variables / if yes, are there any associations? / are there any unusual observations, departing from overall trend? / is the **variance** of y-variable stable when the value of the x-variable changes?

```r
plot(midterm,final, pch = 20)
```

<figure><img src="../../.gitbook/assets/image (12).png" alt="" width="501"><figcaption></figcaption></figure>

## robust statistics

* when the statistical method is insensitive to slight departures from assumptions
  *   **robustness of assumptions (\~ robustness of statistical methods)**

      > a statistical method is said to be robust with respect to a particular assumption if it performs adequately even when the assumption is modestly violated
* robustness can be measured by measures such as breakdown point, influence curve and gross error sensitivity
* example: 95% CI for the population mean
  * A 95% confidence interval for a population mean _μ_ is:
  * $$X̄ ± t_(n-1,0.975) * (s / √n)$$
    * where $$t_{n-1,0.975}$$  corresponds to the 0.975-quantile of a _t_-distribution with _n_-1 degrees of freedom.
  * Assumptions:
    * sample must be obtained by randomisation - not robust
    * distribution of the data must be approximately normal - robust (just need to ensure no extreme outliers or that sample size n is relatively large enough)

### robust estimation of location parameter

in comparison to sample mean

1. **trimmed mean**

* 100α % trimmed mean = discards the lowest 100α % and the highest 100α % and take the arithmetic mean of the remaining data
* recommended α = 0.1 - 0.2

```r
mean(x, trim = 0.2)
```

2. **winsorised mean**

* mean of trimmed and replaced data - by replacing extreme data values with less extreme values
* recommended α = 0.1 - 0.2

```r
winsor<-function(x, alpha = 0.2) { 
	n = length(x)
	xq = n * alpha
	x = sort(x)
	m = x[(round(xq)+1)]
	M = x[(n - round(xq))]
	x[which(x<m)] = m
	x[which(x>M)] = M
	return(c(mean(x),var(x))) 
} 
winsor(x, alpha = 0.2) 

library(MASS)
hubers(x, k= 0.84) # this gives the 20% Winsorized mean
```

3. M-estimates / Huber’s M-estimates

* MLE = maximum likelihood estimation = a method of estimating the parameters of a statistical model given observations, by finding the parameter values that maximize the likelihood of making the observations given the parameters
* $$\sum_{i=1}^{n} \rho(x_i - T)$$ , where T is the the estimator/minimizer and $$\rho$$ is a non-constant function that is meaningful

4. Tukey’s bisquare estimator
5. Humpe’s M-estimator

### robust estimator of scale

Sample standard deviation, _s_, is a commonly used estimator of population scale parameter _σ_.

* It is not robust and is sensitive to outliers and may not remain bounded when a single data point is replaced by an arbitrary number (i.e., become infinite).
* IQR is not a robust estimator of _σ_.
  * For normal distribution \~N(0,1), the SD _σ_ can be estimated by IQR/1.35.

1. **Median Absolute Deviation (MAD)** - most popular robust estimator of scale parameter.

* $$MAD = median_i(|x_i - median_j(x_j)|)$$
  * where the inner median, $$median_j(x_j)$$, is the median of _n_ observations and&#x20;
  * the outer median, $$median_i$$ is the median of the _n_ absolute values of the deviations about the median.
* For normal distribution where _Z_ \~ _N_(0, 1), SD _σ_ can be estimated by 1.4826 × MAD.

```r
median(abs(x - median(x))) #MAD
mad(x) # estimate of \sigma, = 1.4826*MAD
```

2. **Gini's mean difference**

* mean of all absolute mutual differences of any 2 observations of the sample.
* $$G = \frac{1}{\binom{n}{2}} \sum_{i<j} |x_i - x_j|$$
* For normal distribution, SD _σ_ can be estimated by $$f(x) = x * e^{2 pi i \xi x}$$ $$\sqrt{\pi} * \frac{G}{2}$$

&#x20;

