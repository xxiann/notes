# basics of r

note: **R passes arguments by value**

{% embed url="https://blog.penjee.com/wp-content/uploads/2015/02/pass-by-reference-vs-pass-by-value-animation.gif" %}

## logical operators

* AND, &
* OR, |
* NOT, !
* &&, || - only operates on the first elements of the 2 vectors

## relational operators

*   < and > for strings - evaluated based on alphabetical order

    ```r
    "raining" < "raining dogs"
    [1] TRUE
    ```
* TRUE = 1, FALSE = 0
* r comparators can handle vectors - returns in vectors
  * r comparators can also compare vectors to vectors, comparing the elements one by one
  * returns logical vectors

## data types

* vector: sets of elements with the same mode
* matrix: set of elements in row and columns but with the same mode
* dataframe: set of elements in row and columns, but columns can be of different modes
* list: generalisation of a vector/ a collection of data objects

### creating vectors

* using `c()`
  * `number ← c(1,2,3,4,5)`
* using `numeric()`, `character()`, `logical()` etc. - creates a vector with all its elements being 0
  * `number ← numeric(3)`
* using `rep(a,b)` - replicate the item a by times
  * `number ← rep(2,3)`
  * `number ← rep(c(2,1),3:4)`
    * 2, 2, 2, 1, 1, 1, 1
* using `seq(from=a, to=b, by=c, length=d)` # d = no of elements
  * `number ← seq(2, 10, 2)`

### creating matrix

* using dim(v), takes in the values of vector v to fill the matrix by column
  * `dim(v) ← c(2,3) # 2 rows, 3 columns`
  * using `dim(v) = NULL` to convert the matrix back to vector again
* using matrix(v,r,c) - takes in values from vector v to create matrix with r rows by c columns; default fills by column
  * `v ← c(1:6)`
  * `m ← matrix(v, nrow = 2, ncol = 3, byrow = F)`
  * if provided with extra elements - fills the matrix with the first few elements then give warnings for the remaining elements not used
  * if provided with insufficient elements - fills the matrix will all the elements then repeats + warning message
* using `rbind()`, binding rows to matrix; using `cbind()`, binding columns to matrix
  * note: length of row/column should match each other

#### other functions for matrix

* `solve(m)` to return the inverse of m
* `t(m)` to return the transpose of m
* `m %*% n` to multiply the 2 matrix
* `diag(n)` to create diagonal matrix

### creating data frame

* rows = the observations/measurements
* columns = the variables
* `data.frame()`= converts matrix into dataframe
  * `names(df)`- if default, can be used to rename
* reading/importing data
  * `scan()` - read data into vector or list from console/file
  * `read.table()`, read.csv(), read.fwf()
    * file path format “C:/User/xiaox/… ” - backslash must be negated using \\\\
    * to navigate to previous folder using “../”
    * take note of header, sep, col.names, row.names
  * `read.fwf()` - allows the specification of width of data for each columns using width=c(…)
  * reading of other file formats like spss, S-plus, SAS using library(foreign)

#### to assess data frame

* `data ← read.table(…)`
* `attach(data)` - makes the variables in the df accessible by name within the R session
* `data[r, c]`
* using logical test, selects all the rows that fits this condition, `data[Gender == “M”, ]`

#### combining data frame

* using `cbind` and `rbind`, the variable, ie names of columns are not the same or do not have the same length in both dataframe, an error message will be displayed
* using merge, `merge(a, b, by= , all= F)` 2 dataframes by a common variable
  * if all=T, extra rows (unmatched rows) will be added as well

#### sorting dataframe

* in ascending order: `data[order(CA1), ]`
* in descending order: `data[rev(order(CA1)), c(1, 3:4)]` - & selecting only specific rows

## loops

* `while (condition) {expression}`
* `for (<variable> in <range>) {expression}`

## redirecting output

* sink() - useful to keep a copy of output in a file or when the contents of aan object or funciton that may be too big to display on screen

```r
sink("C:/file.txt")
<body>
sink() # all the output are stored in the file
```

* `cat()` prints the output shown in the console unless redirected by sink
  * \n - to start a new line
* `write.table()`, `write.csv()`

## functions

* max, min, sum, mean, median, range, var, cor(x,y), sort, rank
* user-defined functions: `name <- function (parameters) {function body}`
  * default - last line of the function will be returned unless `return()` is included
* `args()` - to find out the arguments required for the fn without going through `help()` or `?`
* r passes argument by value (not by reference)

## packages

* base : automatically installed
* using `install.packages()`
* load package = attach to search list - `search()` - R will go through the list whenever the fn is called
  * the functions are added into R session, and .GlobalEnv (user workspace)
* `require()` - works similarly to library() but returns a warning msg when a library is not installed instead of an error
  * good for loading libraries dynamically in a fn
* library(ggplot2). library(”ggplot2”) will both work

## lapply

* instead of using a for loop that goes through a list/vector for a function
* fn can return R objects of different classes
* `lapply(nyc, class) #applies class() on the list elements of nyc`
  * returns a `list` containing the class corr to the input vector
* to convert this list into a vector use `unlist(lapply())`
* lapply with functions with more than one argument - `lapply(data, fn, args)`
* lapply with user defined function = `lapply(x, function(x){…})`

## sapply

* simplify apply, returns a named **vector**, containing the same info as the vector obtained with lapply; tries to simplify list to array
  * mechanism: applies lapply first that use simplify2array to convert list into array
* `sapply(cities, nchar, USE.NAMES = TRUE)` is default
  * `USE.NAMES = FALSE` to remove the names/original input
* if function returns more than one output - becomes row, the original input are the columns # but the length of o/p must be the same
  * or else, sapply will output similar but named list as lapply
  * using identical()
  * sapply() does not simplify list of NULL into vector bcus it will simply be NULL, which isnt a vector with the same length as input

## vapply

* apply fn over list or vector, explicitly specify output format
* `vapply(X, FUN, FUN.VALUE, …, USE.NAMES = TRUE)`
  * FUN = function
  * FUN.VALUE = the return data type, e.g. numeric(1) - nchar() returns a single numerical value
* prevents the case when sapply fails to convert the list into vector, returns warnings

## useful functions

* rep(seq, times = 1, each = 1), times = how many times the whole sequence relipcates, each = how many times each element in the sequence replicates
* sort(x, decreasing = FALSE) - default ascending order
* is._() , as._()
* unlist()
* str()
* append() - concatenates the vectors
* rev() - reverse

## regular expressions

* sequence of (meta)characters
* pattern existence
* pattern replacement
* pattern extraction
* `grep()`, `grepl()`
  * grepl(pattern = \<regex>, x = \<string>) - returns vector of logicals
  * ^ - starting with
  * $ - ending with
  * ?regex
  * grep() - returns the indices of the TRUE elements
  * which(grepl()) → grep()
* sub(), gsub()
  * sub(pattern = \<regex>, replacement = \<str>, x = \<str>)
    * sub only looks and replace the first one that match in string
    * to replace every occurrence, use gsub()
  * metacharacter: | = OR
*   metacharacter

    * `\\` to escape other metacharacters
    * `.*` - matches any character (.) zero or more times (\*)
    * `\\s`: Match a space. The "s" is normally a character, escaping it (`\\`) makes it a metacharacter.
    * `[0-9]+`: Match the numbers 0 to 9, at least once (+).
    * `([0-9]+)`: The parentheses are used to make parts of the matching string available to define the replacement. The `\\1` in the replacement argument of `sub()` gets set to the string that is captured by the regular expression `[0-9]+`.

    ```r
    awards <- c("Won 1 Oscar.",
      "Won 1 Oscar. Another 9 wins & 24 nominations.",
      "1 win and 2 nominations.",
      "2 wins & 3 nominations.",
      "Nominated for 2 Golden Globes. 1 more win & 2 nominations.",
      "4 wins & 1 nomination.")

    sub(".*\\s([0-9]+)\\snomination.*$", "\\1", awards)

    [1] "Won 1 Oscar." "24"           "2"            "3"            "2"           
    [6] "1"
    ```

## time and dates

* `Sys.Date()` - returns data in YYYY-MM-DD

<figure><img src="../.gitbook/assets/image (13).png" alt="" width="287"><figcaption></figcaption></figure>

* `as.Date()`

<figure><img src="../.gitbook/assets/image (14).png" alt="" width="563"><figcaption></figcaption></figure>

* `as.POSIXct`
* addition and subtraction of class-Date changes by days
* time unit for POSIXct object is seconds
* date object is a mroe advanced representation of the no of days since the first of 1 Jan in 1970
  * unclass(date) - returns the no of days from 1 Jan 1970
  * for POSIXct is the midnight of 1 Jan 1970
* dedicated libraries for dates
  * lubridate
  * zoo
  * xts
* symbols for date
  * `%Y`: 4-digit year (1982)
  * `%y`: 2-digit year (82)
  * `%m`: 2-digit month (01)
  * `%d`: 2-digit day of the month (13)
  * `%A`: weekday (Wednesday)
  * `%a`: abbreviated weekday (Wed)
  * `%B`: month (January)
  * `%b`: abbreviated month (Jan)

```r
# Convert dates to formatted strings
format(date1, "%A")

# Convert the strings to POSIXct objects: time1, time2
time1 <- as.POSIXct(str1, format = "%B %d, '%y hours:%H minutes:%M seconds:%S")
```

* symbols for POSIXct
  * `%H`: hours as a decimal number (00-23)
  * `%I`: hours as a decimal number (01-12)
  * `%M`: minutes as a decimal number
  * `%S`: seconds as a decimal number
  * `%T`: shorthand notation for the typical format `%H:%M:%S`
  * `%p`: AM/PM indicator

