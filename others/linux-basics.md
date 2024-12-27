# linux basics

## checking size

* `du -sh <directory name>`



## awk

[https://www.geeksforgeeks.org/awk-command-unixlinux-examples/](https://www.geeksforgeeks.org/awk-command-unixlinux-examples/)

* awk usage
* -F - indicate separator
* ‘{printf …}’



## storing outputs

* \>file.txt



## copying files

`cp ./DirectoryA_1/README.txt ./DirectoryA_2`

`cp ./DirectoryA_1/*.txt ./DirectoryA_2`



## sed

This works:

```
sed -rne 's/(dbservername)\\s+\\w+/\\1 yyy/gip'

```

(When you use the -r option, you don't have to escape the parens.)

Bit of explanation:

* <mark style="color:red;">`r`</mark> is extended regular expressions - makes a difference to how the regex is written.
* <mark style="color:red;">`n`</mark> does not print unless specified - `sed` prints by default otherwise,
* <mark style="color:red;">`e`</mark> means what follows it is an expression. Let's break the expression down:
  * <mark style="color:red;">`s///`</mark> is the command for search-replace, and what's between the first pair is the regex to match, and the second pair the replacement,
  * <mark style="color:red;">`gip`</mark>, which follows the search replace command; <mark style="color:red;">`g`</mark> means global, i.e., every match instead of just the first will be replaced in a line; <mark style="color:red;">`i`</mark> is case-insensitivity; <mark style="color:red;">`p`</mark> means print when done (remember the `n` flag from earlier!),
  * The brackets represent a match part, which will come up later. So <mark style="color:red;">`dbservername`</mark> is the first match part,
  * <mark style="color:red;">`\\s`</mark> is whitespace, `+` means one or more (vs , zero or more) occurrences,
  * <mark style="color:red;">`\\w`</mark> is a word, that is any letter, digit or underscore,
  * <mark style="color:red;">`\\1`</mark> is a special expression for GNU `sed` that prints the first bracketed match in the accompanying search.

sed is the right tool, try doing :

```bash
var="movie.MOV"
sed -i "1s/.*/$var/" file.txt

```

**explanations**

* <mark style="color:red;">`1`</mark> mean first line
* the rest is the substitution <mark style="color:red;">`s///`</mark> : we substitute everything (.\*) by the <mark style="color:red;">`$var`</mark> variable
* the double shell quotation is mandatory here
  * double quotes are used when substituting variables; single quotes substitutes the literal term

## giving permissions

[https://www.pluralsight.com/blog/it-ops/linux-file-permissions](https://www.pluralsight.com/blog/it-ops/linux-file-permissions)



## checking md5

