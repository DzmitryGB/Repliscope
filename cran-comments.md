## Test environments
* ubuntu 18.04.1, R 3.4.4
* windows 10, R 3.5.2
* MacOS 10.12.6, R 3.3.2
## R CMD check results
There were no ERRORs or WARNINGs or NOTEs. 
3

## Downstream dependencies
I have also run R CMD check on downstream dependencies using devtools::revdep_check(). I have attempted to use revdepcheck package but it seems to only work when there is a pre-existing CRAN package available, and, therefore, unsuitable for the first submission.

I have received the following error:
* devtools: Error: 'sessioninfo' is not in Suggests: for 'devtools'
 The issue has been reported already over two weeks ago:
 https://github.com/r-lib/devtools/issues/1951

Further checks (using skip="devtools", after reset) do not execute due to the same error.
