### Info

This package is not useful in it's own. It was created to test and validate the package
`sfreemap`.


### Package Requirements for Production

You need to have R installed on your system. If you are using a debian/ubuntu based distribution, just type the following command in a terminal.

`sudo apt-get install r-base-core`

### Package Requirements for Development

`sudo apt-get install r-base-core texlive-full`

### Install

```
install.packages('devtools')
install_github('dpasqualin/sfreemap.tests')
```

If you have troubles installing the `devtools` package, try downloading
`sfreemap.tests` and then building and installing it using the following commands:

```
git clone https://github.com/dpasqualin/sfreemap.tests.git
R CMD build sfreemap.tests && R CMD INSTALL sfreemap.tests
```

If you choose to install using the command above, the documentation will be
available in the directory `sfreemap.tests.Rcheck`.

Some people might have problems with package `Briostrings` as well, which is
a dependency of `phangorn`, which is a dependency of `sfreemap`. If you do,
the official Biostrings website suggest the following commands to install
it:

```
## try http:// if https:// URLs are not supported
source('https://bioconductor.org/biocLite.R')
biocLite('Biostrings')
```
