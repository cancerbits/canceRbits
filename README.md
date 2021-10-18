# R Package: canceRbits

Wrapper and helper functions commonly used by the Developmental 
Cancer Genomics group at St. Anna Children's Cancer Research Institute.

## To install from github

If a `GITHUB_PAT` environment variable is already set to your access token
```
remotes::install_github(repo = 'cancerbits/canceRbits')
```

Else, pass your token
```
remotes::install_github(repo = 'cancerbits/canceRbits', auth_token = 'ghp_xxx')
```

## Workflow for updating the package
* Add/edit one or more .R files in R directory
	* Make sure to include roxygen documentation
* run `devtools::document()`
* add required packages to NAMESPACE file (if there are new ones)
* (optional) Edit DESCRIPTION file (e.g. add your name, bump the version number)
* (optional) Edit overall package documentation in `R/canceRbits.R`
* (optional but strongly suggested) Add unit test for new functionality (in `tests/testthat/`)
* check whether package is still OK via `devtools::check(document = FALSE)` or via Rstudio menu: 'Build' > 'Check Package'

## Suggested reading

[R packages by Hadley Wickham and Jenny Bryan](https://r-pkgs.org/index.html)

[useful cheat sheet](https://rawgit.com/rstudio/cheatsheets/master/package-development.pdf)
