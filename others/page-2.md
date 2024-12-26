---
description: Setting up reproducible environments
---

# renv

## Setting up renv

<figure><img src="../.gitbook/assets/image.png" alt="" width="459"><figcaption></figcaption></figure>

Requires you to set up a new Rproject for each project.

```r
renv::init()
renv::snapshot() # updates the lockfile with metadata about the currently-used packages in the project library
renv::restore() # can easily reproduce your current environment by running

renv::install("Seurat") # installs r packages from a variety of sources
```

*   To convert a project to use renv, call [`renv::init()`](https://rstudio.github.io/renv/reference/init.html). It adds three new files and directories to your project:

    * project library, `renv/library` that contains all packages currently used by your project
    * **lockfile**, `renv.lock`, records enough metadata about every package that it can be re-installed on a new machine
    * project R profile, `.Rprofile` , run automatically every time you start R (in that project), and renv uses it to configure your R session to use the project library. This ensures that once you turn on renv for a project, it stays on, until you deliberately turn it off.


* Convenient installing with `renv::install(...)`
  * `pkg`: install latest version of `pkg` from CRAN.
  * `pkg@version`: install specified version of `pkg` from CRAN.
  * `username/repo`: install package from GitHub
  * `bioc::pkg`: install `pkg` from Bioconductor.

### Session info

```r
renv::snapshot()
sessionInfo()
```

## Updating packages <a href="#updating-packages" id="updating-packages"></a>

We recommend that you regularly (at least once a year) or if you’re making major changes to a project that you haven’t worked on for a while use [`renv::update()`](https://rstudio.github.io/renv/reference/update.html) to get the latest versions of all dependencies and renv.&#x20;

* After calling [`renv::update()`](https://rstudio.github.io/renv/reference/update.html), you should run the code in your project and verify that it still works (or make any changes needed to get it working). Then call [`renv::snapshot()`](https://rstudio.github.io/renv/reference/snapshot.html) to record the new versions in the lockfile.&#x20;
* If you get stuck, and can’t get the project to work with the new versions, you can call [`renv::restore()`](https://rstudio.github.io/renv/reference/restore.html) to roll back changes to the project library and revert to the known good state recorded in your lockfile. If you need to roll back to an even older version, take a look at [`renv::history()`](https://rstudio.github.io/renv/reference/history.html) and [`renv::revert()`](https://rstudio.github.io/renv/reference/history.html).

See [`renv::upgrade()`](https://rstudio.github.io/renv/reference/upgrade.html) if you ever want to upgrade just renv, or you need to install a development version from GitHub.

## Uninstalling renv <a href="#uninstalling-renv" id="uninstalling-renv"></a>

* To deactivate renv in a project, use [`renv::deactivate()`](https://rstudio.github.io/renv/reference/activate.html).&#x20;
  * This removes the renv auto-loader from the project `.Rprofile`, but doesn’t touch any other renv files used in the project. If you’d like to later re-activate renv, you can do so with [`renv::activate()`](https://rstudio.github.io/renv/reference/activate.html).
* To completely remove renv from a project, call `renv::deactivate(clean = TRUE)`.&#x20;
  * If you later want to use renv for this project, you’ll need to start from scratch with `renv::init().`

If you want to stop using renv for all your projects, you’ll also want to remove `renv`'s global infrastructure with the following R code:

```r
root <- renv::$root()
(root, recursive = TRUE)
```

* To uninstall the renv package with `utils::remove.packages("renv")`.

