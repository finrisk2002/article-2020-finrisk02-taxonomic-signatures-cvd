# nice tool script from Matti R. 2019
# pkg is a list of package names
# check if installed, if not install with BiocManager
# afterward, require all packages in pkg
is.installed <- function(pkg, myvrs="3.8") {
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg)) {
    BiocManager::install(new.pkg, ask=F, version=myvrs)
  }
  sapply(pkg, require, character.only = TRUE)
}