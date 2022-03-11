require(fs)
require(stringr)

ordered.ls.dir <- function (my.dir, regexp="")
{
    return (dir_ls(my.dir, regexp=regexp)[order(
                      strtoi(
                          str_match(dir_ls(my.dir, regexp=regexp), "_([0-9]+)[_\\.$]")[,2])
                  )]
            )
}
