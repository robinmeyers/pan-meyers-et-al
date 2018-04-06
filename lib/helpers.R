row_correlations <- function(x, ...) {
    aaply(x, 1, cor, y=t(x), ...) %>%
        set_colnames(rownames(x)) %>%
        set_rownames(rownames(x))
}

alphabetize_genes <- function(df, x=Gene.x, y=Gene.y) {
#
    x = enquo(x)
    y = enquo(y)

    new_x <- pmin(df %>% pull(!! x),
                  df %>% pull(!! y))
    new_y <- pmax(df %>% pull(!! x),
                  df %>% pull(!! y))

    df %>%
        mutate(!! quo_name(x) := new_x,
               !! quo_name(y) := new_y) %>%
        arrange(!! x, !! y)
}

de_excelify_genes <- function(df, x) {

    x <- enquo(x)

    df %>%
        mutate(!! quo_name(x) := ifelse(is.na(as.integer(!! x)), !! x,
                           as.Date(as.integer(!! x), origin="1900-01-01") %>%
                               format("%d-%b") %>% str_replace("-0", "-"))) %>%
        mutate(!! quo_name(x) := str_replace(!! x, "(\\d+)-Sep", "SEPT\\1"),
               !! quo_name(x) := str_replace(!! x, "(\\d+)-Mar", "MARCH\\1"))

}
