# Toolbox.R

require(reshape2)

is_integer = function(x) {
    if (all(is.na(as.numeric(x)))) {
        return(FALSE)
    } else {
        unique(as.numeric(x)) %>% .[!is.na(.)] %>% {. %% 1 == 0} %>% all
    }
}

is_continuous = function(x) {
    is.numeric(x) & (!is_integer(x))
}

# check if vector of data is a binary indicator
is_indicator = function(x) {

    x = x[!is.na(x)]
    
    values = x %>% as.character() %>% unique
    
    return(
        setequal(values, c('0', '1')) | setequal(values, c('0')) | setequal(values, c('1'))
    )
}

crosstab = function(D, x, y, fraction = FALSE, overall = TRUE) {
        
    ylevels = unique(D[[y]]) %>% .[!is.na(.)] %>% as.character
    
    if (tolower(x) == 'total') {

        D = D %>% dcast(
                paste0('. ~ ', y), 
                fun.aggregate = length,
                value.var = y
            ) %>%
            dplyr::rename('levels' = '.') %>%
            mutate(
                levels = 'total',
                variable = 'total',
                overall = 'total'
            ) %>%
            select(all_of(c('variable', 'levels', ylevels)))

        if (overall) {
            D = D %>% mutate(overall = sum(.[ylevels]))
        }

        return(D)

    } else if (tolower(y) == 'total') {
        D %>% count(get(x)) %>% 
        dplyr::rename('levels' = 'get(x)', 'total' = 'n') %>%
        mutate(variable = gsub('_', '-', x)) %>%
        select(all_of(c('variable', 'levels', 'total')))
    } else if (is_continuous(D[[y]]) & (is_continuous(D[[x]]) | class(D[[x]]) == 'integer')) {
        NA
    } else if ((is_continuous(D[[x]]) | class(D[[x]]) == 'integer') & !is_indicator(D[[x]])) {

        D[['y']] = D[[y]]

        if (overall) {

            D = rbind(
                D,
                mutate(D, y = 'overall')
            )

            ylevels = c(ylevels, 'overall')
        }

        # x is continuous
        D %>% 
        filter(!is.na(get(x))) %>%
        dcast(
            paste0('. ~ ', 'y'), 
            fun.aggregate = function(x){
                paste0(round(mean(x), 1), ' (', round(sd(x), 1), ')')
            },
            value.var = x
        ) %>%
        dplyr::rename('levels' = '.') %>%
        mutate(levels = 'mean (SD)') %>%
        mutate(variable = gsub('_', '-', x)) %>%
        select(all_of(c('variable', 'levels', ylevels)))

    } else if (is_continuous(D[[y]])) {
        # y is continous
        D %>% dcast(
            paste0('. ~ ', x), 
            fun.aggregate = function(x){
                paste0(round(mean(x), 1), ' (', round(sd(x), 1), ')')
            },
            value.var = y
        ) %>% .[,-1] %>% t %>% 
        as.data.frame %>%
        tibble::rownames_to_column(var = 'levels') %>%
        setNames(c('levels', 'mean (SD)')) %>%
        mutate(variable = gsub('_', '-', x)) %>%
        select(all_of(c('variable', 'levels', 'mean (SD)')))
    } else {
        res = table(D[[x]], D[[y]]) %>%
            as.data.frame.matrix %>% as.data.frame %>%
            tibble::rownames_to_column(var = 'levels') %>%
            mutate(variable = gsub('_', '-', x)) %>%
            select(all_of(c('variable', 'levels', ylevels))) %>%
            mutate(overall = sum(.[ylevels]))

        if (fraction) {
            res = res %>%
                group_by(variable) %>%
                mutate_at(ylevels, function(x){paste0(x, ' (', round(x * 100/sum(x), 1), '%', ')')}) %>%
                ungroup()
        }

        return(res)
    }
}

tab_response = function(D, xs, y, fraction = FALSE) {    
    lapply(
        xs,
        function(x) {crosstab(D, x, y, fraction)}
    ) %>% 
    Reduce(rbind, .)
}