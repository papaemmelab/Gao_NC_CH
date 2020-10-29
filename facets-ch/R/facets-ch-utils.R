suppressMessages(library("argparse", quietly = TRUE))
suppressMessages(library("dplyr", quietly = TRUE))
suppressMessages(library("data.table", quietly = TRUE))
suppressMessages(library("ggplot2", quietly = TRUE))
suppressMessages(library("patchwork", quietly = TRUE))
suppressMessages(library('stringr', quietly = TRUE))
suppressMessages(library('glue', quietly = TRUE))

get_logr = function(pu1, pu2, snp.nbhd, cval, snp_blacklist = c(), het_whitelist = c(), het_bias = NULL, gamma = 0.25, blacklist = TRUE) {
    
    rcmat = inner_join(
        pu1,
        pu2 %>% rename(File2R = File1R, File2A = File1A, File2E = File1E, File2D = File1D),
        by = c("Chromosome", "Position", "Ref", "Alt")
    ) %>%
    mutate(
        NOR.DP = File2R + File2A,
        NOR.RD = File2R,
        TUM.DP = File1R + File1A,
        TUM.RD = File1R
    ) %>%
    mutate(
        baf = File1A / (File1R + File1A),
        het = baf < 0.9 & baf > 0.1
    )
        
    if (blacklist) {
        rcmat = rcmat %>%
            mutate(
                Chromosome = ifelse(Chromosome == 'X', '23', Chromosome),
                Chromosome = ifelse(Chromosome == 'Y', '24', Chromosome)
            ) %>%
            mutate(snp_index = paste0(Chromosome, '_', Position)) %>%
            filter(((!het) & (!snp_index %in% snp_blacklist)) | ((het) & (snp_index %in% het_whitelist))) %>%
            filter(!(het & Chromosome == '6' & Position < 58400000 & snp_index %in% snp_blacklist)) %>%
            filter(!(het & Chromosome == '7' & Position > 61100000 & snp_index %in% snp_blacklist)) %>%
            filter(!(het & Chromosome == '2' & Position > 50e6 & snp_index %in% snp_blacklist)) %>%
#             filter(!(het & Chromosome == '3' & Position > 93200000 & snp_index %in% snp_blacklist)) %>%
#             filter(!(het & Chromosome == '1' & Position > 121100000 & snp_index %in% snp_blacklist)) %>%
            mutate(
                Chromosome = ifelse(Chromosome == '23', 'X', Chromosome),
                Chromosome = ifelse(Chromosome == '24', 'Y', Chromosome)
            )
    }
        
    set.seed(0) # for reproducibility
    # https://github.com/mskcc/facets/issues/76
    xx <- preProcSample(rcmat, ndepth = 50, het.thresh = 0.1, ndepthmax = 1500,
                        unmatched = TRUE, snp.nbhd = snp.nbhd, cval = 25, het_bias = het_bias,
                        gamma = gamma)
    
    oo <- procSample(xx, cval = cval)

    return(oo$jointseg)
}

phis = seq(0, 0.9, by = 0.01)

cnlr = function(phi, m = 0, p = 1) {
  log2((2 * (1 - phi) + (m + p) * phi) / 2)
}

valor = function(phi, m = 0, p = 1) {
  abs(log((m * phi + 1 - phi) / (p * phi + 1 - phi)))
}

del_line = data.frame(
  phi = phis,
  cnlr = sapply(
    phis,
    function(phi){cnlr(phi, m = 0, p = 1)}
  ),
  valor = sapply(
    phis,
    function(phi){valor(phi, m = 0, p = 1)}
  ))

amp_line = data.frame(
  phi = phis,
  cnlr = sapply(
    phis,
    function(phi){cnlr(phi, m = 2, p = 1)}
  ),
  valor = sapply(
    phis,
    function(phi){valor(phi, m = 2, p = 1)}
  ))

loh_line = data.frame(
  phi = phis,
  cnlr = sapply(
    phis,
    function(phi){cnlr(phi, m = 2, p = 0)}
  ),
  valor = sapply(
    phis,
    function(phi){valor(phi, m = 2, p = 0)}
  ))

variant_type = function(cnlr_seg, valor_seg, aberrant = TRUE) {
  
  if ((!aberrant) | is.na(valor_seg)) {
    return(res = c('type' = 'none', 'phi' = NA, 'err' = NA))
  }
  
  del_opt = del_line %>% mutate(err = (cnlr - cnlr_seg)^2 + (valor - valor_seg)^2) %>% arrange(err) %>% .[1,c('err', 'phi')]
  amp_opt = amp_line %>% mutate(err = (cnlr - cnlr_seg)^2 + (valor - valor_seg)^2) %>% arrange(err) %>% .[1,c('err', 'phi')]
  loh_opt = loh_line %>% mutate(err = (cnlr - cnlr_seg)^2 + (valor - valor_seg)^2) %>% arrange(err) %>% .[1,c('err', 'phi')]
      
  min_err = min(c(del_opt[['err']], amp_opt[['err']], loh_opt[['err']]))
  
  if (min_err > 0.1) {
    res = c('type' = 'err', 'phi' = NA, 'err' = min_err)
  } else if(min_err == del_opt[['err']]) {
    res = c('type' = 'del', 'phi' = del_opt[['phi']], 'err' = min_err)
  } else if (min_err == amp_opt[['err']]) {
    res = c('type' = 'amp', 'phi' = amp_opt[['phi']], 'err' = min_err)
  } else {
    res = c('type' = 'loh', 'phi' = loh_opt[['phi']], 'err' = min_err)
  }
  
  res[['phi']] = as.numeric(res[['phi']])
  res[['err']] = as.numeric(res[['err']])
  
  return(res)
}

calc_sem_x = function(sig, seg_id, n_marker) {
        
    cnlrs = sig %>% filter((!aberrant) & (seg != seg_id)) %>% pull(cnlr)
    
    sd(
        sapply(
            1:(length(cnlrs) - n_marker),
            function(i) {
                mean(cnlrs[i:(i+n_marker)])
        })
    )
}

clac_x_bar = function(x, aberrant, n_marker) {
    unlist(lapply(1:length(x), function(i) {mean(x[-i][(!aberrant[-i]) & (n_marker[-i] > 30)])}))
}

call_aberrant = function(sig, seg, p_threshold = 1e-8, k = 0, debug = FALSE, trace = data.frame()) {
    
    num_aberrant = sum(as.integer(seg$aberrant))
    
    seg = seg %>%
        select(-one_of('sd_y', 'y_bar', 'sem_y')) %>%
        left_join(
            sapply(unique(sig$seg), function(s){
                sig_prime = sig %>% filter(seg != s) %>% filter(!aberrant)
                sd_y = sd(na.omit(abs(sig_prime$valor)))
                y_bar = mean(na.omit(abs(sig_prime$valor)))
                return(c('seg' = s, 'sd_y' = sd_y, 'y_bar' = y_bar))
            }) %>% t %>%
            as.data.frame(),
            by = 'seg'
        ) %>%
        mutate(
            x_bar = clac_x_bar(cnlr, aberrant, n_marker)
        ) %>%
        rowwise() %>%
        mutate(
            sem_x = calc_sem_x(sig, seg, n_marker),
            sem_y = sd_y / sqrt(n_het),
            z_x = (cnlr - x_bar) / sem_x,
            p_x = min(1 - pnorm(z_x), pnorm(z_x)),
            z_y = (abs(valor) - y_bar) / sem_y,
            p_y = ifelse(is.na(z_y), 1, 1 - pnorm(z_y)),
            chisq = ifelse(is.na(z_y), z_x^2, z_y^2 + z_x^2),
            p_chisq = 1 - pchisq(chisq, df = 2),
            aberrant = p_chisq < p_threshold & p_y < 0.05
        ) %>%
        ungroup() %>%
        mutate(
            diplogr = median(cnlr[!aberrant]), 
            cnlr_adj = cnlr - diplogr
        )
        
    sig = sig %>%
        select(-one_of('aberrant')) %>%
        left_join(
            seg %>% select(seg, chrom, aberrant),
            by = c('seg', 'chrom')
        )
    
    num_aberrant_new = sum(as.integer(seg$aberrant))
    
    if (debug) {
        print(glue('Iteration {k}. Found {num_aberrant_new} aberrant segments.'))
        trace = rbind(
            trace, 
            seg %>% mutate(k = k)
        )
    }
                
    if (num_aberrant_new == num_aberrant) {
        if (debug) {
            return(list('seg' = seg, 'trace' = trace))
        } else {
            return(seg)
        }
    } else {
        return(call_aberrant(sig, seg, p_threshold, k = k + 1, debug = debug, trace = trace))
    }
}

get_seg = function(sig) {
    
    seg = sig %>% group_by(seg, chrom) %>%
        summarise(
            start_snp = min(snp_i),
            end_snp = max(snp_i),
            cnlr = mean(na.omit(cnlr)),
            valor = mean(abs(na.omit(valor))),
            n_marker = n(),
            n_het = sum(het),
            start = min(maploc),
            end = max(maploc),
            size = max(maploc) - min(maploc)
        ) %>%
        ungroup() %>%
        mutate(aberrant = FALSE)
    
    return(seg)
}

facets_plot = function(sig, seg, title = '', snp_distance = TRUE, return_list = FALSE) {
  
    if (snp_distance) {
        sig = sig %>% mutate(pos = snp_i)
        seg = seg %>% mutate(start = start_snp, end = end_snp)
    } else {
        sig = sig %>% mutate(pos = maploc)
    }

    panel_theme = theme_classic() +
      theme(
          axis.text.x = element_blank(),
          panel.spacing = unit(0, 'mm'),
          axis.ticks.x = element_blank(),
          strip.background = element_blank(),
          panel.border = element_rect(color = 'gray', size = 0.5, fill = NA),
          axis.line = element_blank(),
          legend.position = 'none',
          plot.title = element_text(hjust = 0.5, size = 10),
          plot.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt")
      )

    p_cnlr = ggplot(
            sig,
            aes(x = pos, y = cnlr)
        ) +
        geom_vline(
            data = sig %>% group_by(chrom, arm) %>% 
                summarise(start = min(pos), end = max(pos)) %>%
                group_by(chrom) %>%
                filter(n() > 1) %>%
                filter(arm == 'q'),
            aes(xintercept = start),
            color = 'royalblue',
            linetype = 'dotted',
            size = 0.3,
            alpha = 1
        ) +
        geom_point(
            shape = 21, size = 0.4, color = 'gray60'
        ) +
        geom_segment(
            data = seg,
            aes(x = start, xend = end, y = cnlr, yend = cnlr, color = aberrant),
            size = 0.75,  alpha = 0.9
        ) + 
        facet_grid(.~chrom, scale = 'free', space = 'free_x') +
        panel_theme +
        ylim(-1,1) +
        ylab('cnlr') +
        xlab('') +
        scale_color_manual(values = c('black', 'red')) +
        ggtitle(title)

    p_valor = ggplot(
          sig,
          aes(x = pos, y = valor)
        ) +
        geom_vline(
            data = sig %>% group_by(chrom, arm) %>% 
                summarise(start = min(pos), end = max(pos)) %>%
                group_by(chrom) %>%
                filter(n() > 1) %>%
                filter(arm == 'q'),
            aes(xintercept = start),
            color = 'royalblue',
            linetype = 'dotted',
            size = 0.3,
            alpha = 1
        ) +
        geom_point(
          shape = 21, size = 0.4, color = 'gray60'
        ) +
        geom_point(
          data = sig %>% mutate(valor = -valor),
          shape = 21, size = 0.4, color = 'gray60'
        ) +
        geom_segment(
          data = seg %>% mutate(value = valor, variable = 'valor'),
          aes(x = start, xend = end, y = value, yend = value, color = aberrant),
          size = 0.75, alpha = 0.9
        ) +
        geom_segment(
          data = seg %>% mutate(value = -valor, variable = 'valor'),
          aes(x = start, xend = end, y = value, yend = value, color = aberrant),
          size = 0.75, alpha = 0.9
        ) +
        facet_grid(.~chrom, scale = 'free', space = 'free_x') +
        ylim(-2.25, 2.25) +
        ylab('valor') +
        xlab('') +
        scale_color_manual(values = c('black', 'red')) +
        panel_theme +
        theme(
            strip.text = element_blank()
        )
    if (return_list) {
        return(list('p_cnlr' = p_cnlr, 'p_valor' = p_valor))
    } else {
        return(p_cnlr/p_valor)
    }
}