clr()

y1 = readxl::read_excel("statsforfigure2B-C.xlsx", sheet = 1)

get_geno = function(x) {
  stringr::str_split(x,"\\(")[[1]][[1]] %>% stringr::str_trim() %>% 
    stringr::str_extract("\\d\\s*\\w+$") %>% 
    stringr::str_replace_all("\\d+","") %>% stringr::str_trim()
}

y = y1[1:5,6:17] %>% dplyr::mutate_if(is.character, as.numeric) %>% 
  dplyr::bind_cols(y1[1:5,5]) %>% 
  dplyr::rename(organelle = `Total adj.SPC`) %>%
  tidyr::gather("des", "val", -organelle) %>%
  dplyr::mutate(gen = purrr::map_chr(des, get_geno)) %>%
  dplyr::mutate(rep = stringr::str_extract(des,"rep\\d"))

y %>% dplyr::count(gen, rep) %>% print.data.frame()

# y %>% print.data.frame()

# y %>% dplyr::group_by(organelle,gen) %>% 
#   dplyr::summarise(avg = mean(val), std = sd(val)) %>% 
#   dplyr::ungroup() %>% print.data.frame()

# # for initial results in Feb 2018 (each against WT)
# all_combs = tibble(org = unique(y$organelle), gen1 = "WT") %>%
#   tidyr::crossing(gen2 = purrr::keep(unique(y$gen), ~ .x != "WT"))

# for new results in June 2018 (every pair)
all_combs = utils::combn(unique(y$gen),2) %>% t() %>% tibble::as_tibble() %>% 
  tidyr::crossing(org = unique(y$organelle)) %>% dplyr::select(org, gen1 = V1, gen2 = V2) %>% 
  dplyr::arrange(org, gen1, gen2)

Y = all_combs %>% 
  dplyr::mutate(
    dat = purrr::pmap(
      list(org, gen1, gen2), 
      function(oe,g1,g2) {
        y %>% dplyr::filter(organelle == oe, gen %in% c(g1,g2)) %>% 
          dplyr::select(rep,gen,val) # %>% tidyr::spread(key = gen, value = val)
      }
    )
  ) %>% 
  dplyr::mutate(
    ttN = purrr::map(dat, function(x) t.test(val ~ gen, data = x, var.equal = FALSE, paired = FALSE, alternative = "two.sided")),
    vt = purrr::map(dat, function(x) var.test(val ~ gen, data = x)),
    ttE = purrr::map(dat, function(x) t.test(val ~ gen, data = x, var.equal = TRUE, paired = FALSE, alternative = "two.sided"))
  ) %>% 
  dplyr::mutate(
    ttNP = purrr::map(ttN, function(x) x %>% broom::tidy() %>% dplyr::transmute(statistic,`p.value`,inference = ifelse(`p.value` < 0.05, "NOTEQ", "EQUAL"),method,alternative)),
    vtP = purrr::map(vt, function(x) x %>% broom::tidy() %>% dplyr::transmute(statistic,`p.value`,inference = ifelse(`p.value` < 0.05, "NOTEQ", "EQUAL"),method,alternative)),
    ttEP = purrr::map(ttE, function(x) x %>% broom::tidy() %>% dplyr::transmute(statistic,`p.value`,inference = ifelse(`p.value` < 0.05, "NOTEQ", "EQUAL"),method,alternative))
  )

dfs = list(
  Y %>% dplyr::select(org,gen1,gen2,ttNP) %>% tidyr::unnest(ttNP),
  Y %>% dplyr::select(org,gen1,gen2,vtP) %>% tidyr::unnest(vtP),
  Y %>% dplyr::select(org,gen1,gen2,ttEP) %>% tidyr::unnest(ttEP)
)

dfs %<>% purrr::map(function(x) {
  colnames(x) = toupper(colnames(x))
  x
})

toxl(dfs, sheet_names = c("ttest_uneqVaR","var_test","ttest_eqVar"), open_file = FALSE, out_file = "Fig2B_tests.xlsx")

# Y %>% dplyr::select(org,gen1,gen2,ttNP) %>% tidyr::unnest(ttNP) %>% write_csv("Fig2B_ttest_uneqVar.csv")
# Y %>% dplyr::select(org,gen1,gen2,vtP) %>% tidyr::unnest(vtP) %>% write_csv("Fig2B_vartest.csv")
# Y %>% dplyr::select(org,gen1,gen2,ttEP) %>% tidyr::unnest(ttEP) %>% write_csv("Fig2B_ttest_eqVar.csv")
