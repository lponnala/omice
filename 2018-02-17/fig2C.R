clr()

y2 = readxl::read_excel("statsforfigure2B-C.xlsx", sheet = 2)

y = y2 %>% dplyr::rename(Function = X__1) %>% 
  dplyr::select(Function,contains("Aver")) %>%
  tidyr::gather(key = "des", value = "val", -Function) %>%
  dplyr::mutate(gen = purrr::map_chr(des, function(x) tolower(x) %>% stringr::str_split("average") %>% .[[1]] %>% .[[2]] %>% stringr::str_trim()))

# # for initial results in Feb 2018 (each against WT)
# all_combs = tibble(gen1 = "wt") %>% 
#   tidyr::crossing(gen2 = unique(y$gen) %>% purrr::keep(~ .x != "wt"))

# for new results in June 2018 (every pair)
all_combs = utils::combn(unique(y$gen),2) %>% t() %>% tibble::as_tibble() %>% 
  dplyr::rename(gen1 = V1, gen2 = V2)

Y = all_combs %>%
  dplyr::mutate(dat = purrr::map2(gen1, gen2, function(g1,g2) 
    y %>% dplyr::filter(gen %in% c(g1,g2))
  )) %>%
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
  Y %>% dplyr::select(gen1,gen2,ttNP) %>% tidyr::unnest(ttNP),
  Y %>% dplyr::select(gen1,gen2,vtP) %>% tidyr::unnest(vtP),
  Y %>% dplyr::select(gen1,gen2,ttEP) %>% tidyr::unnest(ttEP)
)

dfs %<>% purrr::map(function(x) {
  colnames(x) = toupper(colnames(x))
  x
})

toxl(dfs, sheet_names = c("ttest_uneqVaR","var_test","ttest_eqVar"), open_file = FALSE, out_file = "Fig2C_tests.xlsx")

# Y %>% dplyr::select(gen1,gen2,ttNP) %>% tidyr::unnest(ttNP) %>% write_csv("Fig2C_ttest_uneqVar.csv")
# Y %>% dplyr::select(gen1,gen2,vtP) %>% tidyr::unnest(vtP) %>% write_csv("Fig2C_vartest.csv")
# Y %>% dplyr::select(gen1,gen2,ttEP) %>% tidyr::unnest(ttEP) %>% write_csv("Fig2C_ttest_eqVar.csv")
