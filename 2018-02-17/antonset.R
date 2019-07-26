clr()

getgen = function(x) {
  x %>% stringr::str_split("\\s") %>% .[[1]] %>% 
    purrr::keep(~ .x != "") %>% .[[2]]
}

y = readxl::read_excel("fort-test-Antonset.xlsx") %>%
  dplyr::slice(1:28) %>% dplyr::select(2:9) %>%
  dplyr::rename(Function = `simple functional name (July 2017)`) %>%
  tidyr::gather(key = "des", value = "val", -Function, -location) %>%
  dplyr::mutate(rep = stringr::str_extract(des,"rep \\d")) %>%
  dplyr::mutate(gen = purrr::map_chr(des, ~ .x %>% getgen()))

y %>% dplyr::count(gen)
y %>% dplyr::count(rep)

Y = y %>% dplyr::group_by(Function, location) %>% tidyr::nest(.key = "dat") %>% 
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
  Y %>% dplyr::select(Function,location,ttNP) %>% tidyr::unnest(ttNP),
  Y %>% dplyr::select(Function,location,vtP) %>% tidyr::unnest(vtP),
  Y %>% dplyr::select(Function,location,ttEP) %>% tidyr::unnest(ttEP)
)

dfs %<>% purrr::map(function(x) {
  colnames(x) = toupper(colnames(x))
  x
})

toxl(dfs, sheet_names = c("ttest_uneqVaR","var_test","ttest_eqVar"), out_file = "Antonset_tests.xlsx")

# Y %>% dplyr::select(Function,location,ttNP) %>% tidyr::unnest(ttNP) %>% write_csv("Antonset_ttest_uneqVar.csv")
# Y %>% dplyr::select(Function,location,vtP) %>% tidyr::unnest(vtP) %>% write_csv("Antonset_vartest.csv")
# Y %>% dplyr::select(Function,location,ttEP) %>% tidyr::unnest(ttEP) %>% write_csv("Antonset_ttest_eqVar.csv")
