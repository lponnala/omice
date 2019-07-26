
# see "notes.txt" for thoughts and explanation
clr()
load("funcs-data.RData")
y %<>% 
  tidyr::gather(key = "gr", value = "val", -func) %>% 
  tidyr::separate(col = "gr", into = c("gen","rep"))
y %>% dplyr::count(func, gen)
y %>% dplyr::count(gen)
y %>% dplyr::count(rep)


# // 2-way anova for testing differences in genotype+function //
y %>% dplyr::group_by(gen,func) %>% dplyr::summarise(m = mean(val)) %>% dplyr::ungroup() %>% dplyr::summarise(cv = sd(m)/mean(m))
# .. do anova ..
ya = y %>% dplyr::mutate(func = as.factor(func), gen = as.factor(gen))
with_interaction = TRUE
if (with_interaction) {
  aov2 = stats::aov(val ~ func*gen, data = ya) # in genotype, function, and their interaction
  summary(aov2) %>% print() # func, gen, interaction are all significant
  aov2 %>% broom::tidy() %>% print()
  out_file = "two-way-anova-with-interaction.xlsx"
} else {
  aov2 = stats::aov(val ~ func + gen, data = ya) # in genotype, function
  summary(aov2) %>% print() # func-effect is very significant, gen-effect too but less so
  aov2 %>% broom::tidy() %>% print()
  out_file = "two-way-anova.xlsx"
}
# .. which pairs of genotype are different? ..
stats::TukeyHSD(aov2, which = "gen") %>% broom::tidy() %>% tibble::as_tibble() %>% dplyr::arrange(`adj.p.value`)
# stats::TukeyHSD(aov2, which = "func") %>% broom::tidy() %>% tibble::as_tibble() %>% dplyr::arrange(`adj.p.value`)
# summary(multcomp::glht(aov2, linfct = multcomp::mcp(gen = "Tukey"))) %>% broom::tidy() %>% tibble::as_tibble() %>% dplyr::arrange(`p.value`)
# stats::pairwise.t.test(ya$val, ya$gen, p.adjust.method = "BH") %>% broom::tidy() %>% tibble::as_tibble() %>% dplyr::arrange(`p.value`)
# .. test validity of ANOVA assumptions ..
# homoegeneity of variance
plot(aov2, 1)
car::leveneTest(val ~ func*gen, data = ya) # rejects null-hyp, so equal-variance is *not* a valid assumption!
# normality
plot(aov2, 2)
shapiro.test(x = residuals(object = aov2)) # rejects assumption of normality
# .. write output ..
list(
  aov2 %>% broom::tidy(), 
  stats::TukeyHSD(aov2, which = "gen") %>% broom::tidy() %>% dplyr::arrange(`adj.p.value`)
) %>% toxl(sheet_names = c("2-way ANOVA","Compare genotypes"), open_file = FALSE, out_file = out_file)


# // 1-way anova for testing differences in genotype //
y %>% dplyr::group_by(gen) %>% dplyr::summarise(m = mean(val)) %>% dplyr::summarise(cv = sd(m)/mean(m))
aov1 = stats::aov(val ~ gen, data = y)
summary(aov1) # no significant "genotype" effect: probably gets muddied due to all the varying levels of "function"
aov1 %>% broom::tidy()


# // t-test of pairwise genotype difference in each function //
all_combs = utils::combn(unique(y$gen),2) %>% t() %>% tibble::as_tibble() %>%
  tidyr::crossing(func = unique(y$func)) %>% dplyr::select(func, gen1 = V1, gen2 = V2) %>%
  dplyr::arrange(func, desc(gen1), gen2)
Y = all_combs %>%
  dplyr::mutate(data = purrr::pmap(list(func,gen1,gen2), function(fn,g1,g2) y %>% dplyr::filter(func == fn, gen %in% c(g1,g2)))) %>%
  dplyr::mutate(
    ttN = purrr::map(data, function(x) t.test(val ~ gen, data = x, var.equal = FALSE, paired = FALSE, alternative = "two.sided")),
    vt = purrr::map(data, function(x) var.test(val ~ gen, data = x)),
    ttE = purrr::map(data, function(x) t.test(val ~ gen, data = x, var.equal = TRUE, paired = FALSE, alternative = "two.sided"))
  ) %>%
  dplyr::mutate(
    ttNP = purrr::map(ttN, function(x) x %>% broom::tidy() %>% dplyr::transmute(statistic,`p.value`,inference = ifelse(`p.value` < 0.05, "NOTEQ", "EQUAL"),method,alternative)),
    vtP = purrr::map(vt, function(x) x %>% broom::tidy() %>% dplyr::transmute(statistic,`p.value`,inference = ifelse(`p.value` < 0.05, "NOTEQ", "EQUAL"),method,alternative)),
    ttEP = purrr::map(ttE, function(x) x %>% broom::tidy() %>% dplyr::transmute(statistic,`p.value`,inference = ifelse(`p.value` < 0.05, "NOTEQ", "EQUAL"),method,alternative))
  )
dfs = list(
  Y %>% dplyr::select(func,gen1,gen2,ttNP) %>% tidyr::unnest(ttNP),
  Y %>% dplyr::select(func,gen1,gen2,vtP) %>% tidyr::unnest(vtP),
  Y %>% dplyr::select(func,gen1,gen2,ttEP) %>% tidyr::unnest(ttEP)
)
dfs %<>% purrr::map(function(x) {
  colnames(x) = toupper(colnames(x))
  x
})
dfs %>% purrr::map(~ .x %>% dplyr::count(INFERENCE))
# toxl(dfs, sheet_names = c("ttest_uneqVaR","var_test","ttest_eqVar"), open_file = FALSE, out_file = "func-tests.xlsx")
toxl(dfs[[1]], sheet_names = "ttest_uneqVaR", open_file = FALSE, out_file = "func-tests.xlsx")
dfs[[1]] %>% dplyr::filter(INFERENCE == "NOTEQ") %>% toxl(open_file = FALSE, out_file = "funcs-diff.xlsx")

# # // [ignore for now] pairwise genotypic difference across functions //
# all_combs = utils::combn(unique(y$gen),2) %>% t() %>% tibble::as_tibble() %>%
#   dplyr::rename(gen1 = V1, gen2 = V2)
# Y = all_combs %>%
#   dplyr::mutate(data = purrr::map2(gen1, gen2, function(g1,g2) y %>% dplyr::filter(gen %in% c(g1,g2))))
# Y$data[[1]] %>% dplyr::count(gen)
# Y$data[[1]] %>% dplyr::count(func)
# Y$data[[1]] %>% dplyr::count(gen, func)

# ----------------------------------------

# # // Prepare the data //
# clr()
# Y = readxl::read_excel("lalit-4-genotypes.xlsx")
# y = Y %>% dplyr::select(1:13)
# # colnames(y) %<>% purrr::map_chr(~ .x %>% stringr::str_replace_all(" ","") %>% stringr::str_replace_all("\\(","_") %>% stringr::str_replace_all("\\)",""))
# colnames(y) %<>% purrr::map_chr(~ ifelse(.x == "function", "func", .x %>% stringr::str_replace_all(" ","") %>% stringr::str_replace_all("\\(","_") %>% stringr::str_replace_all("\\)","")))
# y %>% dplyr::count(func)
# if (nrow(y %>% dplyr::count(func) %>% dplyr::filter(n != 1)) != 0) stop("Functions are not distinct")
# save(y, file = "funcs-data.RData")

# ----------------------------------------
