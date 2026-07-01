knitr::opts_chunk$set(echo = TRUE, message = FALSE)
setwd("D:/18240116syh") 
library(readr)
library(dplyr)
library(TOSTER) # for equivalence tests
library(knitr) # for tables
library(haven) # for reading spss files
library(afex) # for SPSS-style ANOVAs; also loads lme4 and lmerTest for LMEMs
library(tidyverse) # for data wrangling
library(psych) # for Cronbach's alphas

theme_set(theme_minimal(base_size = 12))

rround <- function(x, digits = 3, val = "p", remove0 = TRUE) {
  if (val == "p" & x[1] < .001) return("p < .001")
  
  x <- round(x, digits)
  if (val != "") {
    val <- paste0(sub("%", "%%", val), " = ")
  }
  if (length(x) == 2) {
    txt <- paste0(val, "[%.", digits, "f, %.", digits, "f]") %>%
      sprintf(x[[1]], x[[2]])
  } else {
    txt <- paste0(val, "%.", digits, "f") %>%
      sprintf(x)
  }
  
  if (remove0) {
    gsub("0\\.", ".", txt)
  } else {
    txt
  }
}

t2r <- function(tt, ci = FALSE) {
  if (ci) {
    t <- tt$conf.int/tt$stderr
  } else {
    t <- tt$statistic
  }
  df <- tt$parameter
  r <- sqrt(t^2 / (t^2 + df))
  
  # add sign
  ifelse(t > 0, r, -r)
}

anova_CI <- function(F.value, df.1, df.2, conf.level = .95, digits = 3) {
  ci <- apaTables::get.ci.partial.eta.squared(
    F.value, df.1, df.2, conf.level)
  
  rround(c(ci$LL, ci$UL), digits, "90% CI", FALSE)
}

# report ANOVA stats
report_F <- function(anova_obj, the_term) {
  suppressWarnings(broom::tidy(anova_obj)) %>%
    filter(term == the_term) %>%
    mutate(F = round(statistic, 2) %>% sprintf("%.2f", .),
           p = rround(p.value, 3),
           es = round(pes, 3) %>% sprintf("%.3f", .),
           ci = anova_CI(statistic, num.Df, den.Df, .95),
           txt = paste0("F(", num.Df, ", ", den.Df, ") = ", F,
                        ", $\\eta_{p}^{2}$ ", 
                        ifelse(pes < .001, "< .001", paste0("= ", es)), 
                        ", ", ci, ", ", p)
    ) %>% pull(txt)
}

colspec <- cols(
  response = col_double(),
  .default = col_guess()
)

safe_read <- function(file) {
  tryCatch({
    read_tsv(file, col_types = colspec, show_col_types = FALSE)
  }, error = function(e) {
    warning("文件读取失败:", file)
    return(NULL)
  })
}

data.vu <- list.files(
  path = "rawdata/rawdata/VU", 
  pattern = "dat$", 
  recursive = TRUE, 
  full.names = TRUE
) %>%
  lapply(read_tsv, col_types = colspec, show_col_types = FALSE) %>%
  bind_rows() %>%
  filter(subject != 9999) %>%  ## exclude test run 
  select(subject:stimulusitem5) %>%
  mutate(site = "VU")

data.gl <- list.files(path = "rawdata/rawdata/Glasgow", 
                      pattern = "iqdat$", recursive = TRUE,
                      full.names = TRUE) %>%
  lapply(read_tsv, col_types = colspec, show_col_types = FALSE) %>%
  bind_rows() %>%
  select(subject:stimulusitem5) %>%
  mutate(site = "GL")

colspec <- cols(response = col_double(), .default = col_guess())

mi.excl <- c(9999, 999, 333, 90909, 9994, 999914, 
             999913, 999911, 99998, 999910, 99996, 99999)

data.mi <- list.files(path = "rawdata/rawdata/Michigan", 
                      pattern = "iqdat$", recursive = TRUE,
                      full.names = TRUE) %>%
  lapply(read_tsv, col_type = colspec) %>%
  bind_rows() %>%
  filter(!(subject %in% mi.excl)) %>%
  select(subject:stimulusitem5) %>%
  mutate(site = "MI")

dat.raw <- bind_rows(data.vu, data.gl, data.mi)

trials <- dat.raw %>%
  mutate(
    stim_id = str_replace(stimulusitem3, ".bmp", ""),
    target = str_replace(stimulusitem5, ".bmp", "")
  ) %>%
  select(-(stimulusitem1:stimulusitem5)) %>%
  filter(blockcode != "practice") %>%
  separate(trialcode, c("real", "position", "trial_type", "block")) %>%
  separate(stim_id, c("w", "face_sex", "face_type", "face_n"), sep = c(1, 2,3), remove = F) %>%
  select(-real, -w, -face_n) %>%
  mutate( # circle = 18, square = 23
    correct = case_when(
      target == "circle" & response == 18 ~ 1,
      target == "square" & response == 23 ~ 1,
      TRUE ~ 0
    ),
    face_type = recode(face_type, "d" = "disfigured", "n" = "typical")
  ) %>%
  rename(rt = latency, user_id = subject)

# rater data from qualtrics

raters.vu <- haven::read_sav("rawdata/rawdata/VU/VU Dot Probe Qualtrics.sav") %>%
  mutate(site = "VU")
## fixing a typo, participant number was 534 and not 334
raters.vu$PID[raters.vu$PID == 334] <- 534  

raters.gl <- haven::read_sav("rawdata/rawdata/Glasgow/Glasgow Dot Probe Qualtrics.sav") %>%
  mutate(site = "GL")

raters.mi <- read_sav("rawdata/rawdata/Michigan/Dot Probe Questionnaire - MICHIGAN_September 30, 2019_09.17.sav") %>%
  mutate(site = "MI") %>%
  filter(PID != 9999) # filter out test run

### PVD Subscales

column_labels <- purrr::map(raters.vu, ~attr(., 'label'))

pvd <- column_labels[17:31] %>%
  purrr::map(~str_replace(., "How strongly do you agree or disagree with the following statement\\?\\s*", ""))

tibble(
  "column name" = names(pvd),
  "subscale" = c(1,2,1,1,2,2,1,2,1,2,1,2,1,2,1) %>% 
    recode("1" = "germ aversion", "2" = "perceived infectability"),
  "question" = pvd
) %>%
  knitr::kable()

### NO recoding of the reverse questions, using original questionnaire data
raters <- bind_rows(raters.vu, raters.gl, raters.mi) %>%
  mutate(
    sex = NA_character_,
    age = NA_real_
  ) %>%
  select(
    site, sex, age, 
    user_id = PID,
    recent_cold = DH1,
    recency_1 = DH2,
    recency_2 = DH3,
    recency_3 = DH4,
    recency_4 = DH5,
    PVDpi_1 = PVD2,
    PVDpi_2 = PVD5R,
    PVDpi_3 = PVD6,
    PVDpi_4 = PVD8,
    PVDpi_5 = PVD10,
    PVDpi_6 = PVD12R,
    PVDpi_7 = PVD14R,
    PVDga_1 = PVD1,
    PVDga_2 = PVD3R,
    PVDga_3 = PVD4,
    PVDga_4 = PVD7,
    PVDga_5 = PVD9,
    PVDga_6 = PVD11R,
    PVDga_7 = PVD13R,
    PVDga_8 = PVD15,
    pd_1 = DS1,
    pd_2 = DS2,
    pd_3 = DS3,
    pd_4 = DS4,
    pd_5 = DS5,
    pd_6 = DS6,
    pd_7 = DS7,
    img_1 = A1b___yellow_stuff,
    img_2 = A2b_sickface,
    img_3 = A3b_Crowded_train,
    img_4 = A4b_Yellow_red_towel,
    img_5 = A5b_wound,
    img_6 = A6b_worms,
    img_7 = A7b_mite,
    O1R:H10R,
    WFD1:WMN10
  ) %>%
  rename(WFD10 = WFD10___rating) %>%
  mutate(user_id = as.numeric(user_id))


### Correlation plot for PVD questions

raters %>%
  select(PVDpi_1:PVDga_8) %>%
  cor() %>% # create the correlation matrix
  as.data.frame() %>% # make it a data frame
  rownames_to_column(var = "V1") %>% # set rownames as V1
  gather("V2", "r", PVDpi_1:PVDga_8) %>%
  mutate(V2 = factor(V2),
         V2 = factor(V2, levels = rev(levels(V2)))) %>%
  ggplot(aes(V1, V2, fill=r)) +
  geom_tile(color = "grey") +
  scale_fill_viridis_c() +
  scale_x_discrete(position = "top") +
  labs(x = "", y = "") +
  theme(axis.text.x = element_text(angle = 45, hjust = 0))

### Correlation plot for scores

# average to create scores

raters_scores <- raters %>%
  mutate(
    recent_cold = case_when(
      recent_cold %in% 1:3 ~ "recent",
      recent_cold %in% 4:7 ~ "not recent"
    )
  ) %>%
  gather(var, score, recency_1:WMN10) %>%
  mutate(var = str_replace(var, "WFD", "WFD_0"),
         var = str_replace(var, "WMD", "WMD_0"),
         var = str_replace(var, "WFN", "WFN_0"),
         var = str_replace(var, "WMN", "WMN_0"),
         var = str_replace(var, "(H|E|X|A|C|O)", "hexaco.\\1_0"),
         var = str_replace(var, "_010", "_10"),
         var = str_replace(var, "R$", "")
  ) %>%
  separate(var, c("var", "var_n"), sep = "_") %>%
  group_by(user_id, sex, age, recent_cold, var) %>%
  summarise(score = mean(score)) %>%
  ungroup() %>%
  spread(var, score)

raters_scores %>%
  select(hexaco.A:WMN) %>%
  cor() %>% # create the correlation matrix
  as.data.frame() %>% # make it a data frame
  rownames_to_column(var = "V1") %>% # set rownames as V1
  gather("V2", "r", hexaco.A:WMN) %>%
  mutate(V2 = factor(V2),
         V2 = factor(V2, levels = rev(levels(V2)))) %>%
  ggplot(aes(V1, V2, fill=r)) +
  geom_tile(color = "grey") +
  geom_text(aes(label = round(r, 2)), color = "white") +
  scale_fill_viridis_c() +
  scale_x_discrete(position = "top") +
  labs(x = "", y = "") +
  theme(axis.text.x = element_text(angle = 45, hjust = 0))

### Exclusions

trials_inc <- trials %>%
  filter(correct == 1, trial_type == "diff") %>% 
  group_by(user_id) %>%
  mutate(mean_rt = mean(rt),
         sd_rt = sd(rt),
         n_correct = n()) %>%
  filter(rt <= mean_rt + 3*sd_rt, rt >= mean_rt - 3*sd_rt) %>%
  ungroup() %>%
  group_by() %>%
  filter(n_correct >= mean(n_correct) - 3*sd(n_correct)) %>%
  ungroup() 

### Join Data

# centre and scale continuous scores from only included raters
raters_inc <- raters_scores %>%
  semi_join(trials_inc, by = "user_id") %>%
  mutate(
    recent_cold = as.factor(recent_cold),
    recency.z = (recency - mean(recency)) / sd(recency),
    PVDga.z = (PVDga - mean(PVDga)) / sd(PVDga),
    PVDpi.z = (PVDpi - mean(PVDpi)) / sd(PVDpi),
    img.z = (img - mean(img)) / sd(img),
    pd.z = (pd - mean(pd)) / sd(pd)
  )

data_agg <- trials_inc %>%
  group_by(user_id, face_type, site) %>%
  summarise(rt = mean(rt)) %>%
  ungroup() %>%
  left_join(raters_inc, by = "user_id") %>%
  filter(!is.na(recent_cold)) %>%
  mutate(recent_cold = factor(recent_cold, levels = c("recent", "not recent")))

write_csv(data_agg, "SBVM_data_aggregated.csv")

all.raters <- raters %>% 
  full_join(trials, by = "user_id", 
            suffix = c(".quest", ".exp"))

pre_excl <- all.raters %>% 
  summarise(n = n_distinct(user_id)) %>%
  pull(n)

n_per_site <- all.raters %>% 
  mutate(site = ifelse(is.na(site.quest), site.exp, site.quest)) %>%
  group_by(site) %>% 
  summarise(n = n_distinct(user_id))

pre_excl_gl <- filter(n_per_site, site == "GL") %>% pull(n)
pre_excl_mi <- filter(n_per_site, site == "MI") %>% pull(n)
pre_excl_vu <- filter(n_per_site, site == "VU") %>% pull(n)

error.excl <- trials %>%
  anti_join(trials_inc, by = "user_id") %>%
  group_by(user_id) %>%
  count(user_id)

excl_error <- nrow(error.excl)

raters.excl <- all.raters %>%
  anti_join(data_agg, by = "user_id") %>%
  anti_join(error.excl, by = "user_id")

excl_missing <- raters.excl %>% 
  summarise(n = n_distinct(user_id)) %>% 
  pull(n)

post_excl <- data_agg %>% 
  summarise(n = n_distinct(user_id)) %>%
  pull(n)

post_excl_ill <- data_agg %>% 
  filter(recent_cold == "recent") %>%
  summarise(n = n_distinct(user_id)) %>%
  pull(n)

post_excl_not <- data_agg %>% 
  filter(recent_cold == "not recent") %>%
  summarise(n = n_distinct(user_id)) %>%
  pull(n)

alpha.pi <- select(raters, PVDpi_1:PVDpi_7) %>% alpha()
alpha.ga <- select(raters, PVDga_1:PVDga_8) %>% alpha()
alpha.img <- select(raters, img_1:img_7) %>% alpha()
alpha.pd <- select(raters, pd_1:pd_7) %>% alpha()

a.pi <- round(alpha.pi$total$std.alpha, 2)
a.ga <- round(alpha.ga$total$std.alpha, 2)
a.img <- round(alpha.img$total$std.alpha, 2)
a.pd <- round(alpha.pd$total$std.alpha, 2)

wide_data <- data_agg %>%
  select(user_id, face_type, rt) %>%
  group_by(user_id, face_type) %>%
  summarise(rt = mean(rt), .groups = "drop") %>% # ??????????????????????????????1?????????
  spread(face_type, rt) # ???????????????,???????????????????????????

facetype.cor <- cor.test(~ disfigured + typical, data = wide_data)
                         
print(facetype.cor)


trials_inc %>%
  semi_join(raters, by = "user_id") %>%
  group_by(user_id, face_type) %>%
  summarise(included_trials = n()) %>%
  ungroup() %>%
  group_by(face_type) %>%
  summarise(n = n(),
            mean = mean(included_trials),
            min = min(included_trials),
            max = max(included_trials)
  ) %>%
  mutate_if(is.numeric, round, 2)

trials_inc %>%
  ggplot(aes(rt, color = face_type, group = stim_id)) +
  geom_density() +
  scale_color_manual(values = c("#CC0000", "#006699"), name = "Face Type") +
  theme_minimal()

raters_scores %>%
  gather(variable, score, recency, PVDga, PVDpi) %>%
  ggplot(aes(score, color = variable)) +
  geom_density()

# face type * recent illness mean & sd
ix_desc <- data_agg %>%
  group_by(face_type, recent_cold) %>%
  summarise(n = n(),
            mean = mean(rt),
            sd = sd(rt)
  ) %>%
  ungroup() %>%
  mutate_if(is.numeric, round, 3)


# face type mean & sd
ft_desc <- data_agg %>%
  group_by(face_type) %>%
  summarise(n = n(),
            mean = mean(rt),
            sd = sd(rt)
  ) %>%
  ungroup() %>%
  mutate_if(is.numeric, round, 3)


# recent cold mean & sd
ir_desc <- data_agg %>%
  group_by(recent_cold) %>%
  summarise(n = n()/2,
            mean = mean(rt),
            sd = sd(rt)
  ) %>%
  ungroup() %>%
  mutate_if(is.numeric, round, 3)

ggplot(data_agg, aes(recent_cold, rt, color = face_type)) +
  geom_violin(aes(fill = face_type), alpha = 0.25) +
  stat_summary(geom = "errorbar", width = 0.1,
               position = position_dodge(width = 0.9),
               show.legend = F) +
  stat_summary(geom = "point", fun = "mean", size = 1,
               position = position_dodge(width = 0.9),
               show.legend = F) +
  labs(x = "Illness Recency",
       y = "Mean RT on Incongruent-Location Trials") +
  scale_x_discrete(labels = c("Not Recent", "Recent")) +
  scale_y_continuous(breaks = seq(300, 1800, 300), limits = c(300, 1800)) +
  scale_color_manual(values = c("#CC0000", "#006699"),
                     name = "Face Type", 
                     labels = c("Disfigured", "Typical")) +
  scale_fill_manual(values = c("#CC0000", "#006699"),
                    guide = F) +
  theme_minimal()

ggsave("SBVM_fig1.png", width = 8, height = 6, dpi = 300)

tabledat <- spread(data_agg, face_type, rt) %>% 
  mutate(latency_diff = disfigured - typical,
         recency_cat = recode(recent_cold, 
                              "recent" = 1, 
                              "not recent" = 0)) %>%
  
  select("Illness recency (discrete)" = recency_cat, 
         "Illness recency (continuous)" = recency, 
         "Latencies (difference)" = latency_diff, 
         "Latencies (disfigured)" = disfigured, 
         "Latencies (typical)" = typical, 
         "Pathogen disgust" = pd, 
         "Disgust toward images" = img, 
         "Germ Aversion" = PVDga, 
         "Perceived Infectability" = PVDpi)

means <- summarise_all(tabledat, mean) %>% t()%>% 
  as_tibble(rownames = "A") %>% 
  rename(Mean = V1)
sds <- summarise_all(tabledat, sd) %>% t() %>% 
  as_tibble(rownames = "A") %>% 
  rename(SD = V1)

rhos <- cor(tabledat)
cortable <- rhos
corci <- psych::cor.ci(rhos, n = nrow(tabledat), plot = FALSE)
cortable[upper.tri(rhos)] <- rround(rhos[upper.tri(rhos)], 2, "")
cis <- map2(corci$ci$lower, corci$ci$upper, ~rround(c(.x, .y), 2, "")) %>% unlist()
cortable[lower.tri(cortable)] <- cis
diag(cortable) <- "--"

cortable <- cortable %>% as_tibble(rownames = "A")
names(cortable) <- c("A", 1:9)

# calculate t-test and r from t for discrete illness recency
disc_r <- c()
disc_ci <- c()
for (i in 2:9) {
  ird <- tabledat %>% pull(1)
  contvar <- tabledat %>% pull(i)
  tt <- t.test(contvar[ird==1], contvar[ird==0])
  
  disc_r <- c(disc_r, t2r(tt))
  disc_ci <- c(disc_ci, t2r(tt, ci = TRUE) %>% rround(2, ""))
}

cortable[1, 3:10] <- rround(disc_r, 2, "") %>% as.list()
cortable[2:9, 2] <- disc_ci

descriptives <- left_join(means, sds, by = "A") %>%
  left_join(cortable, by = "A") %>%
  mutate_if(is.numeric, round, 2) %>%
  mutate(A = paste0(row_number(), ". ", A)) %>%
  rename("  " = A)

kable(descriptives, caption="*Table 1* Means, standard deviations, and correlations between study variables. Pearson correlations appear above the diagonal, and 95% confidence intervals appear below the diagonal. Point-biserial correlations derived from t-values are reported between the dichotomous illness recency variable and other variables.", align = "lrrrrrrrrrrr")


cat_ANOVA <- aov_4(rt ~ recent_cold * face_type + 
                     (1 + face_type | user_id), 
                   data = data_agg, factorize = FALSE)

a1a <- anova(cat_ANOVA, es = "pes", MSE = FALSE) %>% print()

recent_ill <- filter(data_agg, recent_cold == "recent")
rt_dis <- filter(recent_ill, face_type == "disfigured")$rt
rt_typ <- filter(recent_ill, face_type == "typical")$rt
a1b <- t.test(rt_dis, rt_typ, paired = TRUE)
a1b_stats <- sprintf("M = %.2f ms, %s, %s",
                     a1b$estimate,
                     rround(a1b$conf.int, 2, "95% CI", FALSE),
                     rround(a1b$p.value, 3, "p"))

a1b

not_recent_ill <- filter(data_agg, recent_cold == "not recent")
rt_dis_healthy <- filter(not_recent_ill, face_type == "disfigured")$rt
rt_typ_healthy <- filter(not_recent_ill, face_type == "typical")$rt
a1c <- t.test(x = rt_dis_healthy, y = rt_typ_healthy, paired = TRUE)

a1c_stats <- sprintf("M = %.2f ms, %s, %s",
                     a1c$estimate,
                     rround(a1c$conf.int, 2, "95% CI", FALSE),
                     rround(a1c$p.value, 3, "p"))

a1c

ft_dis <- ft_desc %>% filter(face_type == "disfigured") %>%
  select(mean, sd) %>% mutate_all(round)

ft_norm <- ft_desc %>% filter(face_type == "typical") %>%
  select(mean, sd) %>% mutate_all(round)

ir_rec <- ir_desc %>% filter(recent_cold == "recent") %>%
  select(mean, sd) %>% mutate_all(round)

ir_not <- ir_desc %>% filter(recent_cold == "not recent") %>%
  select(mean, sd) %>% mutate_all(round)

rc_stats <- report_F(a1a, "recent_cold")
ft_stats <- report_F(a1a, "face_type")
ir_stats <- report_F(a1a, "recent_cold")
ix_stats <- report_F(a1a, "recent_cold:face_type")

stats_cat <- data_agg %>%
  spread(face_type, rt, convert = T) %>%
  mutate(dv = disfigured - typical) %>%
  group_by(recent_cold) %>%
  summarise(n = n(),
            m = mean(dv),
            sd = sd(dv)
  )

n1 <- stats_cat %>% filter(recent_cold == "not recent") %>% pull(n)
n2 <- stats_cat %>% filter(recent_cold == "recent") %>% pull(n)
m1 <- stats_cat %>% filter(recent_cold == "not recent") %>% pull(m)
m2 <- stats_cat %>% filter(recent_cold == "recent") %>% pull(m)
sd1 <- stats_cat %>% filter(recent_cold == "not recent") %>% pull(sd)
sd2 <- stats_cat %>% filter(recent_cold == "recent") %>% pull(sd)
low_eqbound_d <- -.35
high_eqbound_d <- .35

mytost <- TOSTtwo(m1, m2, sd1, sd2, n1, n2, low_eqbound_d, high_eqbound_d)

nhst <- data_agg %>%
  spread(face_type, rt, convert = T) %>%
  mutate(dv = disfigured - typical) %>%
  t.test(dv ~ recent_cold, .)

nhst_null <- sprintf("t(%.1f) = %.2f, %s",
                     round(nhst$parameter, 1),
                     round(nhst$statistic, 2),
                     rround(nhst$p.value, 3))

sdpooled <- sqrt((sd1^2 + sd2^2)/2)
dz <- round((m1 - m2) / sdpooled, 2)

ci <- sprintf("[%.2f - %.2f]", 
              round(mytost$LL_CI_TOST/sdpooled, 2),
              round(mytost$UL_CI_TOST/sdpooled, 2))

lb <- sprintf("t(%.1f) = %.2f, %s",
              round(mytost$TOST_df, 1),
              round(mytost$TOST_t1, 2),
              rround(mytost$TOST_p1, 3))

ub <- sprintf("t(%.1f) = %.2f, %s",
              round(mytost$TOST_df, 1),
              round(mytost$TOST_t2, 2),
              rround(mytost$TOST_p2, 3))

cont_ANOVA <- aov_4(rt ~ recency.z * face_type + 
                      (1 + face_type | user_id), 
                    data = data_agg, factorize = FALSE)

a1a_cont <- anova(cont_ANOVA, es = "pes", MSE = FALSE) %>% print()

irc_stats <- report_F(a1a_cont, "recency.z")
ftc_stats <- report_F(a1a_cont, "face_type")
ixc_stats <- report_F(a1a_cont, "recency.z:face_type")

data_agg_wide <- spread(data_agg, face_type, rt)
a2a <- t.test(data_agg_wide$PVDpi~data_agg_wide$recent_cold) %>% print()

a2b <- t.test(data_agg_wide$PVDga~data_agg_wide$recent_cold) %>% print()

a2c <- t.test(data_agg_wide$pd~data_agg_wide$recent_cold) %>% print()

a2d <- t.test(data_agg_wide$img~data_agg_wide$recent_cold) %>% print()

ill_pi <- sprintf("t(%.1f) = %.2f, %s",
                  round(a2a$parameter, 1),
                  round(a2a$statistic, 2),
                  rround(a2a$p.value, 3))

ill_ga <- sprintf("t(%.1f) = %.2f, %s",
                  round(a2b$parameter, 1),
                  round(a2b$statistic, 2),
                  rround(a2b$p.value, 3))

ill_pd <- sprintf("t(%.1f) = %.2f, %s",
                  round(a2c$parameter, 1),
                  round(a2c$statistic, 2),
                  rround(a2c$p.value, 3))

ill_img <- sprintf("t(%.1f) = %.2f, %s",
                   round(a2d$parameter, 1),
                   round(a2d$statistic, 2),
                   rround(a2d$p.value, 3))


rill_pi <- sprintf("%s, %s, %s",
                   rround(t2r(a2a), 2, "r"),
                   rround(t2r(a2a, ci=TRUE), 2, "CI"),
                   rround(a2a$p.value, 3))

rill_ga <- sprintf("%s, %s, %s",
                   rround(t2r(a2b), 2, "r"),
                   rround(t2r(a2b, ci=TRUE), 2, "CI"),
                   rround(a2b$p.value, 3))

rill_pd <- sprintf("%s, %s, %s",
                   rround(t2r(a2c), 2, "r"),
                   rround(t2r(a2c, ci=TRUE), 2, "CI"),
                   rround(a2c$p.value, 3))

rill_img <- sprintf("%s, %s, %s",
                    rround(t2r(a2d), 2, "r"),
                    rround(t2r(a2d, ci=TRUE), 2, "CI"),
                    rround(a2d$p.value, 3))

a3a <- cor.test(data_agg_wide$recency.z, data_agg_wide$PVDpi.z) %>% print()

a3b <- cor.test(data_agg_wide$recency.z, data_agg_wide$PVDga.z) %>% print()

a3c <- cor.test(data_agg_wide$recency.z, data_agg_wide$pd.z) %>% print()

a3d <- cor.test(data_agg_wide$recency.z, data_agg_wide$img.z) %>% print()

cill_pi <- sprintf("%s, %s, %s",
                   rround(a3a$estimate, 2, "r"),
                   rround(a3a$conf.int, 2, "CI"),
                   rround(a3a$p.value, 3))

cill_ga <- sprintf("%s, %s, %s",
                   rround(a3b$estimate, 2, "r"),
                   rround(a3b$conf.int, 2, "CI"),
                   rround(a3b$p.value, 3))

cill_pd <- sprintf("%s, %s, %s",
                   rround(a3c$estimate, 2, "r"),
                   rround(a3c$conf.int, 2, "CI"),
                   rround(a3c$p.value, 3))

cill_img <- sprintf("%s, %s, %s",
                    rround(a3d$estimate, 2, "r"),
                    rround(a3d$conf.int, 2, "CI"),
                    rround(a3d$p.value, 3))

trials_w <- trials %>%
  filter(correct == 1, trial_type == "diff") %>% 
  group_by(user_id) %>%
  # calculate participant-specific means and SDs
  mutate(mean_rt = mean(rt),
         sd_rt = sd(rt),
         n_correct = n()) %>%
  ungroup() %>%
  # replace >3SD reaction times
  mutate(
    rt = case_when(
      rt > mean_rt + 3*sd_rt ~ mean_rt + 3*sd_rt,
      rt < mean_rt - 3*sd_rt ~ mean_rt - 3*sd_rt,
      TRUE ~ rt
    )
  ) %>%
  group_by() %>%
  filter(n_correct >= mean(n_correct) - 3*sd(n_correct)) %>%
  ungroup() 

# centre and scale continuous scores from only included raters
raters_w <- raters_scores %>%
  semi_join(trials_w, by = "user_id") %>%
  mutate(
    recent_cold = as.factor(recent_cold),
    recency.z = (recency - mean(recency)) / sd(recency),
    PVDga.z = (PVDga - mean(PVDga)) / sd(PVDga),
    PVDpi.z = (PVDpi - mean(PVDpi)) / sd(PVDpi),
    img.z = (img - mean(img)) / sd(img),
    pd.z = (pd - mean(pd)) / sd(pd)
  )

data_agg_w <- trials_w %>%
  group_by(user_id, face_type, site) %>%
  summarise(rt = mean(rt)) %>%
  ungroup() %>%
  left_join(raters_w, by = "user_id") %>%
  filter(!is.na(recent_cold))

cat_ANOVA_w <- aov_4(rt ~ recent_cold * face_type + 
                       (1 + face_type | user_id), 
                     data = data_agg_w, factorize = FALSE)

anova(cat_ANOVA_w, es = "pes", MSE = FALSE)

cont_ANOVA_w <- aov_4(rt ~ recency.z * face_type + 
                        (1 + face_type | user_id), 
                      data = data_agg_w, factorize = FALSE)

anova(cont_ANOVA_w, es = "pes", MSE = FALSE)

cat_ANCOVA <- aov_4(disfigured ~ typical + recent_cold + 
                      (1 | user_id), 
                    data = data_agg_wide, factorize = FALSE)

anc <- anova(cat_ANCOVA, es = "pes", MSE = FALSE) %>% print()

rc <- anc[2, ]
cat_stats <- sprintf("F(%.0f, %.0f) = %.2f, %s",
                     rc$`num Df`, rc$`den Df`, 
                     round(rc$F, 2), 
                     rround(rc$`Pr(>F)`, 3))

# continuous illness variable
cont_ANCOVA <- aov_4(disfigured ~ typical + recency.z + 
                       (1 | user_id), 
                     data = data_agg_wide, factorize = FALSE)

canc <- anova(cont_ANCOVA, es = "pes", MSE = FALSE) %>% print()

crc <- canc[2, ]
cont_stats <- sprintf("F(%.0f, %.0f) = %.2f, %s",
                      crc$`num Df`, crc$`den Df`, 
                      round(crc$F, 2), 
                      rround(crc$`Pr(>F)`, 3))

aov_4(rt ~ recent_cold * face_type + 
        (1 + face_type | user_id), 
      data = filter(data_agg, site == "MI"), 
      factorize = FALSE) %>%
  anova(es = "pes", MSE = FALSE)

aov_4(rt ~ recency.z * face_type + 
        (1 + face_type | user_id), 
      data = filter(data_agg, site == "MI"), 
      factorize = FALSE) %>%
  anova(es = "pes", MSE = FALSE)

aov_4(rt ~ recent_cold * face_type + 
        (1 + face_type | user_id), 
      data = filter(data_agg, site == "VU"), 
      factorize = FALSE) %>%
  anova(es = "pes", MSE = FALSE)

aov_4(rt ~ recency.z * face_type + 
        (1 + face_type | user_id), 
      data = filter(data_agg, site == "VU"), 
      factorize = FALSE) %>%
  anova(es = "pes", MSE = FALSE)

aov_4(rt ~ recent_cold * face_type + 
        (1 + face_type | user_id), 
      data = filter(data_agg, site == "GL"), 
      factorize = FALSE) %>%
  anova(es = "pes", MSE = FALSE)

aov_4(rt ~ recency.z * face_type + 
        (1 + face_type | user_id), 
      data = filter(data_agg, site == "GL"), 
      factorize = FALSE) %>%
  anova(es = "pes", MSE = FALSE)

data_all <- trials_inc %>%
  left_join(raters_inc, by = "user_id") %>%
  filter(!is.na(recent_cold)) %>%
  mutate(recent_cold.e = recode(recent_cold, 
                                "not recent" = -0.5, 
                                "recent" = 0.5),
         face_type.e = recode(face_type, 
                              "disfigured" = -0.5, 
                              "typical" = 0.5)
  )

lmer_cat <- lmer(rt ~ recent_cold.e * face_type.e + 
                   (1 | site) +
                   (1 | user_id) +
                   (1 | stim_id),
                 data = data_all)

summary(lmer_cat)

lmer_cont <- lmer(rt ~ recency.z * face_type.e + 
                    (1 | site) +
                    (1 | user_id) +
                    (1 | stim_id),
                  data = data_all)

summary(lmer_cont)

hexaco <- purrr::map_df(c("H", "E", "X", "A", "C", "O"), function(l) {
  var <- paste0("hexaco.", l)
  t.test(data_agg_wide[[var]]~data_agg_wide$recent_cold) %>%
    broom::tidy() %>%
    mutate(hexaco = l)
}) %>%
  select(hexaco, estimate1:parameter) %>%
  mutate(r = t2r(list(statistic = statistic, parameter = parameter))) %>%
  mutate_if(is.numeric, signif, 3) %>%
  mutate(statistic = -statistic, # switch sign for t-test order
         r = round(r, 2)) %>%
  rename(t = statistic,
         p = p.value,
         df = parameter,
         "not recent" = estimate1,
         "recent" = estimate2)

kable(hexaco, caption="*Table X* Illness recency effects on the six HEXACO factors")

facetype.cor.est <- rround(facetype.cor$estimate, 2, "r")

kable(descriptives, caption="*Table 1* Means, standard deviations, and correlations between study variables. Pearson correlations appear above the diagonal, and 95% confidence intervals appear below the diagonal. Point-biserial correlations derived from t-values are reported between the dichotomous illness recency variable and other variables.", align = "lrrrrrrrrrrr")
library(tidyverse)
library(writexl)
library(broom)
write_csv(data_agg_w, "D:/18240116syh/final_analysis_data.csv")
write_xlsx(data_agg_w, "D:/18240116syh/final_analysis_data.xlsx")
write_csv(trials_w, "D:/18240116syh/cleaned_trials.csv")
write_csv(raters_w, "D:/18240116syh/cleaned_questionnaire.csv")
write_csv(tidy(a2a), "D:/18240116syh/t_test_PVDpi.csv")
write_csv(tidy(a2b), "D:/18240116syh/t_test_PVDga.csv")
write_csv(tidy(a1a), "D:/18240116syh/anova_results.csv")
write_csv(tidy(a1a_cont), "D:/18240116syh/anova_continuous_results.csv")
tost_df <- data.frame(
  TOST_t1 = mytost$TOST_t1, TOST_p1 = mytost$TOST_p1,
  TOST_t2 = mytost$TOST_t2, TOST_p2 = mytost$TOST_p2,
  TOST_df = mytost$TOST_df,
  LL_CI_TOST = mytost$LL_CI_TOST, UL_CI_TOST = mytost$UL_CI_TOST,
  NHST_t = mytost$NHST_t, NHST_p = mytost$NHST_p
)
write_csv(tost_df, "D:/18240116syh/tost_equivalence_result.csv")
write_csv(tidy(ft_dis), "D:/18240116syh/desc_ft_dis.csv")
write_csv(tidy(ft_norm), "D:/18240116syh/desc_ft_norm.csv")
write_csv(tidy(ir_rec), "D:/18240116syh/desc_ir_rec.csv")
write_csv(tidy(ir_not), "D:/18240116syh/desc_ir_not.csv")
save.image("D:/18240116syh/FULL_WORKSPACE.RData")