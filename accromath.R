setwd(this.path::here())
library(dplyr)
library(ggplot2)
library(patchwork)
library(scales)
library(lubridate)
# remotes::install_github("lbelzile/longevity")
library(longevity) # version 1.1.2
theme_set(theme_classic())

# Base de données sur la longévité canadienne
# http://www.bdlc.umontreal.ca/bdlc/pres-donnees.htm
#
# International Database on Longevity
# https://www.supercentenarians.org/fr/

# Données dans la version locale du paquet longévité
data(idl, package = "longevity")
IDL_QC <- idl |> dplyr::filter(country == "QC")
fit <- with(IDL_QC, fit_elife(time = ndays/365.25,
                       ltrunc = cbind(ltrunc1, ltrunc2)/365.25,
                       rtrunc = cbind(rtrunc1, rtrunc2)/365.25,
                       event = 1,
                       thresh = 106,
                       family = "exp",
                       export = TRUE))
# Estimation des paramètres
summary(fit)
coef(fit)
# IC 95% de Wald
round(coef(fit) + qnorm(c(0.025,0.975)) * fit$std.error, 3)

# Vraisemblance incorrecte (assume échantillon aléatoire simple)
fit_inc <- with(
  IDL_QC,
  fit_elife(time = ndays/365.25,
            event = 1,
            thresh = 106,
            family = "exp",
            export = TRUE))
round(coef(fit_inc) + qnorm(c(0.025,0.975)) * fit_inc$std.error, 3)
profile106_inc <- with(
  IDL_QC,
  prof_exp_scale(time = ndays/365.25,
                 event = 1,
                 thresh = 106))
round(profile106_inc,3)

plot(fit, plot.type = "ggplot", which.plot = "qq")

ecdf <- with(IDL_QC, np_elife(time = ndays/365.25,
                       ltrunc = ltrunc1/365.25,
                       rtrunc = rtrunc1/365.25,
                       event = 1,
                       thresh = 106))
plot(ecdf$cdf)



# Graphique 1 - Calcul du nombre de centenaires au fil des années, données de BDLC
BDLC <- read.table("https://www.prdh.umontreal.ca/BDLC/data/que/Population.txt", skip = 2, header = TRUE) |>
  mutate(Age = as.integer(ifelse(Age == "110+", "110", Age)))
BDLC_pop_cent <- BDLC |>
  dplyr::filter(Age >= 100) |>
  group_by(Year) |>
  summarize(hommes = sum(Male), femmes = sum(Female), total = sum(Total))
BDLC_pop_cent_long <- BDLC_pop_cent  |>
  tidyr::pivot_longer(cols = c(hommes, femmes, total), names_to = "type", values_to = "decompte") |>
  rename(annees = Year) |>
  mutate(type = factor(type))

u <- 100
BDLC_pop_prop <- BDLC |>
  group_by(Year) |>
  summarize(total = sum(Total[Age >= u]) / sum(Total[Age < u]),
            hommes = sum(Male[Age >= u]) / sum(Male[Age < u]),
            femmes = sum(Female[Age >= u]) / sum(Female[Age < u])) |>
  rename(annees = Year) |>
  tidyr::pivot_longer(cols = c(hommes, femmes, total), names_to = "type", values_to = "decompte") |>
  mutate(type = factor(type))



g1 <- ggplot(data = BDLC_pop_cent_long, mapping = aes(x = annees, y = decompte, color = type)) +
  geom_line() +
  scale_y_log10() +
  MetBrewer::scale_color_met_d(name = "Hiroshige") +
  labs(x = "année",
       y = "",
       subtitle = "Nombre de centenaires au Québec (échelle log)",
       color = "",
       caption = "source: BDLC, décembre 2023") +
  theme(legend.position = "bottom",
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_line(colour = "gray", size = 0.1),
        panel.grid.minor.y = element_blank())
g1
ggsave(filename = "../Figures/graphique1.pdf", width = 5, height = 3)
ggsave(filename = "../Figures/graphique1.png", dpi = 300, width = 5, height = 3)

g5 <- ggplot(data = BDLC_pop_prop, mapping = aes(x = annees, y = decompte, color = type)) +
  geom_line() +
  scale_y_log10(breaks = c(1e-3, 1e-4, 1e-5, 1e-6),
                labels = trans_format("log10", math_format(10^.x))) +
  MetBrewer::scale_color_met_d(name = "Hiroshige") +
  labs(x = "année",
       y = "",
       subtitle = "Proportion de centenaires au Québec (échelle log)",
       color = "",
       caption = "source: BDLC, décembre 2023") +
  theme(legend.position = "bottom",
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_line(colour = "gray", size = 0.1),
        panel.grid.minor.y = element_blank()
  )
g5
ggsave(filename = "../Figures/graphique5.pdf", width = 5, height = 3)
ggsave(filename = "../Figures/graphique5.png", dpi = 300, width = 5, height = 3)

(g1 +
    labs(caption = "",
         x = "") +
    theme(legend.position = "none")) / g5
ggsave(filename = "../Figures/graphique1et5.pdf", width = 5, height = 7)
ggsave(filename = "../Figures/graphique1et5.png", dpi = 300, width = 5, height = 7)


# Graphique 2 - Calcul des décomptes moyens et comparaison avec la loi exponentielle
th <- 106
IDL_QC_exc <- IDL_QC |> dplyr::filter(ndays > 365.25*th)
nexc <- nrow(IDL_QC_exc)
ecdf <- with(IDL_QC_exc,
             np_elife(time = ndays/365.25 - th,
                      ltrunc = pmax(0, ltrunc1/365.25 - th),
                      rtrunc = rtrunc1/365.25 - th,
                      event = 1))

pts <- 1:7
opt_rate <- 1/coef(fit)
df_g2 <- data.frame(
  x = c(105 + pts, 105 + pts),
  y = nexc*c((ecdf$cdf(pts)-ecdf$cdf(pts-1)),
             pexp(pts, rate = 1/coef(fit)) - pexp(pts-1, rate = 1/coef(fit))),
  type = factor(rep(c("empirique", "exponentiel"), each = length(pts))))

g2 <- ggplot(data = df_g2,
             mapping = aes(x = x,
                           y = y,
                           group = type,
                           fill = type)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_y_continuous(expand = expansion(add = c(0, NA))) +
  scale_x_continuous(breaks = 105+pts) +
  labs(x = "âge",
       y = "",
       subtitle = "Décompte moyen du nombre de morts par tranche d'âge",
       fill = "",
       caption = "source: IDL") +
  MetBrewer::scale_fill_met_d(name = "Hiroshige") +
  theme(legend.position = "bottom",
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_line(colour = "gray", size = 0.1),
        panel.grid.minor.y = element_blank())
g2
ggsave(filename = "../Figures/graphique2.pdf", width = 5, height = 5)
ggsave(filename = "../Figures/graphique2.png", dpi = 300, width = 5, height = 5)


# Test du khi-deux d'ajustement, en combinant les décomptes dans six catégories de longueur 1 an
Ex <- nexc * diff(pexp(c(0:5,Inf), rate = 1/coef(fit)))
Ob <- nexc * diff(ecdf$cdf(c(0:5,Inf)))
# https://www.prdh.umontreal.ca/BDLC/data/can/Deaths_lexis.txt

pchisq(sum((Ob-Ex)^2/Ex), df = length(Ex)-1, lower.tail = FALSE)


# Graphique 3:
# Données simulées, avec les taux d'entrée du nombre de personnes au delà de 106 ans
# On utilise la population canadienne pour avoir une estimation plus fiable du taux
set.seed(202411)
nabove <- read.table("https://www.prdh.umontreal.ca/BDLC/data/que/Population.txt",
                     skip = 2, header = TRUE) |>
  mutate(Age = as.integer(ifelse(Age == "110+", "110", Age))) |>
  dplyr::filter(Age > 105) |>
  group_by(Year) |>
  summarize(total = round(sum(Total)))
nsim <- sum(nabove$total)
# Contenant pour dates simulées et âges simulés
dates_sim <- Date(length = nsim)
age_deces_sim <- numeric(length = nsim)
# Compteur pour nombre
cpt <- 0
# Simuler trajectoires exponentielles
for(i in seq_len(nrow(nabove))){
  # Nombre de personnes atteignant l'âge de 106 ans
  nobs <- nabove$total[i]
  if(nobs > 0){
    # Simuler une trajectoire exponentielle
    age_deces_sim[cpt + 1:nobs] <- rexp(nobs, rate = 1/fit$par)
    # Simuler des dates d'anniversaire uniformes
    dates_sim[cpt + 1:nobs] <- lubridate::ymd(paste0(nabove$Year[i], "-01-01")) +
      sample(1:365, size = nobs, replace = TRUE)
    # Incrémenter le compteur
    cpt <- cpt + nobs
  }
}
# Diagramme de Lexis
df_g3 <- data.frame(date = dates_sim,
                    age = age_deces_sim) |>
  mutate(mort = pmin(ymd("2019-12-31"), date + lubridate::days(round(365.25*age))),
         fenetre = factor(ifelse(
           mort >= ymd("1985-01-01") & mort <= ymd("2010-01-01"),
           "observé",
           "fantôme")))
g3 <- ggplot(data = df_g3 |>
         dplyr::filter(date > lubridate::ymd("1975-01-01"))) +
  geom_segment(mapping = aes(x = date,
                             xend = mort,
                             y = 106,
                             yend = pmin(106 + age, 106 + as.numeric(ymd("2019-12-31") - date)/365.25),
                             col = fenetre),
               # alpha = 0.5,
               show.legend = FALSE) +
  geom_vline(xintercept = ymd(c("1985-01-01","2009-12-31"))) +
  scale_color_grey(start = 0.8, end = 0.2) +
  scale_x_date(limits = ymd(c("1975-01-01", "2020-01-01")),
               minor_breaks = NULL,
               oob = oob_keep) +
  scale_y_continuous(breaks = 106:115,
                     limits = c(106, NA),
                     minor_breaks = NULL,
                     labels = 106:115,
                     expand = expansion(add = 0),
                     oob = oob_keep) +
  labs(x = "année", y = "", subtitle = "Diagramme de Lexis: âge jusqu'au décès") +
  theme_classic()
g3
ggsave(filename = "../Figures/graphique3.pdf", width = 5, height = 3)
ggsave(filename = "../Figures/graphique3.png", dpi = 300, width = 5, height = 3)

# Graphique 4: Impact de la troncature sur les maximums observés

IDL_QC_max <- idl |>
  dplyr::filter(country == "QC" & ddate >= ymd("1985-01-01")) |>
  mutate(byear = year(bdate)) |>
  group_by(byear) |>
  summarize(max_age = max(ndays)/365.25)


g4 <- ggplot() +
  geom_point(data = IDL_QC |>
               dplyr::filter(ddate >= ymd("1985-01-01")),
             mapping = aes(y = ndays/365.25, x = year(bdate)),
             shape = 4, color = "grey") +
  geom_point(data = IDL_QC_max,
             mapping = aes(y = max_age, x = byear),
             shape = 4) +
  geom_abline(slope = -1, intercept = 2010) +
  geom_abline(slope = -1, intercept = 1984) +
  scale_y_continuous(minor_breaks = 105:115, breaks = seq(106L, 114L, by = 2), oob = oob_keep) +
  scale_x_continuous(breaks = seq(1875L, 1905L, by= 5L), minor_breaks = NULL, limits = c(1875,1906), oob = oob_keep) +
  labs(x = "année de naissance",
       y = "",
       subtitle = "Âge au décès",
       caption = "source: IDL (INSQ)") +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_line(colour = "gray", size = 0.1),
        panel.grid.minor.y = element_blank())
g4
ggsave(filename = "../Figures/graphique4.pdf", width = 5, height = 3)
ggsave(filename = "../Figures/graphique4.png", dpi = 300, width = 5, height = 3)

# Calcul de l'intervalle de confiance (basé sur le rapport de vraisemblance)
profile106 <- with(IDL_QC, prof_exp_scale(time = ndays/365.25,
                                          ltrunc = ltrunc1/365.25,
                                          rtrunc = rtrunc1/365.25,
                                          event = 1,
                                          thresh = 106))
profile_courbe <- with(IDL_QC, prof_exp_scale(time = ndays/365.25,
                                              ltrunc = ltrunc1/365.25,
                                              rtrunc = rtrunc1/365.25,
                                              event = 1,
                                              thresh = 106, confint = FALSE))

g6 <- ggplot() +
  geom_line(data = with(profile_courbe, data.frame(x = psi, y = pll)),
            mapping = aes(x = x, y = y)) +
  geom_hline(yintercept = -qchisq(0.95,1)/2, linetype = "dashed", alpha = 0.1) +
  geom_vline(xintercept = as.numeric(unlist(profile106)), linetype = "dashed", alpha = 0.1) +
  geom_rug(sides = "b",
           data = data.frame(x = as.numeric(unlist(profile106))),
           mapping = aes(x = x)) +

  scale_x_continuous(limits = c(1.25, 2), breaks = seq(1.2, 2, by = 0.1)) +
  scale_y_continuous(limits = c(-7, 0), breaks = seq(-6L, 0, by = 2L), oob = oob_keep) +
  labs(x = expression(theta), y = "", subtitle = "Fonction de log vraisemblance") +
  theme_classic()
g6
ggsave(filename = "../Figures/graphique6.pdf", width = 5, height = 3)
ggsave(filename = "../Figures/graphique6.png", dpi = 300, width = 5, height = 3)


g6b <- ggplot() +
  geom_line(data = with(profile_courbe, data.frame(x = psi, y = pll)),
            mapping = aes(x = x, y = y)) +
  # geom_hline(yintercept = -qchisq(0.95,1)/2, linetype = "dashed", alpha = 0.1) +
  geom_rug(sides = "b",
           data = data.frame(x = as.numeric(unlist(profile106))),
           mapping = aes(x = x)) +
  scale_y_continuous(limits = c(-7, 0), breaks = seq(-10L, 0, by = 2L), expand = expansion(add = 0), oob = oob_keep) +
  scale_x_continuous(limits = c(1.25, 2), breaks = seq(1.2, 2, by = 0.1)) +
  labs(x = expression(theta), y = "", subtitle = "Fonction de log vraisemblance") +
  theme_classic()
  g6b
ggsave(filename = "../Figures/graphique6b.pdf", width = 5, height = 3)
ggsave(filename = "../Figures/graphique6b.png", dpi = 300, width = 5, height = 3)




# Probabilité de vivre au moins un an de plus (pièce de monnaie) pour les supercentenaires
exp(-1/profile106)
# Fonction de risque
1/profile106
# Moyenne (âge en excès moyen)
profile106 + 106
# Calcul de l'intervalle de confiance (basé sur le rapport de vraisemblance)
profile110 <- with(IDL_QC, prof_exp_scale(time = ndays/365.25,
                                          ltrunc = ltrunc1/365.25,
                                          rtrunc = rtrunc1/365.25,
                                          event = 1,
                                          thresh = 110))

# Probabilité de vivre un an ou moins (pièce de monnaie) pour les supercentenaires
1-exp(-1/profile110)
# Fonction de risque
1/profile110
# Moyenne (âge en excès moyen)
profile110 + 110
