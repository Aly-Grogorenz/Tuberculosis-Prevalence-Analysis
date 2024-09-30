######Loading Packages#####
set.seed(1234)
library(Rmisc)
library(tidyverse)
library(Amelia)
library(GGally)
library(lme4)
library(lmerTest)
library(sjstats)
library(MuMIn)
library(sjPlot)
library(performance)
library(caret)
library(flextable)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)

######Read and clean up individual dataframes#######
population_country <- read_csv("population_density.csv")
incidence_country <- read_csv("incidence-by-country.csv")
funding_country <- read_csv("healthfund-by-country.csv", skip = 1)
gdp_country <- read_csv("gdp-by-country.csv")

# clean up data
pop_c <- population_country
pop_c <- pop_c %>% 
  select(-"Series Name", -"Series Code", -"Time Code",
         -"Upper middle income [UMC]", -"Low income [LIC]", 
         -"Lower middle income [LMC]", -"High income [HIC]") %>% 
  gather(Country, Population, -Time) %>% 
  mutate(Population = as.numeric(Population)) %>% 
  filter(!is.na(Time) | !is.na(Population)) %>% 
  dplyr::rename(Year = Time)
# how many countries
n_distinct(pop_c$Country)

gdp_c <- gdp_country
gdp_c <- gdp_c %>% 
  select(-"Series Name", -"Series Code") %>% 
  gather(Year, GDP, -"Country Name", -"Country Code") %>% 
  dplyr::rename(Country = "Country Name", Code = "Country Code") %>% 
  mutate(GDP = as.numeric(GDP)) %>% 
  filter(!is.na(GDP) | !is.na(Country))
n_distinct(gdp_c$Country)

fund_c <- funding_country
fund_c <- fund_c %>% 
  dplyr::rename(Country = "Countries, territories and areas") %>% 
  gather(Year, Funding, -Country) %>% 
  mutate(Year = as.numeric(Year))
n_distinct(fund_c$Country)

inc_c <- incidence_country
inc_c <- inc_c %>% 
  dplyr::rename(Country = "Countries, territories and areas", 
                Incidence = "Incidence of tuberculosis (per 100 000 population per year)") %>% 
  select(Year, Country, Incidence)
n_distinct(inc_c$Country)

# add country codes
# make conversion table
count_code <- substr(pop_c$Country, nchar(pop_c$Country)-3, nchar(pop_c$Country)-1)
countries <- substr(pop_c$Country, 1, nchar(pop_c$Country)-6)

conversion <- tibble(Code=unique(count_code), Country=unique(countries))

# combare fund_c with conversion
anti_join(fund_c, conversion, by="Country")
# add countries not in original conversion table
conversion <- conversion %>% 
  add_row(Code = "COK", Country = "Cook Islands") %>%   #made-up code
  add_row(Code = "NUE", Country = "Niue")  #made-up code
# duplicate conversions: different names of countries already in the list
dup_conversion <- conversion %>% 
  add_row(Code = "BHS", Country = "Bahamas") %>% 
  add_row(Code = "BOL", Country = "Bolivia (Plurinational State of)") %>% 
  add_row(Code = "COG", Country = "Congo") %>% 
  add_row(Code = "COD", Country = "Democratic Republic of the Congo") %>% 
  add_row(Code = "EGY", Country = "Egypt") %>% 
  add_row(Code = "GMB", Country = "Gambia") %>% 
  add_row(Code = "IRN", Country = "Iran (Islamic Republic of)") %>% 
  add_row(Code = "KGZ", Country = "Kyrgyzstan") %>% 
  add_row(Code = "LAO", Country = "Lao People's Democratic Republic") %>% 
  add_row(Code = "FSM", Country = "Micronesia (Federated States of)") %>% 
  add_row(Code = "NLD", Country = "Netherlands (Kingdom of the)") %>% 
  add_row(Code = "KOR", Country = "Republic of Korea") %>% 
  add_row(Code = "MDA", Country = "Republic of Moldova") %>% 
  add_row(Code = "KNA", Country = "Saint Kitts and Nevis") %>% 
  add_row(Code = "LCA", Country = "Saint Lucia") %>% 
  add_row(Code = "VCT", Country = "Saint Vincent and the Grenadines") %>% 
  add_row(Code = "SVK", Country = "Slovakia") %>% 
  add_row(Code = "GBR", Country = "United Kingdom of Great Britain and Northern Ireland") %>% 
  add_row(Code = "TZA", Country = "United Republic of Tanzania") %>% 
  add_row(Code = "USA", Country = "United States of America") %>% 
  add_row(Code = "VEN", Country = "Venezuela (Bolivarian Republic of)") %>%
  add_row(Code = "YEM", Country = "Yemen") 
# still mismatch?
anti_join(fund_c, dup_conversion, by="Country")

# test on inc_c:
anti_join(inc_c, dup_conversion, by="Country")
dup_conversion <- dup_conversion %>% 
  add_row(Code = "PRK", Country = "Democratic People's Republic of Korea")
anti_join(inc_c, dup_conversion, by="Country")

# test on gdp_c:
anti_join(gdp_c, dup_conversion, by = "Code")

# clean up pop_c: split Country into Country & Code:
pop_c <- pop_c %>% 
  mutate(Code = substr(pop_c$Country, nchar(pop_c$Country)-3, nchar(pop_c$Country)-1),
         Country = countries <- substr(pop_c$Country, 1, nchar(pop_c$Country)-6))

# add country code to funding & incidence
fund_c <- left_join(fund_c, dup_conversion, by = "Country")
inc_c <- left_join(inc_c, dup_conversion, by="Country")

# get rid off [] brackets
# a function that takes a string and removes anything after a space
rm_brackets <- function(string){
  s <- ""   
  for (l in 1:nchar(string)){   # go through every letter
    letter <- substr(string, l, l)   #l-th character
    if (letter == " ") {   #if letter is a space, stop there
      return(s)
    } else {
      s <- paste(s, letter, sep="")
    }
  }
}

# incidence
for (i in 1:length(inc_c$Incidence)){
  inc_c$Incidence[i] <- rm_brackets(inc_c$Incidence[i])
}
inc_c <- mutate(inc_c, Incidence = as.numeric(Incidence))

# gdp:
for (i in 1:length(gdp_c$Year)){
  gdp_c$Year[i] <- rm_brackets(gdp_c$Year[i])
}
gdp_c <- mutate(gdp_c, Year = as.numeric(Year))


###### Combine data into one dataframe #######
# left join of all other data, because data is irrelevant if there
# is no incidence in the country
data <- inc_c %>% 
  select(Year, Code, Incidence) %>% 
  left_join(select(fund_c, -Country), by = c("Code", "Year")) %>% 
  left_join(select(gdp_c, -Country), by = c("Code", "Year")) %>% 
  left_join(select(pop_c, -Country), by = c("Code", "Year")) %>% 
  left_join(conversion, by = "Code") %>%  # add country names
  dplyr::rename(PopDensity = Population)

# add regions:
# make vectors for each WHO region containing all relevant countries
africa <- c("Algeria", "Angola", "Benin", "Botswana", "Burkina Faso", 
            "Burundi", "Cameroon", "Cabo Verde", "Central African Republic", 
            "Chad", "Comoros", "Cote d'Ivoire", 
            "Democratic Republic of the Congo", "Equatorial Guinea", 
            "Eritrea", "Ethiopia", "Gabon", "Gambia", "Ghana", "Guinea", 
            "Guinea-Bissau", "Kenya", "Lesotho", "Liberia", "Madagascar", 
            "Malawi", "Mali", "Mauritania", "Mauritius", "Mozambique", 
            "Namibia", "Niger", "Nigeria", "Congo", 
            "Rwanda", "Sao Tome and Principe", "Senegal", "Seychelles", 
            "Sierra Leone", "South Africa", "South Sudan", "Eswatini", 
            "Togo", "Uganda", "Tanzania", "Zambia", "Zimbabwe")
americas <- c("Antigua and Barbuda", "Argentina", "Bahamas", "Barbados", 
              "Belize", "Bolivia", "Brazil", "Canada", "Chile", "Colombia", 
              "Costa Rica", "Cuba", "Dominica", "Dominican Republic", "Ecuador", 
              "El Salvador", "Grenada", "Guatemala", "Guyana", "Haiti", 
              "Honduras", "Jamaica", "Mexico", "Nicaragua", "Panama", "Paraguay", 
              "Peru", "Saint Kitts and Nevis", "Saint Lucia", 
              "Saint Vincent and the Grenadines", "Suriname", 
              "Trinidad and Tobago", "United States", "Uruguay", "Venezuela, RB")
seasia <- c("Bangladesh", "Bhutan", "Korea, Dem. People's Rep.", "India", "Indonesia", 
            "Maldives", "Myanmar", "Nepal", "Sri Lanka", "Thailand", 
            "Timor-Leste")
europe <- c("Albania", "Andorra", "Armenia", "Austria", "Azerbaijan", 
            "Belarus", "Belgium", "Bosnia and Herzegovina", "Bulgaria", 
            "Croatia", "Cyprus", "Czechia", "Denmark", "Estonia", "Finland", 
            "France", "Georgia", "Germany", "Greece", "Hungary", "Iceland", 
            "Ireland", "Israel", "Italy", "Kazakhstan", "Kyrgyzstan", "Latvia", 
            "Lithuania", "Luxembourg", "Malta", "Moldova", "Monaco", "Montenegro", 
            "Netherlands", "North Macedonia", "Norway", "Poland", "Portugal", 
            "Romania", "Russian Federation", "San Marino", "Serbia", "Slovakia", "Slovenia", 
            "Spain", "Sweden", "Switzerland", "Tajikistan", "Turkiye", 
            "Turkmenistan", "Ukraine", "United Kingdom", "Uzbekistan")
eastmediterranean <- c("Afghanistan", "Bahrain", "Djibouti", "Egypt", 
                       "Iran, Islamic Rep.", "Iraq", "Jordan", 
                       "Kuwait", "Lebanon", "Libya", "Morocco", "Oman", 
                       "Pakistan", "Qatar", "Saudi Arabia", "Somalia", "Sudan", 
                       "Syrian Arab Republic", "Tunisia", "United Arab Emirates", "Yemen")
westpacific <- c("Australia", "Brunei Darussalam", "Cambodia", "China", "Cook Islands", 
                 "Fiji", "Japan", "Kiribati", "Lao PDR", "Malaysia", 
                 "Marshall Islands", "Micronesia, Fed. Sts.", "Mongolia", "Nauru", 
                 "New Zealand", "Niue", "Palau", "Papua New Guinea", 
                 "Philippines", "Samoa", "Singapore", "Solomon Islands", 
                 "Republic of Korea", "Tonga", "Tuvalu", "Vanuatu", "Viet Nam")

#make label vector
reg <- c(rep("Africa", length.out = length(africa)), 
         rep("Americas", length.out = length(americas)),
         rep("South-East Asia", length.out = length(seasia)), 
         rep("Europe", length.out = length(europe)), 
         rep("Eastern Mediterranean", length.out = length(eastmediterranean)),
         rep("Western Pacific", length.out = length(westpacific)))
# make country vector
cs <- append(africa, americas) %>% 
  append(seasia) %>% 
  append(europe) %>% 
  append(eastmediterranean) %>% 
  append(westpacific)

# make regions vector:
regions <- tibble(reg, cs) %>% 
  dplyr::rename(Region = "reg", Country = "cs")

# is every country in regions in the conversion table?
anti_join(regions, dup_conversion, by = "Country")

# add country code to regions table
regions <- left_join(regions, dup_conversion, by = "Country")

# add region to data:
data <- left_join(data, select(regions, -Country), by = "Code")

# make relevant variables factors under final variable
tub <- data %>% 
  mutate_if(is.character, as.factor) %>% 
  select(-Code) %>% 
  filter(Year < 2022)
summary(tub)

# all values have massive ranges -> recode as log?
log_tb <- tub %>% 
  mutate(LogInc = log(Incidence),
         LogGDP = log(GDP),
         LogFund = log(Funding),
         LogDen = log(PopDensity)) %>% 
  select(-Incidence, -GDP, -Funding, -PopDensity)
summary(log_tb)
# -Inf in Incidence:
# +1 added to all Incidence
log_tb <- tub %>% 
  mutate(LogInc = log(Incidence+1),
         LogGDP = log(GDP),
         LogFund = log(Funding),
         LogDen = log(PopDensity)) %>% 
  select(-Incidence, -GDP, -Funding, -PopDensity)
summary(log_tb)  
# seems okay now


### Missing Data and Multiple Imputations ####
# where is missing data?
filter(tub, is.na(Incidence))
filter(tub, is.na(Year))
filter(tub, is.na(Country))
filter(tub, is.na(Region))
# all complete

# missing GDP data:
filter(log_tb, is.na(LogGDP)) %>% 
  group_by(Country) %>% 
  count()

# missing PopDensity data:
filter(tub, is.na(PopDensity)) %>% 
  group_by(Country) %>% 
  count()

# missing health expenditure data:
filter(tub, is.na(Funding)) %>% 
  group_by(Country) %>% 
  count()

# make subset of data without NAs:
log_tb_na <- log_tb %>% 
  filter(!is.na(LogGDP)) %>% 
  filter(!is.na(LogFund)) %>% 
  filter(!is.na(LogDen)) 

# how many countries in log_tb_na?
n_distinct(log_tb_na$Country)  # 188

# impute missing data:
imp_tub <- amelia(x = as.data.frame(log_tb), m = 5, idvars = c("Country", "Region"), ts = "Year", p2s = 0)
# mean of imputations (do they make sense)
# consider Somalia as an example if the imputations are any good:
ggplot() +
  geom_line(data = imp_tub$imputations$imp1[imp_tub$imputations$imp1$Country == "Somalia", ], mapping = aes(x=Year, y=LogGDP, col = "Imputation 1")) +
  geom_line(data = imp_tub$imputations$imp2[imp_tub$imputations$imp2$Country == "Somalia", ], mapping = aes(x=Year, y=LogGDP, col = "Imputation 2")) +
  geom_line(data = imp_tub$imputations$imp3[imp_tub$imputations$imp3$Country == "Somalia", ], mapping = aes(x=Year, y=LogGDP, col = "Imputation 3")) +
  geom_line(data = imp_tub$imputations$imp4[imp_tub$imputations$imp4$Country == "Somalia", ], mapping = aes(x=Year, y=LogGDP, col = "Imputation 4")) +
  geom_line(data = imp_tub$imputations$imp5[imp_tub$imputations$imp5$Country == "Somalia", ], mapping = aes(x=Year, y=LogGDP, col = "Imputation 5")) +
  geom_line(data = log_tb[log_tb$Country == "Somalia", ], mapping = aes(x=Year, y=LogGDP, col = "Data"), col = "black")  +
  theme_bw() +
  labs(title = "Comparison of Imputed Values for the log(GDP) of Somalia", y = "Log-scaled GDP per capita in USD", col = "")
# imputations are all over the place

########by region######
# try by region:
tb_africa <- filter(log_tb, Region == "Africa") %>% 
  select(-Region)
imp_africa <- amelia(x = as.data.frame(tb_africa), m = 5,
                     idvars = "Country", ts = "Year", 
                     p2s = 0)
par(mfrow=c(2,3))
ggplot(tb_africa, aes(x=Year, y = LogGDP, col=Country)) +
  geom_line()
ggplot(tb_africa, aes(x=Year, y = LogFund, col=Country)) +
  geom_line()
ggplot(tb_africa, aes(x=Year, y = LogDen, col=Country)) +
  geom_line()

ggplot(imp_africa$imputations$imp1, aes(x=Year, y=LogGDP, col=Country)) +
  geom_line()
ggplot(imp_africa$imputations$imp1, aes(x=Year, y=LogFund, col=Country)) +
  geom_line()
ggplot(imp_africa$imputations$imp1, aes(x=Year, y=LogDen, col=Country)) +
  geom_line()
# this doesn't work
par(mfrow = c(1,1))

### Exploratory Data Analysis ####
# relationship between different variables
ggpairs(select(tub, -Country),
        lower=list(continuous="smooth"),
        diag=list(continuous = "barDiag"),
        axisLabels = "show")
ggpairs(select(log_tb, -Country),
        lower=list(continuous="smooth"),
        diag=list(continuous = "barDiag"),
        axisLabels = "show")
# variables are all correlated

# distribution of data
ggplot(tub, aes(x = Incidence)) +
  geom_histogram() +
  facet_wrap(~Region)
ggplot(tub, aes(x = PopDensity)) +
  geom_histogram() +
  facet_wrap(~Region)
ggplot(tub, aes(x = GDP)) +
  geom_histogram() +
  facet_wrap(~Region)
ggplot(tub, aes(x = Funding)) +
  geom_histogram() +
  facet_wrap(~Region)
# potential clustering based on WHO Region

# did log-transformation make sense? -> yes
par(mfrow = c(3,2))
ggplot(tub, aes(x = GDP, y = Incidence)) +
  geom_point()
ggplot(log_tb, aes(x = LogGDP, y = LogInc)) +
  geom_point()
ggplot(tub, aes(x = Funding, y = Incidence)) +
  geom_point()
ggplot(log_tb, aes(x = LogFund, y = LogInc)) +
  geom_point()
ggplot(tub, aes(x = PopDensity, y = Incidence)) +
  geom_point()
ggplot(log_tb, aes(x = LogDen, y = LogInc)) +
  geom_point()

### Spatial Visualisation #####
# get world map
world_data <- ne_countries(returnclass = "sf")

# get log tb incidence numbers for 2021
inc_data <- log_tb %>% 
  filter(Year == 2021) %>% 
  select(LogInc, Country)
# add country code to inc_data:
inc_data <- left_join(inc_data, dup_conversion, by = "Country")

# rename country code in incidence data set to aid merging
inc_data <- dplyr::rename(inc_data, adm0_a3 = Code)
# any significant countries missing because of different country code?
#View(anti_join(world_data, inc_data, by = "adm0_a3")) 
# South Sudan's country code is SDS in world data, SSD in incidence data
inc_data$adm0_a3[inc_data$adm0_a3 == "SSD"] = "SDS"

# join into one data frame:
df_world <- left_join(world_data, inc_data, by = "adm0_a3")

# Plot Spatial Visualisation of Incidence:
ggplot(data = df_world) +
  geom_sf(aes(fill = LogInc)) +
  scale_fill_viridis_c(option = "H") +
  theme_bw() +
  labs(fill = "Incident \ntuberculosis \ncases per \n100,000 population \n(log-scaled)",
       title = "Worldwide tuberculosis cases in 2021")


### Linear Mixed Effects Models #####

#### One Variable ####
# Understanding if nested models are necessary
# GDP, no random slope, only region
m_g_s <- lmer(LogInc ~ LogGDP + (1|Region), log_tb_na)
summary(m_g_s)   # significant
a <- model_performance(m_g_s)$AIC
aa <- r.squaredGLMM(m_g_s)[1,2]
ai <- performance::icc(m_g_s)[1,1]
names(aa) <- NULL

# with random slope:
m_g <- lmer(LogInc ~ LogGDP + (LogGDP|Region), log_tb_na, control = lmerControl(optimizer = "Nelder_Mead"))
summary(m_g)    # significant
c <- model_performance(m_g)$AIC
cc <- r.squaredGLMM(m_g)[1,2]
ci <- performance::icc(m_g)[1,1]
names(cc) <- NULL

# is random slope significant?
anova(m_g_s, m_g)   # random slope is significant
model_performance(m_g)    # AIC 10465

# GDP, no random slope, nested
mn_g_s <- lmer(LogInc ~ LogGDP + (1|Region) + (1|Region:Country), log_tb_na)
summary(mn_g_s) # signif
b <- model_performance(mn_g_s)$AIC
bb <- r.squaredGLMM(mn_g_s)[1,2]
bi <- performance::icc(mn_g_s)[1,1]
names(bb) <- NULL

# with random slope:
mn_g <- lmer(LogInc ~ LogGDP + (LogGDP|Region) + (LogGDP|Region:Country), log_tb_na)
summary(mn_g)   # signif
d <- model_performance(mn_g)$AIC
dd <- r.squaredGLMM(mn_g)[1,2]
di <- performance::icc(mn_g)[1,1]
names(dd) <- NULL
# random slope significant?
anova(mn_g, mn_g_s)   # significant
model_performance(mn_g)   # AIC 1220


# Funding, no random slope, only region
m_f_s <- lmer(LogInc ~ LogFund + (1|Region), log_tb_na)
summary(m_f_s)   # significant
e <- model_performance(m_f_s)$AIC
ee <- r.squaredGLMM(m_f_s)[1,2]
ei <- performance::icc(m_f_s)[1,1]
names(ee) <- NULL

# with random slope:
m_f <- lmer(LogInc ~ LogFund + (LogFund|Region), log_tb_na)
summary(m_f)    # significant
g <- model_performance(m_f)$AIC
gg <- r.squaredGLMM(m_f)[1,2]
gi <- performance::icc(m_f)[1,1]
names(gg) <- NULL
# is random slope significant?
anova(m_f_s, m_f)   # random slope is significant
model_performance(m_f)    # AIC 10508

# Funding, no random slope, nested
mn_f_s <- lmer(LogInc ~ LogFund + (1|Region) + (1|Region:Country), log_tb_na)
summary(mn_f_s) # signif
f <- model_performance(mn_f_s)$AIC
ff <- r.squaredGLMM(mn_f_s)[1,2]
fi <- performance::icc(mn_f_s)[1,1]
names(ff) <- NULL
# with random slope:
mn_f <- lmer(LogInc ~ LogFund + (LogFund|Region) + (LogFund|Region:Country), log_tb_na)
summary(mn_f)   # signif
h <- model_performance(mn_f)$AIC
hh <- r.squaredGLMM(mn_f)[1,2]
hi <- performance::icc(mn_f)[1,1]
names(hh) <- NULL
# random slope significant?
anova(mn_f, mn_f_s)   # significant
model_performance(mn_f)   # AIC 1081


# Density, no random slope, only region
m_d_s <- lmer(LogInc ~ LogDen + (1|Region), log_tb_na)
summary(m_d_s)   # significant
i <- model_performance(m_d_s)$AIC
ii <- r.squaredGLMM(m_d_s)[1,2]
ic <- performance::icc(m_d_s)[1,1]
names(ii) <- NULL

# with random slope:
m_d <- lmer(LogInc ~ LogDen + (LogDen|Region), log_tb_na)
summary(m_d)    # significant
k <- model_performance(m_d)$AIC
kk <- r.squaredGLMM(m_d)[1,2]
ki <- performance::icc(m_d)[1,1]
names(kk) <- NULL
# is random slope significant?
anova(m_d_s, m_d)   # random slope is significant
model_performance(m_d)    # AIC 12093

# Funding, no random slope, nested
mn_d_s <- lmer(LogInc ~ LogDen + (1|Region) + (1|Region:Country), log_tb_na)
summary(mn_d_s) # signif
j <- model_performance(mn_d_s)$AIC
jj <- r.squaredGLMM(mn_d_s)[1,2]
ji <- performance::icc(mn_d_s)[1,1]
names(jj) <- NULL
# with random slope:
mn_d <- lmer(LogInc ~ LogDen + (LogDen|Region) + (LogDen|Region:Country), log_tb_na, control = lmerControl(optimizer ="Nelder_Mead"))
summary(mn_d)   # model failed to converge, not signif
l <- model_performance(mn_d)$AIC
ll <- r.squaredGLMM(mn_d)[1,2]
li <- performance::icc(mn_d)[1,1]
names(ll) <- NULL
# random slope significant?
#anova(mn_f, mn_f_s)   # significant
model_performance(mn_f_s)   # AIC 2033 (no slope)

# best model: mn_g (only GDP, random slope in nested model)
# AIC 1209
r.squaredGLMM(mn_g)    # R = 0.977

m1 <- mn_g
m2 <- mn_f

#### Two Variables ####
# only test nested models, since they have better AICs compared to their non-nested counterparts

# GDP & Funding interaction, no slope
m_gf_s <- lmer(LogInc ~ LogGDP*LogFund + (1|Region) + (1|Region:Country), log_tb_na)
summary(m_gf_s)   # not significant
# slope on GDP:
m_gf_g <- lmer(LogInc ~ LogGDP*LogFund + (LogGDP|Region) + (LogGDP|Region:Country), log_tb_na, control = lmerControl(optimizer = "Nelder_Mead"))
summary(m_gf_g)   # Intercept not significant
# slope on Funding:
m_gf_f <- lmer(LogInc ~ LogGDP*LogFund + (LogFund|Region) + (LogFund|Region:Country), log_tb_na, control = lmerControl(optimizer = "Nelder_Mead"))
summary(m_gf_f)  # LogGDP & LogFund not significant

# no interaction term, no slope
m_gfi_s <- lmer(LogInc ~ LogGDP + LogFund + (1|Region) + (1|Region:Country), log_tb_na)
summary(m_gfi_s)  # significant
# random slope GDP:
m_gfi_g <- lmer(LogInc ~ LogGDP + LogFund + (LogGDP|Region) + (LogGDP|Region:Country), log_tb_na)
summary(m_gfi_g)  # LogGDP not significant
# random slope Funding
m_gfi_f <- lmer(LogInc ~ LogGDP + LogFund + (LogFund|Region) + (LogFund|Region:Country), log_tb_na, control = lmerControl(optimizer = "Nelder_Mead"))
summary(m_gfi_f)  # terms significant
# slope significant?
anova(m_gfi_s, m_gfi_f)   # random slope is significant
model_performance(m_gfi_f)   # AIC 1082  - already better AIC than single variable models
r.squaredGLMM(m_gfi_f)   # R = 0.978


# GDP & Density interaction, no slope
m_gd_s <- lmer(LogInc ~ LogGDP*LogDen + (1|Region) + (1|Region:Country), log_tb_na)
summary(m_gd_s)   # not significant
# slope on GDP:
m_gd_g <- lmer(LogInc ~ LogGDP*LogDen + (LogGDP|Region) + (LogGDP|Region:Country), log_tb_na)
summary(m_gd_g)   # not significant
# slope on Density:
m_gd_d <- lmer(LogInc ~ LogGDP*LogDen + (LogDen|Region) + (LogDen|Region:Country), log_tb_na, control = lmerControl(optimizer ="Nelder_Mead"))
summary(m_gd_d)  # not signif, model didn't converge

# no interaction term, no slope
m_gdi_s <- lmer(LogInc ~ LogGDP + LogDen + (1|Region) + (1|Region:Country), log_tb_na)
summary(m_gdi_s)  # significant
# random slope GDP:
m_gdi_g <- lmer(LogInc ~ LogGDP + LogDen + (LogGDP|Region) + (LogGDP|Region:Country), log_tb_na)
summary(m_gdi_g)  # LogGDP not significant
# random slope Density
m_gdi_d <- lmer(LogInc ~ LogGDP + LogDen + (LogDen|Region) + (LogDen|Region:Country), log_tb_na, control = lmerControl(optimizer = "Nelder_Mead"))
summary(m_gdi_d)  # LogDen not significant, singularity
# test random intercepts interaction model:
model_performance(m_gdi_s)   # AIC 2004
r.squaredGLMM(m_gdi_s)  # R 0.967


# Funding & Density interaction, no slope
m_fd_s <- lmer(LogInc ~ LogFund*LogDen + (1|Region) + (1|Region:Country), log_tb_na)
summary(m_fd_s)   # not significant
# slope on Fund:
m_fd_f <- lmer(LogInc ~ LogFund*LogDen + (LogFund|Region) + (LogFund|Region:Country), log_tb_na)
summary(m_fd_f)   # significant
# is slope significant?
anova(m_fd_s, m_fd_f)   # random slope is significant
model_performance(m_fd_f)   # AIC 909  - better

# slope on Density:
m_fd_d <- lmer(LogInc ~ LogFund*LogDen + (LogDen|Region) + (LogDen|Region:Country), log_tb_na, control = lmerControl(optimizer ="Nelder_Mead"))
summary(m_fd_d)  # not signif, model didn't converge

# no interaction term, no slope
m_fdi_s <- lmer(LogInc ~ LogFund + LogDen + (1|Region) + (1|Region:Country), log_tb_na)
summary(m_fdi_s)  # significant
model_performance(m_fdi_s)   # AIC 2020
# random slope Fund:
m_fdi_f <- lmer(LogInc ~ LogFund + LogDen + (LogFund|Region) + (LogFund|Region:Country), log_tb_na)
summary(m_fdi_f)  # LogFund not significant
# random slope Density
m_fdi_d <- lmer(LogInc ~ LogFund + LogDen + (LogDen|Region) + (LogDen|Region:Country), log_tb_na, control = lmerControl(optimizer = "Nelder_Mead"))
summary(m_fdi_d)  # LogDen not significant, singularity

# best models:
m3 <- m_gfi_s
m4 <- m_gfi_f
m5 <- m_gdi_s
m6 <- m_fd_f
m7 <- m_fdi_s

#### Three Variables #####
# linear effect of all 3 Variables, no random slope
m3_i_s <- lmer(LogInc ~ LogGDP + LogFund + LogDen + (1|Region) + (1|Region:Country), log_tb_na)
summary(m3_i_s)  # all significant
model_performance(m3_i_s)   # AIC 1943
# slope on GDP:
m3_i_g <- lmer(LogInc ~ LogGDP + LogFund + LogDen + (LogGDP|Region) + (LogGDP|Region:Country), log_tb_na)
summary(m3_i_g)  # LogGDP not significant
# slope on LogFund
m3_i_f <- lmer(LogInc ~ LogGDP + LogFund + LogDen + (LogFund|Region) + (LogFund|Region:Country), log_tb_na)
summary(m3_i_f)   # LogGDP + LogFund not significant
# slope on LogDen
m3_i_d <- lmer(LogInc ~ LogGDP + LogFund + LogDen + (LogDen|Region) + (LogDen|Region:Country), log_tb_na, control = lmerControl(optimizer ="Nelder_Mead"))
summary(m3_i_d)   # LogDen not signif, singularity



# interaction GDP & Funding, no slope
m3_gf_s <- lmer(LogInc ~ LogGDP*LogFund + LogDen + (1|Region) + (1|Region:Country), log_tb_na)
summary(m3_gf_s)   # LogGDP not signif
# slope GDP
m3_gf_g <- lmer(LogInc ~ LogGDP*LogFund + LogDen + (LogGDP|Region) + (LogGDP|Region:Country), log_tb_na)
summary(m3_gf_g)   # signif
model_performance(m3_gf_g)  # AIC 947  - better

# slope Funding
m3_gf_f <- lmer(LogInc ~ LogGDP*LogFund + LogDen + (LogFund|Region) + (LogFund|Region:Country), log_tb_na)
summary(m3_gf_f)  # signif
model_performance(m3_gf_f)   # AIC 907  - better

# slope Density
m3_gf_d <- lmer(LogInc ~ LogGDP*LogFund + LogDen + (LogDen|Region) + (LogDen|Region:Country), log_tb_na)
summary(m3_gf_f)  # signif
model_performance(m3_gf_d)  # AIC 1080
r.squaredGLMM(m3_gf_d)  # random effect too small, R2 doesn't account for RE


# interaction GDP & Density, no slope
m3_gd_s <- lmer(LogInc ~ LogGDP*LogDen + LogFund + (1|Region) + (1|Region:Country), log_tb_na)
summary(m3_gd_s)   # not signif
# slope GDP
m3_gd_g <- lmer(LogInc ~ LogGDP*LogDen + LogFund + (LogGDP|Region) + (LogGDP|Region:Country), log_tb_na)
summary(m3_gd_g)   # not signif
# slope Funding
m3_gd_f <- lmer(LogInc ~ LogGDP*LogDen + LogFund + (LogFund|Region) + (LogFund|Region:Country), log_tb_na)
summary(m3_gd_f)  # not signif
# slope Density
m3_gd_d <- lmer(LogInc ~ LogGDP*LogDen + LogFund + (LogDen|Region) + (LogDen|Region:Country), log_tb_na)
summary(m3_gd_d)  # not signif


# interaction Funding & Density, no slope
m3_fd_s <- lmer(LogInc ~ LogFund*LogDen + LogGDP + (1|Region) + (1|Region:Country), log_tb_na)
summary(m3_fd_s)   # not signif
# slope GDP
m3_fd_g <- lmer(LogInc ~ LogFund*LogDen + LogGDP + (LogGDP|Region) + (LogGDP|Region:Country), log_tb_na)
summary(m3_fd_g)   # not signif
# slope Funding
m3_fd_f <- lmer(LogInc ~ LogFund*LogDen + LogGDP + (LogFund|Region) + (LogFund|Region:Country), log_tb_na)
summary(m3_fd_f)  # not signif
# slope Density
m3_fd_d <- lmer(LogInc ~ LogFund*LogDen + LogGDP + (LogDen|Region) + (LogDen|Region:Country), log_tb_na, control = lmerControl(optimizer ="Nelder_Mead"))
summary(m3_fd_d)  # not signif



# No Interaction GDP & Funding, no slope
m3_gdfd_s <- lmer(LogInc ~ LogGDP*LogDen + LogFund*LogDen + (1|Region) + (1|Region:Country), log_tb_na)
summary(m3_gdfd_s)   # not signif
# slope GDP:
m3_gdfd_g <- lmer(LogInc ~ LogGDP*LogDen + LogFund*LogDen + (LogGDP|Region) + (LogGDP|Region:Country), log_tb_na)
summary(m3_gdfd_g)   # not signif
# slope Funding:
m3_gdfd_f <- lmer(LogInc ~ LogGDP*LogDen + LogFund*LogDen + (LogFund|Region) + (LogFund|Region:Country), log_tb_na)
summary(m3_gdfd_f)   # not signif
# slope Density
m3_gdfd_d <- lmer(LogInc ~ LogGDP*LogDen + LogFund*LogDen + (LogDen|Region) + (LogDen|Region:Country), log_tb_na)
summary(m3_gdfd_d)  # not signif


# No interaction GDP & Density, no slope
m3_gfdf_s <- lmer(LogInc ~ LogGDP*LogFund + LogFund*LogDen + (1|Region) + (1|Region:Country), log_tb_na)
summary(m3_gfdf_s)   # not signif
# slope GDP:
m3_gfdf_g <- lmer(LogInc ~ LogGDP*LogFund + LogFund*LogDen + (LogGDP|Region) + (LogGDP|Region:Country), log_tb_na)
summary(m3_gfdf_g)   # signif
model_performance(m3_gfdf_g)  # AIC 943  - good 

# slope Funding:
m3_gfdf_f <- lmer(LogInc ~ LogGDP*LogFund + LogFund*LogDen + (LogFund|Region) + (LogFund|Region:Country), log_tb_na)
summary(m3_gfdf_f)   # not signif
# slope Density
m3_gfdf_d <- lmer(LogInc ~ LogGDP*LogFund + LogFund*LogDen + (LogDen|Region) + (LogDen|Region:Country), log_tb_na)
summary(m3_gfdf_d)  # not signif


# no interaction Funding & Density
m3_gfgd_s <- lmer(LogInc ~ LogGDP*LogFund + LogGDP*LogDen + (1|Region) + (1|Region:Country), log_tb_na)
summary(m3_gfgd_s)   # not signif
# slope GDP:
m3_gfgd_g <- lmer(LogInc ~ LogGDP*LogFund + LogGDP*LogDen + (LogGDP|Region) + (LogGDP|Region:Country), log_tb_na)
summary(m3_gfgd_g)   # not signif
# slope Funding:
m3_gfgd_f <- lmer(LogInc ~ LogGDP*LogFund + LogGDP*LogDen + (LogFund|Region) + (LogFund|Region:Country), log_tb_na)
summary(m3_gfgd_f)   # not signif
# slope Density
m3_gfgd_d <- lmer(LogInc ~ LogGDP*LogFund + LogGDP*LogDen + (LogDen|Region) + (LogDen|Region:Country), log_tb_na)
summary(m3_gfgd_d)  # not signif


# Full interaction model, no slope
m3_full_s <- lmer(LogInc ~ LogGDP*LogFund*LogDen + (1|Region) + (1|Region:Country), log_tb_na)
summary(m3_full_s)  # not signif
# slope GDP
m3_full_g <- lmer(LogInc ~ LogGDP*LogFund*LogDen + (LogGDP|Region) + (LogGDP|Region:Country), log_tb_na, control = lmerControl(optimizer ="Nelder_Mead"))
summary(m3_full_g)  # not signif, model didn't converge
# slope Funding
m3_full_f <- lmer(LogInc ~ LogGDP*LogFund*LogDen + (LogFund|Region) + (LogFund|Region:Country), log_tb_na, control = lmerControl(optimizer ="Nelder_Mead"))
summary(m3_full_f)  # not signif
# slope Density
m3_full_d <- lmer(LogInc ~ LogGDP*LogFund*LogDen + (LogDen|Region) + (LogDen|Region:Country), log_tb_na, control = lmerControl(optimizer ="Nelder_Mead"))
summary(m3_full_d) # not signif


# candidate models:
m8 <- m3_i_g
m9 <- m3_gf_g
m10 <- m3_gf_f
m11 <- m3_gfdf_g


#### Tables
# one variable models
models1 <- data.frame(Predictor = c(rep("GDP",4), rep("HE",4), rep("Density",4)),
                      Type = c(rep(c("intercept", "intercept", "slope", "slope"),3)),
                      Random_Effect = c(rep(c("Region", "Region:Country"),6)),
                      Significant = c(rep("yes", 11), "no"),
                      ICC = c(ai, bi, ci, di, ei, fi, gi, hi, ic, ji, ki, li),
                      AIC = c(a,b,c,d,e,f,g,h,i,j,k,l),
                      cR = c(aa,bb,cc,dd,ee,ff,gg,hh,ii,jj,kk,ll))

table1 <- flextable(models1) %>% 
  # rename header row:
  set_header_labels(Predictor = "Predictor (Fixed Effect)",
                    Random_Effect = "Random Effect",
                    Significant = "Significant?",
                    cR = "cR²") %>% 
  merge_v(~Predictor) %>% 
  merge_v(~Type) %>% 
  fix_border_issues() %>% 
  valign(valign = 'top', j = c(1,2), part = 'body') %>%  # align Predictor column
  colformat_double(j = 6, digits = 0) %>% 
  colformat_double(j = c(5,7), digits = 3) %>% 
  hline(i = c(4, 8), border = fp_border_default()) %>%  # add horizontal line
  fontsize(size = 9, part = 'all') %>% 
  autofit()

table1

# two and three variable models:
# get AICs & R^2:
aaa <- model_performance(m3)$AIC
bbb <- model_performance(m4)$AIC
ccc <- model_performance(m5)$AIC
ddd <- model_performance(m6)$AIC
eee <- model_performance(m7)$AIC
fff <- model_performance(m8)$AIC
ggg <- model_performance(m9)$AIC
hhh <- model_performance(m10)$AIC
iii <- model_performance(m11)$AIC

aaaa <- r.squaredGLMM(m3)[1,2]
names(aaaa) <- NULL
bbbb <- r.squaredGLMM(m4)[1,2]
names(bbbb) <- NULL
cccc <- r.squaredGLMM(m5)[1,2]
names(cccc) <- NULL
dddd <- r.squaredGLMM(m6)[1,2]
names(dddd) <- NULL
eeee <- r.squaredGLMM(m7)[1,2]
names(eeee) <- NULL
ffff <- r.squaredGLMM(m8)[1,2]
names(ffff) <- NULL
gggg <- r.squaredGLMM(m9)[1,2]
names(gggg) <- NULL
hhhh <- r.squaredGLMM(m10)[1,2]
names(hhhh) <- NULL
iiii <- r.squaredGLMM(m11)[1,2]
names(iiii) <- NULL

# get ICCs:
aii <- performance::icc(m3)[1,1]
bii <- performance::icc(m4)[1,1]
cii <- performance::icc(m5)[1,1]
dii <- performance::icc(m6)[1,1]
eii <- performance::icc(m7)[1,1]
fii <- performance::icc(m8)[1,1]
gii <- performance::icc(m9)[1,1]
hii <- performance::icc(m10)[1,1]
iic <- performance::icc(m11)[1,1]

models2 <- data.frame(Name = c("m3", "m4", "m5", "m6", "m7", "m8", "m9", "m10", "m11"),
                      Predictors = c("GDP, HE", "GDP, HE", "GDP, Density",
                                     "HE, Density, HE:Density", "HE, Density",
                                     "GDP, HE, Density", 
                                     "GDP, HE, Density, GDP:HE",
                                     "GDP, HE, Density, GDP:HE", 
                                     "GDP, HE, Density, GDP:HE, HE:Density"),
                      RE = c("intercept only", "slope (HE)", "intercept only",
                             "slope (HE)", "intercept only", "intercept only",
                             "slope (GDP)", "slope (HE)", "slope (GDP)"),
                      ICC = c(aii, bii, cii, dii, eii, fii, gii, hii, iic),
                      AIC = c(aaa,bbb,ccc,ddd,eee,fff,ggg,hhh,iii),
                      cR = c(aaaa,bbbb,cccc,dddd,eeee,ffff,gggg,hhhh,iiii))


table2 <- flextable(models2) %>% 
  # rename header row:
  set_header_labels(Name = "Model",
                    Predictors = "Fixed Effects",
                    RE = "Random Effect type",
                    cR = "cR²") %>% 
  fix_border_issues() %>% 
  colformat_double(j = 5, digits = 0) %>% 
  colformat_double(j = c(4,6), digits = 3) %>% 
  #  hline(i = c(4, 8), border = fp_border_default()) %>%  # add horizontal line
  fontsize(size = 9, part = 'all') %>% 
  autofit()
table2

### Assessing Model Fit #####
# Q-Q-Plots
par(mfrow = c(2,3))
qqnorm(residuals(m1), xlab = "Theoretical Quantiles (m1)")
qqnorm(residuals(m2), xlab = "Theoretical Quantiles (m2)")
qqnorm(residuals(m6), xlab = "Theoretical Quantiles (m6)")
qqnorm(residuals(m9), xlab = "Theoretical Quantiles (m9)")
qqnorm(residuals(m10), xlab = "Theoretical Quantiles (m10)")
qqnorm(residuals(m11), xlab = "Theoretical Quantiles (m11)")
par(mfrow = c(1,1))

# check lowest level residuals:
# make data frame containing all residuals, separated by model
incdiag <- data.frame(Residuals = c(resid(m1), resid(m2), resid(m3), 
                                    resid(m4), resid(m5), resid(m6)),
                      Region = rep(log_tb_na$Region, 6),
                      Fitted = c(fitted(m1), fitted(m2), fitted(m3),
                                 fitted(m4), fitted(m5), fitted(m6)),
                      Model = c(rep("m01", length(resid(m1))), 
                                rep("m02", length(resid(m2))), 
                                rep("m06", length(resid(m6))),
                                rep("m09", length(resid(m9))), 
                                rep("m10", length(resid(m10))), 
                                rep("m11", length(resid(m11))))
)
ggplot(data=incdiag, aes(x=Fitted, y=Residuals, col=Model)) +
  geom_point(alpha = 0.5, size = 0.5) +
  facet_wrap(~Region) +
  theme_bw()

### K-Fold Cross Validation ####
set.seed(1234)
# make 10 partitions:
folds <- log_tb_na$LogInc %>% 
  createFolds()
# make partitions 
p1 <- log_tb_na[folds$Fold01, ]
p2 <- log_tb_na[folds$Fold02, ]
p3 <- log_tb_na[folds$Fold03, ]
p4 <- log_tb_na[folds$Fold04, ]
p5 <- log_tb_na[folds$Fold05, ]
p6 <- log_tb_na[folds$Fold06, ]
p7 <- log_tb_na[folds$Fold07, ]
p8 <- log_tb_na[folds$Fold08, ]
p9 <- log_tb_na[folds$Fold09, ]
p10 <- log_tb_na[folds$Fold10, ]

# rest of data:
p_1 <- log_tb_na[-folds$Fold01, ]
p_2 <- log_tb_na[-folds$Fold02, ]
p_3 <- log_tb_na[-folds$Fold03, ]
p_4 <- log_tb_na[-folds$Fold04, ]
p_5 <- log_tb_na[-folds$Fold05, ]
p_6 <- log_tb_na[-folds$Fold06, ]
p_7 <- log_tb_na[-folds$Fold07, ]
p_8 <- log_tb_na[-folds$Fold08, ]
p_9 <- log_tb_na[-folds$Fold09, ]
p_10 <- log_tb_na[-folds$Fold10, ]

# model m1: LogInc ~ LogGDP + (LogGDP | Region) + (LogGDP | Region:Country)
err1 <- c()
t_err1 <- c()
for (n in 1:10){
  train.inc <- NULL
  test.inc <- NULL
  if (n == 1){
    # first partition is test data
    train.inc <- p_1
    test.inc <- p1
  } else if (n == 2){
    train.inc <- p_2
    test.inc <- p2
  } else if (n == 3){
    train.inc <- p_3
    test.inc <- p3
  } else if (n == 4){
    train.inc <- p_4
    test.inc <- p4
  } else if (n == 5){
    train.inc <- p_5
    test.inc <- p5
  } else if (n == 6){
    train.inc <- p_6
    test.inc <- p6
  } else if (n == 7){
    train.inc <- p_7
    test.inc <- p7
  } else if (n == 8){
    train.inc <- p_8
    test.inc <- p8
  } else if (n == 9){
    train.inc <- p_9
    test.inc <- p9
  } else if (n == 10){
    train.inc <- p_10
    test.inc <- p10
  }
  
  model <- lmer(LogInc ~ LogGDP + (LogGDP | Region) + (LogGDP | Region:Country), 
                train.inc, control = lmerControl(optimizer = "Nelder_Mead"))
  pred <- model %>% predict(test.inc)
  pred_train <- model %>% predict(train.inc)
  names(pred) <- NULL
  # prediction error (root MSE normalised): - validation
  pred_err <- RMSE(pred, test.inc$LogInc)/mean(test.inc$LogInc)
  train_err <- RMSE(pred_train, train.inc$LogInc)/mean(train.inc$LogInc)
  err1 <- cbind(err1, pred_err)
  t_err1 <- cbind(t_err1, train_err)
}

# model m2: LogInc ~ LogFund + (LogFund | Region) + (LogFund | Region:Country)
err2 <- c()
t_err2 <- c()
for (n in 1:10){
  train.inc <- NULL
  test.inc <- NULL
  if (n == 1){
    # first partition is test data
    train.inc <- p_1
    test.inc <- p1
  } else if (n == 2){
    train.inc <- p_2
    test.inc <- p2
  } else if (n == 3){
    train.inc <- p_3
    test.inc <- p3
  } else if (n == 4){
    train.inc <- p_4
    test.inc <- p4
  } else if (n == 5){
    train.inc <- p_5
    test.inc <- p5
  } else if (n == 6){
    train.inc <- p_6
    test.inc <- p6
  } else if (n == 7){
    train.inc <- p_7
    test.inc <- p7
  } else if (n == 8){
    train.inc <- p_8
    test.inc <- p8
  } else if (n == 9){
    train.inc <- p_9
    test.inc <- p9
  } else if (n == 10){
    train.inc <- p_10
    test.inc <- p10
  }
  
  model <- lmer(LogInc ~ LogFund + (LogFund | Region) + (LogFund | Region:Country), 
                train.inc, control = lmerControl(optimizer = "Nelder_Mead"))
  pred <- model %>% predict(test.inc)
  pred_train <- model %>% predict(train.inc)
  names(pred) <- NULL
  # prediction error (root MSE normalised): - validation
  pred_err <- RMSE(pred, test.inc$LogInc)/mean(test.inc$LogInc)
  train_err <- RMSE(pred_train, train.inc$LogInc)/mean(train.inc$LogInc)
  err2 <- cbind(err2, pred_err)
  t_err2 <- cbind(t_err2, train_err)
}

# model m6: LogInc ~ LogFund * LogDen + (LogFund | Region) + (LogFund | Region:Country)
err6 <- c()
t_err6 <- c()
for (n in 1:10){
  train.inc <- NULL
  test.inc <- NULL
  if (n == 1){
    # first partition is test data
    train.inc <- p_1
    test.inc <- p1
  } else if (n == 2){
    train.inc <- p_2
    test.inc <- p2
  } else if (n == 3){
    train.inc <- p_3
    test.inc <- p3
  } else if (n == 4){
    train.inc <- p_4
    test.inc <- p4
  } else if (n == 5){
    train.inc <- p_5
    test.inc <- p5
  } else if (n == 6){
    train.inc <- p_6
    test.inc <- p6
  } else if (n == 7){
    train.inc <- p_7
    test.inc <- p7
  } else if (n == 8){
    train.inc <- p_8
    test.inc <- p8
  } else if (n == 9){
    train.inc <- p_9
    test.inc <- p9
  } else if (n == 10){
    train.inc <- p_10
    test.inc <- p10
  }
  
  model <- lmer(LogInc ~ LogFund * LogDen + (LogFund | Region) + (LogFund | Region:Country), 
                train.inc, control = lmerControl(optimizer = "Nelder_Mead"))
  pred <- model %>% predict(test.inc)
  pred_train <- model %>% predict(train.inc)
  names(pred) <- NULL
  # prediction error (root MSE normalised): - validation
  pred_err <- RMSE(pred, test.inc$LogInc)/mean(test.inc$LogInc)
  train_err <- RMSE(pred_train, train.inc$LogInc)/mean(train.inc$LogInc)
  err6 <- cbind(err6, pred_err)
  t_err6 <- cbind(t_err6, train_err)
}

# model m9: LogInc ~ LogGDP * LogFund + LogDen + (LogGDP | Region) + (LogGDP|Region:Country)
err9 <- c()
t_err9 <- c()
for (n in 1:10){
  train.inc <- NULL
  test.inc <- NULL
  if (n == 1){
    # first partition is test data
    train.inc <- p_1
    test.inc <- p1
  } else if (n == 2){
    train.inc <- p_2
    test.inc <- p2
  } else if (n == 3){
    train.inc <- p_3
    test.inc <- p3
  } else if (n == 4){
    train.inc <- p_4
    test.inc <- p4
  } else if (n == 5){
    train.inc <- p_5
    test.inc <- p5
  } else if (n == 6){
    train.inc <- p_6
    test.inc <- p6
  } else if (n == 7){
    train.inc <- p_7
    test.inc <- p7
  } else if (n == 8){
    train.inc <- p_8
    test.inc <- p8
  } else if (n == 9){
    train.inc <- p_9
    test.inc <- p9
  } else if (n == 10){
    train.inc <- p_10
    test.inc <- p10
  }
  
  model <- lmer(LogInc ~ LogGDP * LogFund + LogDen + (LogGDP | Region) + (LogGDP|Region:Country), 
                train.inc, control = lmerControl(optimizer = "Nelder_Mead"))
  pred <- model %>% predict(test.inc)
  pred_train <- model %>% predict(train.inc)
  names(pred) <- NULL
  # prediction error (root MSE normalised): - validation
  pred_err <- RMSE(pred, test.inc$LogInc)/mean(test.inc$LogInc)
  train_err <- RMSE(pred_train, train.inc$LogInc)/mean(train.inc$LogInc)
  err9 <- cbind(err9, pred_err)
  t_err9 <- cbind(t_err9, train_err)
}

# model m10: LogInc ~ LogGDP * LogFund + LogDen + (LogFund | Region) + (LogFund |      Region:Country)
err10 <- c()
t_err10 <- c()
for (n in 1:10){
  train.inc <- NULL
  test.inc <- NULL
  if (n == 1){
    # first partition is test data
    train.inc <- p_1
    test.inc <- p1
  } else if (n == 2){
    train.inc <- p_2
    test.inc <- p2
  } else if (n == 3){
    train.inc <- p_3
    test.inc <- p3
  } else if (n == 4){
    train.inc <- p_4
    test.inc <- p4
  } else if (n == 5){
    train.inc <- p_5
    test.inc <- p5
  } else if (n == 6){
    train.inc <- p_6
    test.inc <- p6
  } else if (n == 7){
    train.inc <- p_7
    test.inc <- p7
  } else if (n == 8){
    train.inc <- p_8
    test.inc <- p8
  } else if (n == 9){
    train.inc <- p_9
    test.inc <- p9
  } else if (n == 10){
    train.inc <- p_10
    test.inc <- p10
  }
  
  model <- lmer(LogInc ~ LogGDP * LogFund + LogDen + (LogFund | Region) + (LogFund |      Region:Country), 
                train.inc, control = lmerControl(optimizer = "Nelder_Mead"))
  pred <- model %>% predict(test.inc)
  pred_train <- model %>% predict(train.inc)
  names(pred) <- NULL
  # prediction error (root MSE normalised): - validation
  pred_err <- RMSE(pred, test.inc$LogInc)/mean(test.inc$LogInc)
  train_err <- RMSE(pred_train, train.inc$LogInc)/mean(train.inc$LogInc)
  err10 <- cbind(err10, pred_err)
  t_err10 <- cbind(t_err10, train_err)
}

# model m11: LogInc ~ LogGDP * LogFund + LogFund * LogDen + (LogGDP | Region) +  (LogGDP | Region:Country)
err11 <- c()
t_err11 <- c()
for (n in 1:10){
  train.inc <- NULL
  test.inc <- NULL
  if (n == 1){
    # first partition is test data
    train.inc <- p_1
    test.inc <- p1
  } else if (n == 2){
    train.inc <- p_2
    test.inc <- p2
  } else if (n == 3){
    train.inc <- p_3
    test.inc <- p3
  } else if (n == 4){
    train.inc <- p_4
    test.inc <- p4
  } else if (n == 5){
    train.inc <- p_5
    test.inc <- p5
  } else if (n == 6){
    train.inc <- p_6
    test.inc <- p6
  } else if (n == 7){
    train.inc <- p_7
    test.inc <- p7
  } else if (n == 8){
    train.inc <- p_8
    test.inc <- p8
  } else if (n == 9){
    train.inc <- p_9
    test.inc <- p9
  } else if (n == 10){
    train.inc <- p_10
    test.inc <- p10
  }
  
  model <- lmer(LogInc ~ LogGDP * LogFund + LogFund * LogDen + (LogGDP | Region) +  (LogGDP | Region:Country), 
                train.inc, control = lmerControl(optimizer = "Nelder_Mead"))
  pred <- model %>% predict(test.inc)
  pred_train <- model %>% predict(train.inc)
  names(pred) <- NULL
  # prediction error (root MSE normalised): - validation
  pred_err <- RMSE(pred, test.inc$LogInc)/mean(test.inc$LogInc)
  train_err <- RMSE(pred_train, train.inc$LogInc)/mean(train.inc$LogInc)
  err11 <- cbind(err11, pred_err)
  t_err11 <- cbind(t_err11, train_err)
}

# display results in a table:
ktable <- data.frame(Model =c("m1", "m2", "m6", "m9", "m10", "m11"),
                     train = c(mean(t_err1), mean(t_err2), mean(t_err6), 
                               mean(t_err9), mean(t_err10), mean(t_err11)),
                     val = c(mean(err1), mean(err2), mean(err6), mean(err9), 
                             mean(err10), mean(err11))
)

table3 <- flextable(ktable) %>% 
  set_header_labels(train = "Training Error",
                    val = "Validation Error") %>% 
  colformat_double(j = c(2,3), digits = 3) %>% 
  fontsize(size = 10, part = 'all') %>% 
  autofit()

table3
# model m6 is best

### Final Model Parameters ####
summary(m6)
r.squaredGLMM(m6)
performance::icc(m6)
model_performance(m6)  # AIC = 909.496

### Model Predictions ####
mp1 <- plot_model(m6, type = "pred", terms = c("LogFund", "LogDen"), title = "") +
  geom_point() +
  labs(x = "Log-scaled \nhealth expenditure \nper capita in USD", y = "Log-scaled incident TB cases per 100,000 population", col = "Log-scaled \npopulation \ndensity \nper square km") +
  theme_bw()

mp2 <- plot_model(m6, type = "pred", terms = c("LogDen", "LogFund"), title = "") +
  geom_point() +
  labs(x = "Log-scaled \npopulation density \nper square km", y = "", col = "Log-scaled \nhealth expenditure \nper capita \nin USD") +
  theme_bw()

multiplot(mp1, mp2, cols = 2)

### Random effects #####
# random effect of Region
pm1 <- plot_model(m6, type = "re", title = "Random effects of WHO region")[[2]]
pm1$data$facet <- factor(pm1$data$facet)
levels(pm1$data$facet) <- c("log(Health expenditure)", "Region (Intercept)")
pm1 <- pm1 + theme_bw()
pm1

pm2 <- plot_model(m6, type = "re", title = "Random effects of country")[[1]] 
pm2$data$facet <- factor(pm2$data$facet)
levels(pm2$data$facet) <- c("log(Health expenditure)", "Country (Intercept)")
pm2 <- pm2 + theme_bw()
pm2

### Marginal effects
plot_model(m6, type="eff", terms = "LogFund") +
  geom_point(aes(x=LogFund, y = LogInc), data = log_tb_na)
plot_model(m6, type = "eff", terms = "LogDen") +
  geom_point(aes(x=LogDen, y = LogInc), data = log_tb_na)
# these plots aren't too informative, keep out of report