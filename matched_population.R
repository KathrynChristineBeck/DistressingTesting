### create matched study population ######

# load packages
library(tidyverse)
library(arsenal)
library(ggsci)
library(Matching)
library(MatchIt)
library(cobalt)
library(fastDummies)
select <- dplyr::select

# load study df of (alive) students in 3rd year of upper secondary

## load data for immigration background groupings

land<-land %>% 
  select(ID, landbak3gen_todeling)

population<-left_join(population, land, by="ID")

#create groupings of invkat
population <- population %>% 
  mutate(imcat = case_when((invkat=="A")|(invkat=="G")~"AG",
                           (invkat=="B")~"B",
                           (invkat=="C")~"C",
                           (invkat=="E")|(invkat=="F")~"EF",
                           TRUE~NA_character_))

population$immback<-paste(population$imcat, population$landbak3gen_todeling, sep = "")


## create a dummy variable for treat (grade 1) and control (grade 2)

population <- population %>% 
  mutate(TREAT_SKR =  case_when((NOR_SKR == 1) ~ 1,
                                       (NOR_SKR == 2) ~ 0))


## Subset for written Norwegian hovedmål exam
# exclude those with missing STP and gsp

population <- population %>%
  filter((!is.na(NOR_SKR)))  %>%
  # exclude those taking a transfer year from vocational track using NUS codes to select specific course
  filter((FAGKODE=="NOR1211")|(FAGKODE=="VG4000")) %>%  
  # exclude those missing teacher grade
  filter(!is.na(STP)) %>%
  # exclude those missing GPA
  filter(!is.na(gspstd))

#select those in academic track
population<-population %>% 
  filter((STUDRETN == 62) |
           (STUDRETN == 61) | (STUDRETN == 63) |
           (STUDRETN == 60) | (STUDRETN == 21) |
           (STUDRETN == 22) | (STUDRETN == 23)) %>%
  # subset to study years (2005/2006 - 2017/2018)
  filter((SKOLEAR!="20032004")&(SKOLEAR!="20042005")) 

# duplicates check
n_distinct(population$ID)

#drop few duplicates
population <- population %>%
  dplyr::distinct(ID, .keep_all = T)

# create exam year variable
population$examyear <- substr(population$SKOLEAR, 5, 8)
population$examyear<- as.numeric(as.character(population$examyear))

# create dummy variable for gender
population<-population %>% 
  mutate(gender = case_when((kjoenn == 1) ~ 0, (kjoenn == 2) ~ 1))


# create age (in years) variable
population<- population %>% 
  mutate(age = examyear-foedselsaar) %>% 
  # exclude those 22 and older at exam year
  filter(age < 22)

# this leaves entire sample with grades from 1-6

#keep those with a 1 or 2
population<- population %>% 
  filter(!is.na(TREAT_SKR)) 

# replace NA for mother and fathers education with missing value to include separate missing category in matching
population <- population %>% 
  mutate(mor_edu = replace_na(mor_edu, 9)) %>%
  mutate(far_edu = replace_na(far_edu, 9))

# select variables
unmatched_df <- population %>% 
  select(
    ID,
    mor_edu,
    far_edu,
    STP,
    gspstd,
    STUDRETN,
    gender,
    TREAT_SKR,
    examyear,
    age,
    immback
  )


#create dummy variables from categorical vars for PS estimation
unmatched_df <- dummy_cols(unmatched_df, select_columns = "examyear")
unmatched_df <- dummy_cols(unmatched_df, select_columns = "STP")
unmatched_df <- dummy_cols(unmatched_df, select_columns = "mor_edu")
unmatched_df <- dummy_cols(unmatched_df, select_columns = "far_edu")
unmatched_df <- dummy_cols(unmatched_df, select_columns = "immback")


## function of names for covariate balance check
v<-data.frame(old=c("gender", "examyear_2006", "examyear_2007", "examyear_2008", "examyear_2009", "examyear_2010", "examyear_2011",
                    "examyear_2012", "examyear_2013", "examyear_2014", "examyear_2015", "examyear_2016", "examyear_2017", "examyear_2018",
                    "STP_1", "STP_2", "STP_3", "STP_4", "STP_5", "STP_6", "mor_edu_1", "mor_edu_2", "mor_edu_3", "mor_edu_4", 
                    "mor_edu_9", "far_edu_1", "far_edu_2", "far_edu_3", "far_edu_4", "far_edu_9", "poly(gspstd, 3)_1", "poly(gspstd, 3)_2",
                    "poly(gspstd, 3)_3", "age","immback_AG0","immback_AG1","immback_AG2","immback_B1","immback_B2","immback_C1","immback_C2",
                    "immback_EF1","immback_EF2"),
              new=c("Gender", "2006", "2007", "2008", "2009", "2010", "2011","2012", "2013", "2014", "2015", "2016", "2017", "2018",
                    "Teacher Grade: 1", "Teacher Grade: 2", "Teacher Grade: 3", "Teacher Grade: 4", "Teacher Grade: 5", "Teacher Grade: 6", 
                    "Mother's Edu: Lower secondary", "Mother's Edu: Upper secondary", 
                    "Mother's Edu: Bachelor", "Mother's Edu: Master", "Mother's Edu: Missing", "Father's Edu: Lower secondary", 
                    "Father's Edu: Upper secondary", "Father's Edu: Bachelor", 
                    "Father's Edu: Master", "Father's Edu: Missing","Standardized GPA", "GPA Polynomial: 2nd","GPA Polynomial: 3rd", "Age",
                    "2 Norwegian Parents: Norway","2 Norwegian Parents: Western", "2 Norwegian Parents: Non-Western",
                    "Immigrant: Western","Immigrant: Non-Western", "2 Immigrant parents: Western", "2 Immigrant parents: Non-Western",
                    "1 Norwegian Parent: Western", "1 Norwegian Parent: Non-Western"))

#calculate propensity scores#
psmodel <- glm(TREAT_SKR~gender+examyear_2006+examyear_2007+ examyear_2008+examyear_2009+
                     examyear_2010+ examyear_2011+ examyear_2012+ examyear_2013+ 
                     examyear_2014+ examyear_2015+ examyear_2016 + examyear_2017 +
                     examyear_2018+ STP_1+ STP_2 + STP_3 + STP_4 + STP_5 + STP_6 +
                    mor_edu_1+ mor_edu_2+ mor_edu_3 + mor_edu_4 + 
                     mor_edu_9+ far_edu_1 + far_edu_2 + far_edu_3 + far_edu_4 +
                     far_edu_9+ poly(gspstd,3) + age + immback_AG0 + immback_AG1 + immback_AG2 +
                     immback_B1 + immback_B2 + immback_C1 + immback_C2 + immback_EF1 + immback_EF2, 
                   data = unmatched_df, family = binomial())
summary(psmodel)
unmatched_df$p <- predict(psmodel, newdata = unmatched_df, type = "response")
unmatched_df$ps <- log((1 - unmatched_df$p) / (unmatched_df$p))



##### Nearest Neighbor matching with 0.25*SD caliper, without replacement
set.seed(100)
match_out <-
  matchit(TREAT_SKR~gender+examyear_2006+examyear_2007+ examyear_2008+examyear_2009+
            examyear_2010+ examyear_2011+ examyear_2012+ examyear_2013+ 
            examyear_2014+ examyear_2015+ examyear_2016 + examyear_2017 +
            examyear_2018+ STP_1+ STP_2 + STP_3 + STP_4 + STP_5 + STP_6 +
            mor_edu_1+ mor_edu_2+ mor_edu_3 + mor_edu_4 + 
            mor_edu_9+ far_edu_1 + far_edu_2 + far_edu_3 + far_edu_4 +
            far_edu_9+ poly(gspstd,3)+ age + immback_AG0 + immback_AG1 + immback_AG2 +
            immback_B1 + immback_B2 + immback_C1 + immback_C2 + immback_EF1 + immback_EF2,
          data = unmatched_df,
          method = "nearest",
          distance = unmatched_df$ps,
          m.order = "random",
          caliper = 0.25,
          replace = FALSE
  )

# summary
match_out
# create matched dataframe
match_df<-MatchIt::match.data(match_out)





## Boxplot figure of common support region using linear propensity scores##

#unmatched
unmatched_df %>%
  mutate(t = factor(TREAT_SKR, labels = c("Passed", "Failed"))) %>%
  ggplot(aes(
    x = t,
    y = ps,
    color = t,
    fill = t
  )) +
  theme_classic() +
  geom_boxplot(alpha = 0.7) +
  labs(x = "Exam Result", y = "Predicted Probability", 
       title = "Estimated Linear Propensity Scores",
       subtitle = "Unmatched Sample") +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5, face ="bold"),
        plot.subtitle = element_text(hjust = 0.5, size = 10)) +
  scale_color_jama() +
  scale_fill_jama()

#matched
match_df %>%
  mutate(t = factor(TREAT_SKR, labels = c("Passed", "Failed"))) %>%
  ggplot(aes(
    x = t,
    y = ps,
    color = t,
    fill = t
  )) +
  theme_classic() +
  geom_boxplot(alpha = 0.7) +
  labs(x = "Exam Result", y = "Predicted Probability", 
       title = "Estimated Linear Propensity Scores",
       subtitle = "Matched Sample") +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5, face ="bold"),
        plot.subtitle = element_text(hjust = 0.5, size = 10)) +
  scale_color_jama() +
  scale_fill_jama()


# love plot for covariate balance

# SMD
love.plot(
  match_out,
  binary = "std",
  stats = c("mean.diffs"),
  thresholds = c(.1),
  line = TRUE,
  var.order = "unadjusted",
  drop.distance=TRUE,
  var.names=v,
  sample.names = c("Unmatched", "Matched")) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        legend.title = element_blank(),
        plot.subtitle = element_text(hjust = 0.5, size = 10)) +
  scale_color_jama()

# Variance ratio
love.plot(
  match_out,
  binary = "std",
  stats = c("variance.ratios"),
  thresholds = c(2),
  line = TRUE,
  var.order = "unadjusted",
  drop.distance=TRUE,
  var.names=v,
  sample.names = c("Unmatched", "Matched")) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        legend.title = element_blank(),
        plot.subtitle = element_text(hjust = 0.5, size = 10)) +
  scale_color_jama()



# balance table
bal.tab(match_out, m.threshold=0.1)



