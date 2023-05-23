################## ALTERNATIVE PSM TECHNIQUES ###############################
library(tidyverse)
library(arsenal)
library(data.table)
library(ggsci)
library(Matching)
library(MatchIt)
library(cobalt)
library(fastDummies)
library(optmatch)
library(jtools)
library(WeightIt)
library(MASS)
select <- dplyr::select

# load sample df
# for analyses, use same code from regression_analyses file with new matched df


################ NN with replacement 1:2  ----
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

## function for covariate balance check
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


##### Nearest Neighbor matching with 0.25*SD caliper, with replacement (2 control per treat)
# and covariate balance

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
          replace = TRUE,
          ratio = 2
  )

match_out
match_df<-MatchIt::match.data(match_out)

# covariate balance
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



########################### NN with replacement 1:3 ----
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

## function for covariate balance check
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


##### Nearest Neighbor matching with 0.25*SD caliper, with replacement (3 control per treat)
# and covariate balance

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
          replace = TRUE,
          ratio = 3
  )

match_out
match_df<-MatchIt::match.data(match_out)

# covariate balance
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



################################# NN  without replacement 1:2 ----
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

## function for covariate balance check
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


##### NN  matching with 0.25*SD caliper, without replacement, 1:2
# and covariate balance

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
          replace = F,
          ratio = 2
  )

match_out
match_df<-MatchIt::match.data(match_out)

# covariate balance
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



################################# NN without replacement  1:3 ----
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

## function for covariate balance check
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


##### NN  matching with 0.25*SD caliper, without replacement, 1:3
# and covariate balance

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
          replace = F,
          ratio = 3
  )

match_out
match_df<-MatchIt::match.data(match_out)

# covariate balance
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



############# Propensity score weighing using ATT weights ----

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

## function for covariate balance check
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

d.weights <- unmatched_df %>%
  mutate(ate_w = ifelse(TREAT_SKR == 0, 1/(1-p), 1/p),
         att_w = ifelse(TREAT_SKR == 0, p/(1-p), 1))



# ATT weight
weight.out <- WeightIt::weightit(TREAT_SKR~gender+examyear_2006+examyear_2007+ examyear_2008+examyear_2009+
                                   examyear_2010+ examyear_2011+ examyear_2012+ examyear_2013+ 
                                   examyear_2014+ examyear_2015+ examyear_2016 + examyear_2017 +
                                   examyear_2018+ STP_1+ STP_2 + STP_3 + STP_4 + STP_5 + STP_6 +
                                   mor_edu_1+ mor_edu_2+ mor_edu_3 + mor_edu_4 + 
                                   mor_edu_9+ far_edu_1 + far_edu_2 + far_edu_3 + far_edu_4 +
                                   far_edu_9+ poly(gspstd,3) + age + immback_AG0 + immback_AG1 + immback_AG2 +
                                   immback_B1 + immback_B2 + immback_C1 + immback_C2 + immback_EF1 + immback_EF2,
                                 data = unmatched_df,
                                 ps = unmatched_df$p,
                                 estimand = "ATT")
summary(weight.out)

# Visual balance checks may be carried out with:
cobalt::bal.tab(weight.out, m.threshold = .05, disp.v.ratio = TRUE)


## for regression analyses use the following code after finding psychological diagnoses with "regression_analysis" file to use ATT weights

#Negative Binomial
nb <-
  glm.nb(
    num ~ TREAT_SKR + gender + STP + examyear + mor_edu + far_edu + poly(gspstd, 3) + age + immback,
    data = unmatched_df,
    link = log,
    weights = weight.out$weights
  )

summ(
  nb,
  confint = TRUE,
  ci.width = 0.95,
  robust = TRUE,
  exp = T
)


# logistic regression
log <-
  glm(
    zero ~ TREAT_SKR + gender + STP + examyear + mor_edu + far_edu + poly(gspstd, 3) + age + immback,
    data = unmatched_df,
    family = binomial(link = "logit"),
    weights =  weight.out$weights
  )

summ(
  log,
  confint = TRUE,
  ci.width = 0.95,
  robust = TRUE,
  exp = T
)



####### Match with additional vars #############

# create lowest income quintile indicator
population = population %>%  
  mutate(low_inc_indicator = case_when(income_quintile==1~1, TRUE~0))

# create categorical variables for class size and number of children in household
population = population %>% 
  mutate(class_size_cat = cut(n, breaks = c(0, 24,50,100,Inf),
                              labels = c("0_24", "25_49", "50_99", "100"))) %>% 
  mutate(children_cat = cut(num_children, breaks=c( 0.9, 1.9, 2.9, 3.9, Inf),
                            labels= c( "1","2","3","4")))

# giving those with no peer gpa average value (standaridized, so mean = 0)
population <- population %>% 
  mutate(jackknife_gpa_peers = replace_na(jackknife_gpa_peers, 0))

# replace missing centrality index with lowest centrality (equal to centrality index = 0)
population <- population %>% 
  mutate(index = replace_na(index, 6))


population$class_size_cat <- as.character(population$class_size_cat)
population$children_cat <- as.character(population$children_cat)


# replace na values with "9" for number of children in household and "999" for class size - will be separate "missing" category
population <- population %>% 
  mutate(children_cat = case_when(is.na(children_cat)~"9", TRUE~children_cat)) %>% 
  mutate(class_size_cat = case_when(is.na(class_size_cat)~"999", TRUE~class_size_cat))



# select vars
unmatched_df <-
  population %>% 
  select(
    ID,
    SKOLEAR,
    mor_edu,
    far_edu,
    STP,
    gspstd,
    gender,
    TREAT_SKR,
    examyear,
    age,
    immback,
    centrality,
    region,
    index,
    jackknife_gpa_peers,
    low_inc_indicator,
    class_size_cat,
    children_cat
  )

#create dummy variables from categorical vars for PS estimation
unmatched_df <- dummy_cols(unmatched_df, select_columns = "examyear")
unmatched_df <- dummy_cols(unmatched_df, select_columns = "STP")
unmatched_df <- dummy_cols(unmatched_df, select_columns = "mor_edu")
unmatched_df <- dummy_cols(unmatched_df, select_columns = "far_edu")
unmatched_df <- dummy_cols(unmatched_df, select_columns = "immback")
unmatched_df <- dummy_cols(unmatched_df, select_columns = "region")
unmatched_df <- dummy_cols(unmatched_df, select_columns = "children_cat")
unmatched_df <- dummy_cols(unmatched_df, select_columns = "class_size_cat")

## function for covariate balance check
v<-data.frame(old=c("gender", "examyear_2006", "examyear_2007", "examyear_2008", "examyear_2009", "examyear_2010", "examyear_2011",
                    "examyear_2012", "examyear_2013", "examyear_2014", "examyear_2015", "examyear_2016", "examyear_2017", "examyear_2018",
                    "STP_1", "STP_2", "STP_3", "STP_4", "STP_5", "STP_6", "mor_edu_1", "mor_edu_2", "mor_edu_3", "mor_edu_4", 
                    "mor_edu_9", "far_edu_1", "far_edu_2", "far_edu_3", "far_edu_4", "far_edu_9","poly(gspstd, 3)_1", "poly(gspstd, 3)_2","poly(gspstd, 3)_3", "age",
                    "immback_AG0","immback_AG1","immback_AG2","immback_B1","immback_B2","immback_C1","immback_C2","immback_EF1","immback_EF2","index",
                    "region_vest", "region_ost","region_nord","region_sor","region_trondelag", "children_cat_1","children_cat_2","children_cat_3",
                    "children_cat_4", "children_cat_9", "class_size_cat_0_24", "class_size_cat_25_49", "class_size_cat_50_99", "class_size_cat_100", "class_size_cat_999", "low_inc_indicator", "jackknife_gpa_peers"),
              new=c("Gender", "2006", "2007", "2008", "2009", "2010", "2011","2012", "2013", "2014", "2015", "2016", "2017", "2018",
                    "Teacher Grade: 1", "Teacher Grade: 2", "Teacher Grade: 3", "Teacher Grade: 4", "Teacher Grade: 5", "Teacher Grade: 6", 
                     "Mother's Edu: Lower secondary", "Mother's Edu: Upper secondary", 
                    "Mother's Edu: Bachelor", "Mother's Edu: Master", "Mother's Edu: Missing", "Father's Edu: Lower secondary", 
                    "Father's Edu: Upper secondary", "Father's Edu: Bachelor", 
                    "Father's Edu: Master", "Father's Edu: Missing", "Standardized GPA", "GPA Polynomial: 2nd","GPA Polynomial: 3rd", "Age",
                    "2 Norwegian Parents: Norway","2 Norwegian Parents: Western", "2 Norwegian Parents: Non-Western",
                    "Immigrant: Western","Immigrant: Non-Western", "2 Immigrant parents: Western", "2 Immigrant parents: Non-Western",
                    "1 Norwegian Parent: Western", "1 Norwegian Parent: Non-Western","Centrality Index","Western Norway","Eastern Norway","Northern Norway","Southern Norway",
                    "Trøndelag region", "1 child in household","2 children in household","3 children in household", "4 or more children in household", "Children in household missing",
                    "Class size: 0-24", "Class size: 25-49", "Class size: 50-99", "Class size: 100+", "Class size missing", "Lowest income quintile", "Standardized GPA of peers"))


#calculate propensity scores#
psmodel <- glm(TREAT_SKR~gender+examyear_2006+examyear_2007+ examyear_2008+examyear_2009+
                     examyear_2010+ examyear_2011+ examyear_2012+ examyear_2013+ 
                     examyear_2014+ examyear_2015+ examyear_2016 + examyear_2017 +
                     examyear_2018+ STP_1+ STP_2 + STP_3 + STP_4 + STP_5 + STP_6 +
                     mor_edu_1+ mor_edu_2+ mor_edu_3 + mor_edu_4 + 
                     mor_edu_9+ far_edu_1 + far_edu_2 + far_edu_3 + far_edu_4 +
                     far_edu_9+ poly(gspstd,3) + age + immback_AG0 + immback_AG1 + immback_AG2 +
                     immback_B1 + immback_B2 + immback_C1 + immback_C2 + immback_EF1 + immback_EF2+index+region_vest+region_ost+ 
                     region_nord+region_sor+region_trondelag+children_cat_1+children_cat_2+children_cat_3+
                     children_cat_4+children_cat_9+class_size_cat_0_24+class_size_cat_25_49+class_size_cat_50_99+class_size_cat_100+class_size_cat_999+
                     low_inc_indicator+jackknife_gpa_peers, 
                   data = unmatched_df, family = binomial())
summary(psmodel)
unmatched_df$p <- predict(psmodel, newdata = unmatched_df, type = "response")
unmatched_df$ps <- log((1 - unmatched_df$p) / (unmatched_df$p))


##### Nearest Neighbor matching with 0.25*SD caliper, without replacement
# and covariate balance

set.seed(100)
match_out <-
  matchit(TREAT_SKR~gender+examyear_2006+examyear_2007+ examyear_2008+examyear_2009+
            examyear_2010+ examyear_2011+ examyear_2012+ examyear_2013+ 
            examyear_2014+ examyear_2015+ examyear_2016 + examyear_2017 +
            examyear_2018+ STP_1+ STP_2 + STP_3 + STP_4 + STP_5 + STP_6 +
            mor_edu_1+ mor_edu_2+ mor_edu_3 + mor_edu_4 + 
            mor_edu_9+ far_edu_1 + far_edu_2 + far_edu_3 + far_edu_4 +
            far_edu_9+ poly(gspstd,3)+ age + immback_AG0 + immback_AG1 + immback_AG2 +
            immback_B1 + immback_B2 + immback_C1 + immback_C2 + immback_EF1 + immback_EF2+index+region_vest+region_ost+ 
            region_nord+region_sor+region_trondelag+children_cat_1+children_cat_2+children_cat_3+
            children_cat_4+children_cat_9+class_size_cat_0_24+class_size_cat_25_49+class_size_cat_50_99+class_size_cat_100+class_size_cat_999+
            low_inc_indicator+jackknife_gpa_peers,
          data = unmatched_df,
          method = "nearest",
          distance = unmatched_df$ps,
          m.order = "random",
          caliper = 0.25,
          replace = FALSE
  )

match_out
match_df<-MatchIt::match.data(match_out)

# covariate balance
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


######## Sensitivity for PSM model variables ########

#### without gender ####

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


## function for covariate balance check
v<-data.frame(old=c("examyear_2006", "examyear_2007", "examyear_2008", "examyear_2009", "examyear_2010", "examyear_2011",
                    "examyear_2012", "examyear_2013", "examyear_2014", "examyear_2015", "examyear_2016", "examyear_2017", "examyear_2018",
                    "STP_1", "STP_2", "STP_3", "STP_4", "STP_5", "STP_6", "mor_edu_1", "mor_edu_2", "mor_edu_3", "mor_edu_4", 
                    "mor_edu_9", "far_edu_1", "far_edu_2", "far_edu_3", "far_edu_4", "far_edu_9", "poly(gspstd, 3)_1", "poly(gspstd, 3)_2",
                    "poly(gspstd, 3)_3", "age","immback_AG0","immback_AG1","immback_AG2","immback_B1","immback_B2","immback_C1","immback_C2",
                    "immback_EF1","immback_EF2"),
              new=c("2006", "2007", "2008", "2009", "2010", "2011","2012", "2013", "2014", "2015", "2016", "2017", "2018",
                    "Teacher Grade: 1", "Teacher Grade: 2", "Teacher Grade: 3", "Teacher Grade: 4", "Teacher Grade: 5", "Teacher Grade: 6", 
                    "Mother's Edu: Lower secondary", "Mother's Edu: Upper secondary", 
                    "Mother's Edu: Bachelor", "Mother's Edu: Master", "Mother's Edu: Missing", "Father's Edu: Lower secondary", 
                    "Father's Edu: Upper secondary", "Father's Edu: Bachelor", 
                    "Father's Edu: Master", "Father's Edu: Missing","Standardized GPA", "GPA Polynomial: 2nd","GPA Polynomial: 3rd", "Age",
                    "2 Norwegian Parents: Norway","2 Norwegian Parents: Western", "2 Norwegian Parents: Non-Western",
                    "Immigrant: Western","Immigrant: Non-Western", "2 Immigrant parents: Western", "2 Immigrant parents: Non-Western",
                    "1 Norwegian Parent: Western", "1 Norwegian Parent: Non-Western"))

#calculate propensity scores#
psmodel <- glm(TREAT_SKR~examyear_2006+examyear_2007+ examyear_2008+examyear_2009+
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
  matchit(TREAT_SKR~examyear_2006+examyear_2007+ examyear_2008+examyear_2009+
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

match_out
match_df<-MatchIt::match.data(match_out)


# covariate balance
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


###### without age ####

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


## function for covariate balance check
v<-data.frame(old=c("gender", "examyear_2006", "examyear_2007", "examyear_2008", "examyear_2009", "examyear_2010", "examyear_2011",
                    "examyear_2012", "examyear_2013", "examyear_2014", "examyear_2015", "examyear_2016", "examyear_2017", "examyear_2018",
                    "STP_1", "STP_2", "STP_3", "STP_4", "STP_5", "STP_6", "mor_edu_1", "mor_edu_2", "mor_edu_3", "mor_edu_4", 
                    "mor_edu_9", "far_edu_1", "far_edu_2", "far_edu_3", "far_edu_4", "far_edu_9", "poly(gspstd, 3)_1", "poly(gspstd, 3)_2",
                    "poly(gspstd, 3)_3", "immback_AG0","immback_AG1","immback_AG2","immback_B1","immback_B2","immback_C1","immback_C2",
                    "immback_EF1","immback_EF2"),
              new=c("Gender", "2006", "2007", "2008", "2009", "2010", "2011","2012", "2013", "2014", "2015", "2016", "2017", "2018",
                    "Teacher Grade: 1", "Teacher Grade: 2", "Teacher Grade: 3", "Teacher Grade: 4", "Teacher Grade: 5", "Teacher Grade: 6", 
                    "Mother's Edu: Lower secondary", "Mother's Edu: Upper secondary", 
                    "Mother's Edu: Bachelor", "Mother's Edu: Master", "Mother's Edu: Missing", "Father's Edu: Lower secondary", 
                    "Father's Edu: Upper secondary", "Father's Edu: Bachelor", 
                    "Father's Edu: Master", "Father's Edu: Missing","Standardized GPA", "GPA Polynomial: 2nd","GPA Polynomial: 3rd", 
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
                 far_edu_9+ poly(gspstd,3)  + immback_AG0 + immback_AG1 + immback_AG2 +
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
            far_edu_9+ poly(gspstd,3) + immback_AG0 + immback_AG1 + immback_AG2 +
            immback_B1 + immback_B2 + immback_C1 + immback_C2 + immback_EF1 + immback_EF2,
          data = unmatched_df,
          method = "nearest",
          distance = unmatched_df$ps,
          m.order = "random",
          caliper = 0.25,
          replace = FALSE
  )

match_out
match_df<-MatchIt::match.data(match_out)


# covariate balance
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


#### without GPA #######

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


## function for covariate balance check
v<-data.frame(old=c("gender", "examyear_2006", "examyear_2007", "examyear_2008", "examyear_2009", "examyear_2010", "examyear_2011",
                    "examyear_2012", "examyear_2013", "examyear_2014", "examyear_2015", "examyear_2016", "examyear_2017", "examyear_2018",
                    "STP_1", "STP_2", "STP_3", "STP_4", "STP_5", "STP_6", "mor_edu_1", "mor_edu_2", "mor_edu_3", "mor_edu_4", 
                    "mor_edu_9", "far_edu_1", "far_edu_2", "far_edu_3", "far_edu_4", "far_edu_9", "age","immback_AG0","immback_AG1",
                    "immback_AG2","immback_B1","immback_B2","immback_C1","immback_C2",
                    "immback_EF1","immback_EF2"),
              new=c("Gender", "2006", "2007", "2008", "2009", "2010", "2011","2012", "2013", "2014", "2015", "2016", "2017", "2018",
                    "Teacher Grade: 1", "Teacher Grade: 2", "Teacher Grade: 3", "Teacher Grade: 4", "Teacher Grade: 5", "Teacher Grade: 6", 
                    "Mother's Edu: Lower secondary", "Mother's Edu: Upper secondary", 
                    "Mother's Edu: Bachelor", "Mother's Edu: Master", "Mother's Edu: Missing", "Father's Edu: Lower secondary", 
                    "Father's Edu: Upper secondary", "Father's Edu: Bachelor", 
                    "Father's Edu: Master", "Father's Edu: Missing", "Age",
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
                 far_edu_9 + age + immback_AG0 + immback_AG1 + immback_AG2 +
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
            far_edu_9+  age + immback_AG0 + immback_AG1 + immback_AG2 +
            immback_B1 + immback_B2 + immback_C1 + immback_C2 + immback_EF1 + immback_EF2,
          data = unmatched_df,
          method = "nearest",
          distance = unmatched_df$ps,
          m.order = "random",
          caliper = 0.25,
          replace = FALSE
  )

match_out
match_df<-MatchIt::match.data(match_out)


# covariate balance
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


####### without teacher-assessed grade #######


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


## function for covariate balance check
v<-data.frame(old=c("gender", "examyear_2006", "examyear_2007", "examyear_2008", "examyear_2009", "examyear_2010", "examyear_2011",
                    "examyear_2012", "examyear_2013", "examyear_2014", "examyear_2015", "examyear_2016", "examyear_2017", "examyear_2018",
                    "mor_edu_1", "mor_edu_2", "mor_edu_3", "mor_edu_4", 
                    "mor_edu_9", "far_edu_1", "far_edu_2", "far_edu_3", "far_edu_4", "far_edu_9", "poly(gspstd, 3)_1", "poly(gspstd, 3)_2",
                    "poly(gspstd, 3)_3", "age","immback_AG0","immback_AG1","immback_AG2","immback_B1","immback_B2","immback_C1","immback_C2",
                    "immback_EF1","immback_EF2"),
              new=c("Gender", "2006", "2007", "2008", "2009", "2010", "2011","2012", "2013", "2014", "2015", "2016", "2017", "2018",
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
                 examyear_2018+
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
            examyear_2018+
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

match_out
match_df<-MatchIt::match.data(match_out)


# covariate balance
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

#### without exam year #########


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


## function for covariate balance check
v<-data.frame(old=c("gender",
                    "STP_1", "STP_2", "STP_3", "STP_4", "STP_5", "STP_6", "mor_edu_1", "mor_edu_2", "mor_edu_3", "mor_edu_4", 
                    "mor_edu_9", "far_edu_1", "far_edu_2", "far_edu_3", "far_edu_4", "far_edu_9", "poly(gspstd, 3)_1", "poly(gspstd, 3)_2",
                    "poly(gspstd, 3)_3", "age","immback_AG0","immback_AG1","immback_AG2","immback_B1","immback_B2","immback_C1","immback_C2",
                    "immback_EF1","immback_EF2"),
              new=c("Gender",
                    "Teacher Grade: 1", "Teacher Grade: 2", "Teacher Grade: 3", "Teacher Grade: 4", "Teacher Grade: 5", "Teacher Grade: 6", 
                    "Mother's Edu: Lower secondary", "Mother's Edu: Upper secondary", 
                    "Mother's Edu: Bachelor", "Mother's Edu: Master", "Mother's Edu: Missing", "Father's Edu: Lower secondary", 
                    "Father's Edu: Upper secondary", "Father's Edu: Bachelor", 
                    "Father's Edu: Master", "Father's Edu: Missing","Standardized GPA", "GPA Polynomial: 2nd","GPA Polynomial: 3rd", "Age",
                    "2 Norwegian Parents: Norway","2 Norwegian Parents: Western", "2 Norwegian Parents: Non-Western",
                    "Immigrant: Western","Immigrant: Non-Western", "2 Immigrant parents: Western", "2 Immigrant parents: Non-Western",
                    "1 Norwegian Parent: Western", "1 Norwegian Parent: Non-Western"))

#calculate propensity scores#
psmodel <- glm(TREAT_SKR~gender + STP_1+ STP_2 + STP_3 + STP_4 + STP_5 + STP_6 +
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
  matchit(TREAT_SKR~gender + STP_1+ STP_2 + STP_3 + STP_4 + STP_5 + STP_6 +
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

match_out
match_df<-MatchIt::match.data(match_out)


# covariate balance
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


######## without immigration background #########


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


## function for covariate balance check
v<-data.frame(old=c("gender", "examyear_2006", "examyear_2007", "examyear_2008", "examyear_2009", "examyear_2010", "examyear_2011",
                    "examyear_2012", "examyear_2013", "examyear_2014", "examyear_2015", "examyear_2016", "examyear_2017", "examyear_2018",
                    "STP_1", "STP_2", "STP_3", "STP_4", "STP_5", "STP_6", "mor_edu_1", "mor_edu_2", "mor_edu_3", "mor_edu_4", 
                    "mor_edu_9", "far_edu_1", "far_edu_2", "far_edu_3", "far_edu_4", "far_edu_9", "poly(gspstd, 3)_1", "poly(gspstd, 3)_2",
                    "poly(gspstd, 3)_3", "age"),
              new=c("Gender", "2006", "2007", "2008", "2009", "2010", "2011","2012", "2013", "2014", "2015", "2016", "2017", "2018",
                    "Teacher Grade: 1", "Teacher Grade: 2", "Teacher Grade: 3", "Teacher Grade: 4", "Teacher Grade: 5", "Teacher Grade: 6", 
                    "Mother's Edu: Lower secondary", "Mother's Edu: Upper secondary", 
                    "Mother's Edu: Bachelor", "Mother's Edu: Master", "Mother's Edu: Missing", "Father's Edu: Lower secondary", 
                    "Father's Edu: Upper secondary", "Father's Edu: Bachelor", 
                    "Father's Edu: Master", "Father's Edu: Missing","Standardized GPA", "GPA Polynomial: 2nd","GPA Polynomial: 3rd", "Age"))

#calculate propensity scores#
psmodel <- glm(TREAT_SKR~gender+examyear_2006+examyear_2007+ examyear_2008+examyear_2009+
                 examyear_2010+ examyear_2011+ examyear_2012+ examyear_2013+ 
                 examyear_2014+ examyear_2015+ examyear_2016 + examyear_2017 +
                 examyear_2018+ STP_1+ STP_2 + STP_3 + STP_4 + STP_5 + STP_6 +
                 mor_edu_1+ mor_edu_2+ mor_edu_3 + mor_edu_4 + 
                 mor_edu_9+ far_edu_1 + far_edu_2 + far_edu_3 + far_edu_4 +
                 far_edu_9+ poly(gspstd,3) + age, 
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
            far_edu_9+ poly(gspstd,3)+ age,
          data = unmatched_df,
          method = "nearest",
          distance = unmatched_df$ps,
          m.order = "random",
          caliper = 0.25,
          replace = FALSE
  )

match_out
match_df<-MatchIt::match.data(match_out)

# covariate balance
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


##### without mother's education ######


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


## function for covariate balance check
v<-data.frame(old=c("gender", "examyear_2006", "examyear_2007", "examyear_2008", "examyear_2009", "examyear_2010", "examyear_2011",
                    "examyear_2012", "examyear_2013", "examyear_2014", "examyear_2015", "examyear_2016", "examyear_2017", "examyear_2018",
                    "STP_1", "STP_2", "STP_3", "STP_4", "STP_5", "STP_6","far_edu_1", "far_edu_2", "far_edu_3", "far_edu_4", "far_edu_9", "poly(gspstd, 3)_1", "poly(gspstd, 3)_2",
                    "poly(gspstd, 3)_3", "age","immback_AG0","immback_AG1","immback_AG2","immback_B1","immback_B2","immback_C1","immback_C2",
                    "immback_EF1","immback_EF2"),
              new=c("Gender", "2006", "2007", "2008", "2009", "2010", "2011","2012", "2013", "2014", "2015", "2016", "2017", "2018",
                    "Teacher Grade: 1", "Teacher Grade: 2", "Teacher Grade: 3", "Teacher Grade: 4", "Teacher Grade: 5", "Teacher Grade: 6", 
                     "Father's Edu: Lower secondary", 
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
                  far_edu_1 + far_edu_2 + far_edu_3 + far_edu_4 +
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
            far_edu_1 + far_edu_2 + far_edu_3 + far_edu_4 +
            far_edu_9+ poly(gspstd,3)+ age + immback_AG0 + immback_AG1 + immback_AG2 +
            immback_B1 + immback_B2 + immback_C1 + immback_C2 + immback_EF1 + immback_EF2,
          data = unmatched_df,
          method = "nearest",
          distance = unmatched_df$ps,
          m.order = "random",
          caliper = 0.25,
          replace = FALSE
  )

match_out
match_df<-MatchIt::match.data(match_out)

# covariate balance
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


#### without father's education #######

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


## function for covariate balance check
v<-data.frame(old=c("gender", "examyear_2006", "examyear_2007", "examyear_2008", "examyear_2009", "examyear_2010", "examyear_2011",
                    "examyear_2012", "examyear_2013", "examyear_2014", "examyear_2015", "examyear_2016", "examyear_2017", "examyear_2018",
                    "STP_1", "STP_2", "STP_3", "STP_4", "STP_5", "STP_6", "mor_edu_1", "mor_edu_2", "mor_edu_3", "mor_edu_4", 
                    "mor_edu_9",  "poly(gspstd, 3)_1", "poly(gspstd, 3)_2",
                    "poly(gspstd, 3)_3", "age","immback_AG0","immback_AG1","immback_AG2","immback_B1","immback_B2","immback_C1","immback_C2",
                    "immback_EF1","immback_EF2"),
              new=c("Gender", "2006", "2007", "2008", "2009", "2010", "2011","2012", "2013", "2014", "2015", "2016", "2017", "2018",
                    "Teacher Grade: 1", "Teacher Grade: 2", "Teacher Grade: 3", "Teacher Grade: 4", "Teacher Grade: 5", "Teacher Grade: 6", 
                    "Mother's Edu: Lower secondary", "Mother's Edu: Upper secondary", 
                    "Mother's Edu: Bachelor", "Mother's Edu: Master", "Mother's Edu: Missing", "Standardized GPA", "GPA Polynomial: 2nd","GPA Polynomial: 3rd", "Age",
                    "2 Norwegian Parents: Norway","2 Norwegian Parents: Western", "2 Norwegian Parents: Non-Western",
                    "Immigrant: Western","Immigrant: Non-Western", "2 Immigrant parents: Western", "2 Immigrant parents: Non-Western",
                    "1 Norwegian Parent: Western", "1 Norwegian Parent: Non-Western"))

#calculate propensity scores#
psmodel <- glm(TREAT_SKR~gender+examyear_2006+examyear_2007+ examyear_2008+examyear_2009+
                 examyear_2010+ examyear_2011+ examyear_2012+ examyear_2013+ 
                 examyear_2014+ examyear_2015+ examyear_2016 + examyear_2017 +
                 examyear_2018+ STP_1+ STP_2 + STP_3 + STP_4 + STP_5 + STP_6 +
                 mor_edu_1+ mor_edu_2+ mor_edu_3 + mor_edu_4 + 
                 mor_edu_9+ poly(gspstd,3) + age + immback_AG0 + immback_AG1 + immback_AG2 +
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
            mor_edu_9+ poly(gspstd,3)+ age + immback_AG0 + immback_AG1 + immback_AG2 +
            immback_B1 + immback_B2 + immback_C1 + immback_C2 + immback_EF1 + immback_EF2,
          data = unmatched_df,
          method = "nearest",
          distance = unmatched_df$ps,
          m.order = "random",
          caliper = 0.25,
          replace = FALSE
  )

match_out
match_df<-MatchIt::match.data(match_out)


# covariate balance
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


########### Placebo treatment matching ##############

### matching 2 (treat) vs 3 (control) ####

# load data with entire sample of students who received grades 1-6

## create a new dummy variable for treat (grade 2) and control (grade 3)

population <- population %>% 
  mutate(TREAT_SKR =  case_when((SKR==2) ~ 1,
                                (SKR == 3) ~ 0))

#keep those with a 2 or 3
population<- population %>% 
  filter(!is.na(TREAT_SKR)) 

# replace NA values for mother and fathers education with missing value to create separate "missing" category in PSM
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


## function for covariate balance check
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

match_out
match_df<-MatchIt::match.data(match_out)


# covariate balance
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



### matching 3 (treat) vs 4 (control) ####

# load data with entire sample of students who received grades 1-6

## create a new dummy variable for treat (grade 3) and control (grade 4)

population <- population %>% 
  mutate(TREAT_SKR =  case_when((SKR==3) ~ 1,
                                (SKR == 4) ~ 0))

#keep those with a 3 or 4
population<- population %>% 
  filter(!is.na(TREAT_SKR)) 

# replace NA values for mother and fathers education with missing value to create separate "missing" category in PSM
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


## function for covariate balance check
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

match_out
match_df<-MatchIt::match.data(match_out)


# covariate balance
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



###### General code for placebo outcomes analyses (non-psychological diagnoses) #####

# load general physician visits data (specific to each placebo outcome) and matched data

## create an end date which is 1.5 yrs after exam date
# then merge basic exam dates with kuhr data

match_df$exammonth <- 06
match_df$examday <- 30
match_df$examyear <- as.numeric(as.character(match_df$examyear))

match_df <- match_df %>%
  mutate(examdate = make_datetime(examyear, exammonth, examday))

# remove unnecessary variables
match_df <- select(match_df,-c("exammonth":"examday"))

# create date which is 1.5 years after exam date
match_df$enddate<-ymd(match_df$examdate) + days(549)

# temporary df with exam and end date to merge with health data
temp_df<-match_df %>% select(ID, examdate, enddate)

# select variables from health data
kuhr<-kuhr %>% select(ID, DATO, DIAGNOSER, chap)

# merge dates and health data
temp_df<-left_join(temp_df, kuhr, by="ID")


#keep only kuhr data from within the exam date and end date 

temp_df<-temp_df %>% filter((DATO>=examdate)&(DATO<=enddate))

#remerge with matched df
examkuhr<-left_join(match_df, temp_df, by="ID")


##Create event variable if diagnosis
examkuhr <- examkuhr %>%
  mutate(event = case_when(is.na(chap) ~ 0,
                           (!is.na(chap) ~ 1)))

#count total # of events by individual
count <-
  aggregate(examkuhr$event,
            by = list(ID = examkuhr$ID),
            FUN = sum)

# select variables for analyses
temp_df1<-match_df %>% dplyr:: select(ID, TREAT_SKR,gender,STP, gspstd, age, mor_edu, far_edu, examyear, immback, inc_quant)

# merge with count data
count<-left_join(count, temp_df1, by="ID")


# clean up for analysis
count <- count %>% 
  dplyr::rename("num" = x)


## make binary variable if at least one 
count<-count %>% 
  mutate(zero=case_when(num==0~0, num>0~1))

# set variables as factors
count$TREAT_SKR<-as.factor(count$TREAT_SKR)
count$STP<-as.factor(count$STP)
count$gender<-as.factor(count$gender)
count$examyear<-as.factor(count$examyear)
count$mor_edu<-as.factor(count$mor_edu)
count$far_edu<-as.factor(count$far_edu)



# Negative Binomial
nb <-
  glm.nb(
    num ~ TREAT_SKR + gender + STP + examyear + mor_edu + far_edu +
      poly(gspstd, 3) + age + immback,
    data = count,
    link = log
  )
summary(nb)
summ(
  nb,
  confint = TRUE,
  ci.width = 0.95,
  robust = TRUE,
  exp = T
)

# logistic regression
log <-
  glm(
    zero ~ TREAT_SKR + gender + STP + examyear + mor_edu + far_edu +
      poly(gspstd, 3) + age + immback,
    data = count,
    family = binomial(link = "logit")
  )
summary(log)

summ(
  log,
  confint = TRUE,
  ci.width = 0.95,
  robust = TRUE,
  exp = T
)
