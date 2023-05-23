 ############## Analyses ##############
 
 # load libraries
 library(tidyverse)
 library(lubridate)
 library(MASS)
 library(pscl)
 library(countreg)
 library(jtools)
 select <- dplyr::select
 
 # load general physician visits data (psychological diagnoses) and matched data
 
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
 
 
 #keep only kuhr data from within the exam date and end date and remove
 #diagnoses which are unrelated
 
 temp_df<-temp_df %>% filter((DATO>=examdate)&(DATO<=enddate)) %>% 
   filter(DIAGNOSER!="P80") %>% 
   filter(DIAGNOSER!="P09") %>% 
   filter(DIAGNOSER!="P12") %>% 
   filter(DIAGNOSER!="P72") %>% 
   filter(DIAGNOSER!="P05") %>%
   filter(DIAGNOSER!="P11") %>% 
   filter(DIAGNOSER!="P13") %>%
   filter(DIAGNOSER!="P70") %>% 
   filter(DIAGNOSER!="P85")
 
 #remerge with matched df
 examkuhr<-left_join(match_df, temp_df, by="ID")
 
 
 ##Create event variable if P diagnosis
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

####### subgroups ########


####### by gender ########

# boys
boys <- count %>% 
  filter(gender==0)

# girls
girls <- count %>% 
  filter(gender==1)

# boys regression

# negative binomial
nb_boy <-
  glm.nb(
    num ~ TREAT_SKR + STP + examyear + mor_edu + far_edu + poly(gspstd, 3) + age + immback,
    data = boys,
    link = log
  )
summary(nb_boy)
summ(
  nb_boy,
  confint = TRUE,
  ci.width = 0.95,
  robust = TRUE,
  exp = T
)

# logistic regression
log_boy <-
  glm(
    zero ~ TREAT_SKR + STP + examyear + mor_edu + far_edu + poly(gspstd, 3) + age + immback,
    data = boys,
    family = binomial(link = "logit")
  )
summary(log_boy)

summ(
  log_boy,
  confint = TRUE,
  ci.width = 0.95,
  robust = TRUE,
  exp = T
)

# girls regression

# negative binomial
nb_girl <-
  glm.nb(
    num ~ TREAT_SKR + STP + examyear + mor_edu + far_edu + poly(gspstd, 3) + age + immback,
    data = girls,
    link = log
  )
summary(nb_girl)
summ(
  nb_girl,
  confint = TRUE,
  ci.width = 0.95,
  robust = TRUE,
  exp = T
)

# logistic regression
log_girl <-
  glm(
    zero ~ TREAT_SKR + STP + examyear + mor_edu + far_edu + poly(gspstd, 3) + age + immback,
    data = girls,
    family = binomial(link = "logit")
  )
summary(log_girl)

summ(
  log_girl,
  confint = TRUE,
  ci.width = 0.95,
  robust = TRUE,
  exp = T,
)



######## teacher grade ########
count$STP <- as.numeric(as.character(count$STP))

# high
high <- count %>% 
  filter(STP>3)

# low
low <- count %>% 
  filter(STP<=3)

# high grade regression

# negative binomial
nb_high <-
  glm.nb(
    num ~ TREAT_SKR + gender +  examyear + mor_edu + far_edu + poly(gspstd, 3) + age + immback,
    data = high,
    link = log
  )
summary(nb_high)
summ(
  nb_high,
  confint = TRUE,
  ci.width = 0.95,
  robust = TRUE,
  exp = T
)


# logistic regression
log_high <-
  glm(
    zero ~ TREAT_SKR + gender + examyear + mor_edu + far_edu + poly(gspstd, 3) + age + immback,
    data = high,
    family = binomial(link = "logit")
  )
summary(log_high)

summ(
  log_high,
  confint = TRUE,
  ci.width = 0.95,
  robust = TRUE,
  exp = T
)



# low grade regression

# negative binomial
nb_low <-
  glm.nb(
    num ~ TREAT_SKR + gender + examyear + mor_edu + far_edu + poly(gspstd, 3) + age + immback,
    data = low,
    link = log
  )
summary(nb_low)
summ(
  nb_low,
  confint = TRUE,
  ci.width = 0.95,
  robust = TRUE,
  exp = T
)

# logistic regression
log_low <-
  glm(
    zero ~ TREAT_SKR + gender + examyear + mor_edu + far_edu + poly(gspstd, 3) + age + immback,
    data = low,
    family = binomial(link = "logit")
  )
summary(log_low)

summ(
  log_low,
  confint = TRUE,
  ci.width = 0.95,
  robust = TRUE,
  exp = T
)



######## immigration background ###########

###### at least 1 minority parent #########
# subgroup

minority_parent <- count %>% 
  filter(immback == "B1" | immback == "B2" | immback == "C1" |immback == "C2" | immback == "EF1" |immback == "EF2")

# at least one minority parent
nb_min <-
  glm.nb(
    num ~ TREAT_SKR + gender + STP + as.numeric(as.character(examyear)) +
      mor_edu + far_edu + poly(gspstd, 3) + age ,
    data = minority_parent,
    link = log
  )
summary(nb_min)
summ(
  nb_min,
  confint = TRUE,
  ci.width = 0.95,
  robust = TRUE,
  exp = T
)


# logistic regression
log_min <-
  glm(
    zero ~ TREAT_SKR + gender + STP + as.numeric(as.character(examyear)) +
      mor_edu + far_edu + poly(gspstd, 3) + age,
    data = minority_parent,
    family = binomial(link = "logit")
  )
summary(log_min)
summ(
  log_min,
  confint = TRUE,
  ci.width = 0.95,
  robust = TRUE,
  exp = T
)


#### two Norwegian born parents ######
norske_parents <- count %>% 
  filter(immback == "AG0" | immback == "AG1" | immback == "AG2")

#  negative binomial
nb_norske <-
  glm.nb(
    num ~ TREAT_SKR + gender + STP + as.numeric(as.character(examyear)) + mor_edu +
      far_edu + poly(gspstd, 3) + age,
    data = norske_parents,
    link = log
  )
summary(nb_norske)
summ(
  nb_norske,
  confint = TRUE,
  ci.width = 0.95,
  robust = TRUE,
  exp = T
)


# logistic regression
log_norske <-
  glm(
    zero ~ TREAT_SKR + gender + STP + gspstd + as.numeric(as.character(examyear)) +
      mor_edu + far_edu + poly(gspstd, 3) + age,
    data = norske_parents,
    family = binomial(link = "logit")
  )
summary(log_norske)
summ(
  log_norske,
  confint = TRUE,
  ci.width = 0.95,
  robust = TRUE,
  exp = T
)



######## household income ####### 

# below median
# subgroup
bottom <- count %>% 
  filter(inc_quant<=2)

# negative binomial
nb_below <-
  glm.nb(
    num ~ TREAT_SKR + gender + STP + as.numeric(as.character(examyear)) +
      mor_edu + far_edu + poly(gspstd, 3) + age + immback,
    data = bottom,
    link = log
  )
summary(nb_below)
summ(
  nb_below,
  confint = TRUE,
  ci.width = 0.95,
  robust = TRUE,
  exp = T
)

# logistic regression
log_below <-
  glm(
    zero ~ TREAT_SKR + gender + STP + as.numeric(as.character(examyear)) +
      mor_edu + far_edu + poly(gspstd, 3) + age + immback,
    data = bottom,
    family = binomial(link = "logit")
  )
summary(log_below)
summ(
  log_below,
  confint = TRUE,
  ci.width = 0.95,
  robust = TRUE,
  exp = T
)


# top 50%
top <- count %>% 
  filter(inc_quant>=3)

# negative binomial
nb_above <-
  glm.nb(
    num ~ TREAT_SKR + gender + STP + as.numeric(as.character(examyear)) + mor_edu +
      far_edu + poly(gspstd, 3) + age + immback,
    data = top,
    link = log
  )
summary(nb_above)
summ(
  nb_above,
  confint = TRUE,
  ci.width = 0.95,
  robust = TRUE,
  exp = T
)


# logistic regression
log_above <-
  glm(
    zero ~ TREAT_SKR + gender + STP + as.numeric(as.character(examyear)) + mor_edu +
      far_edu + poly(gspstd, 3) + age + immback,
    data = top,
    family = binomial(link = "logit")
  )
summary(log_above)
summ(
  log_above,
  confint = TRUE,
  ci.width = 0.95,
  robust = TRUE,
  exp = T
)







#### interactions ####

# gender
nb_gender <-
  glm.nb(
    num ~ TREAT_SKR * gender + as.character(as.numeric(STP)) + examyear + mor_edu +
      far_edu + poly(gspstd, 3) + age + immback,
    data = count,
    link = log
  )
summary(nb_gender)
summ(
  nb_gender,
  confint = F,
  robust = TRUE,
  exp = F,
  digits = getOption("jtools-digits", 3)
)


# STP 
count <- count %>%
  mutate(STP_group = case_when(STP <= 3 ~ 0, STP >= 4 ~ 1))

nb_STP <-
  glm.nb(
    num ~ TREAT_SKR * STP_group + gender + examyear + mor_edu + far_edu + poly(gspstd, 3) + age + immback,
    data = count,
    link = log
  )
summary(nb_STP)
summ(
  nb_STP,
  confint = F,
  robust = TRUE,
  exp = F,
  digits = getOption("jtools-digits", 3)
)


#immigration background
count <- count %>% 
  mutate(imm_parents = case_when((immback == "B1" | immback == "B2" | immback == "C1" |immback == "C2" | immback == "EF1" |immback == "EF2")~1,
                                 (immback == "AG0" | immback == "AG1" | immback == "AG2")~0))           

nb_imm <-
  glm.nb(
    num ~ TREAT_SKR * imm_parents + gender + as.character(as.numeric(STP)) +
      examyear + mor_edu + far_edu + poly(gspstd, 3) + age,
    data = count,
    link = log
  )
summary(nb_imm)
summ(
  nb_imm,
  confint = F,
  robust = TRUE,
  exp = F,
  digits = getOption("jtools-digits", 3)
)





# income
count <- count %>% 
  mutate(inc_med = case_when(inc_quant<=2~0,inc_quant>=3~1))

nb_income <-
  glm.nb(
    num ~ TREAT_SKR * inc_med + gender + as.character(as.numeric(STP)) +
      examyear + mor_edu + far_edu + poly(gspstd, 3) + age + immback,
    data = count,
    link = log
  )
summary(nb_income)
summ(
  nb_income,
  confint = F,
  robust = TRUE,
  exp = F,
  digits = getOption("jtools-digits", 3)
)




######### graduate from upper secondary ######

# load educational attainment and enrollment data and merge to matched data (called "education" -  educational data is in long format by year)

# graduation, on time
df <- education %>%
  group_by(ID) %>%
  mutate(grad = case_when((year == examyear) &
                            (bu_igang == "BU") & ((NUS_code == 4)|(NUS_code == 5)) ~ 1,
                          TRUE ~ 0)) %>%
  mutate(graduation = case_when(any(grad == 1) ~ 1, TRUE ~ 0)) %>%
  ungroup()

# keep one obs. per person
df <- df %>%
  distinct(ID, .keep_all = T)


# graduation within 1 year
df1 <- education %>%
  # subset to exam years 2017 and earlier
  filter(examyear <= 2017) %>%
  group_by(ID) %>%
  mutate(grad = case_when((year >= examyear) &
                            (year <= examyear + 1) & (bu_igang == "BU") & ((NUS_code == 4)|(NUS_code == 5)) ~ 1,
                          TRUE ~ 0
  )) %>%
  mutate(graduation = case_when(any(grad == 1) ~ 1, TRUE ~ 0)) %>%
  ungroup()

# keep one obs. per person
df1 <- df1 %>%
  distinct(ID, .keep_all = T)



# graduation within 5 years
df5 <- education %>%
  # subset to exam years 2013 and earlier
  filter(examyear <= 2013) %>%
  group_by(ID) %>%
  mutate(grad = case_when((year >= examyear) &
                            (year <= examyear + 5) & (bu_igang == "BU") & ((NUS_code == 4)|(NUS_code == 5)) ~ 1,
                          TRUE ~ 0
  )) %>%
  mutate(graduation = case_when(any(grad == 1) ~ 1, TRUE ~ 0)) %>%
  ungroup()


# keep one obs. per person
df5 <- df5 %>%
  distinct(ID, .keep_all = T)


#### Regressions

# On-time
df_reg <-
  glm(
    graduation ~ TREAT_SKR + gender + STP + examyear + mor_edu + far_edu + poly(gspstd, 3) + age + immback,
    data = df,
    family = binomial(link = "logit")
  )
summary(df_reg)
summ(
  df_reg,
  confint = TRUE,
  ci.width = 0.95,
  robust = TRUE,
  exp = T
)


# One year
df1_reg <-
  glm(
    graduation ~ TREAT_SKR + gender + STP + examyear + mor_edu + far_edu + poly(gspstd, 3) + age + immback,
    data = df1,
    family = binomial(link = "logit")
  )
summary(df1_reg)

summ(
  df1_reg,
  confint = TRUE,
  ci.width = 0.95,
  robust = TRUE,
  exp = T
)


# Five years
df5_reg <-
  glm(
    graduation ~ TREAT_SKR + gender + STP + examyear + mor_edu + far_edu + poly(gspstd, 3) + age + immback,
    data = df5,
    family = binomial(link = "logit")
  )
summary(df5_reg)

summ(
  df5_reg,
  confint = TRUE,
  ci.width = 0.95,
  robust = TRUE,
  exp = T
)





############### Enrolled in tertiary education ###################################

# enrolled, same year
df_enroll <- education %>%
  group_by(ID) %>%
  mutate(enroll = case_when((year == examyear) &
                              (bu_igang == "igang") & (NUS_code >= 6) & (NUS_code != 9) ~ 1,
                            TRUE ~ 0
  )) %>%
  mutate(enrolled = case_when(any(enroll == 1) ~ 1, TRUE ~ 0)) %>%
  ungroup()

# keep one obs. per person
df_enroll <- df_enroll %>%
  distinct(ID, .keep_all = T)


# enrolled within 1 year
df1_enroll <- education %>%
  # subset to years <= 2017
  filter(examyear <= 2017) %>%
  group_by(ID) %>%
  mutate(enroll = case_when((year >= examyear) &
                              (year <= examyear + 1) &
                              (bu_igang == "igang") & (NUS_code >= 6) & (NUS_code != 9) ~ 1,
                            TRUE ~ 0
  )) %>%
  mutate(enrolled = case_when(any(enroll == 1) ~ 1, TRUE ~ 0)) %>%
  ungroup()

# keep one obs. per person
df1_enroll <- df1_enroll %>%
  distinct(ID, .keep_all = T)



# enrolled within 5 years
df5_enroll <- education %>%
  # subset to years <= 2013
  filter(examyear <= 2013) %>%
  group_by(ID) %>%
  mutate(enroll = case_when((year >= examyear) &
                              (year <= examyear + 5) &
                              (bu_igang == "igang") & (NUS_code >= 6) & (NUS_code != 9) ~ 1,
                            TRUE ~ 0
  )) %>%
  mutate(enrolled = case_when(any(enroll == 1) ~ 1, TRUE ~ 0)) %>%
  ungroup()

# keep one obs. per person
df5_enroll <- df5_enroll %>%
  distinct(ID, .keep_all = T)



#### Regressions

# On-time
dfenroll_reg <-
  glm(
    enrolled ~ TREAT_SKR + gender + STP + examyear + mor_edu + far_edu + poly(gspstd, 3) + age + immback,
    data = df_enroll,
    family = binomial(link = "logit")
  )
summary(dfenroll_reg)

summ(
  dfenroll_reg,
  confint = TRUE,
  ci.width = 0.95,
  robust = TRUE,
  exp = T
)


# One year
df1enroll_reg <-
  glm(
    enrolled ~ TREAT_SKR + gender + STP + examyear + mor_edu + far_edu + poly(gspstd, 3) + age + immback,
    data = df1_enroll,
    family = binomial(link = "logit")
  )
summary(df1enroll_reg)

summ(
  df1enroll_reg,
  confint = TRUE,
  ci.width = 0.95,
  robust = TRUE,
  exp = T
)


# Five years
df5enroll_reg <-
  glm(
    enrolled ~ TREAT_SKR + gender + STP + examyear + mor_edu + far_edu + poly(gspstd, 3) + age + immback,
    data = df5_enroll,
    family = binomial(link = "logit")
  )
summary(df5enroll_reg)

summ(
  df5enroll_reg,
  confint = TRUE,
  ci.width = 0.95,
  robust = TRUE,
  exp = T
)


