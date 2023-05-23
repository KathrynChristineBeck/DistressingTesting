##### pre and post trends analysis ############

# load libraries
library(tidyverse)
library(arsenal)
library(data.table)
library(ggsci)
library(lubridate)
library(jtools)
library(MASS)
select <- dplyr::select

# load matched sample and restrict to school years 2008/2009-2016/2017 and 
# those registered as living in Norway the 3 years prior to exam

# then load psychological diagnosis data



##### pre and post trends analysis ORs ############


# follow same procedure as main analysis but this time separate to diagnoses each year for 3 
# years prior to two years after
temp_df <- match_df


# create exam date
temp_df$exammonth <- 06
temp_df$examday <- 30
temp_df$examyear <- as.numeric(as.character(temp_df$examyear))

temp_df <- temp_df %>%
  mutate(examdate = make_datetime(examyear, exammonth, examday))

# remove unnecessary variables
temp_df <- select(temp_df,-c("exammonth":"examday"))

temp_df$pre3yr<-ymd(temp_df$examdate) - days(1098)
temp_df$pre2yr<-ymd(temp_df$examdate) - days(732)
temp_df$pre1yr<-ymd(temp_df$examdate) - days(366)
temp_df$time0<-ymd(temp_df$examdate)
temp_df$post1yr<-ymd(temp_df$examdate) + days(366)
temp_df$post2yr<-ymd(temp_df$examdate) + days(732)


# temporary df with exam dates to merge with health data
temp_df1<-temp_df %>% select(ID, "examdate":"post2yr")

# select variables from health data
kuhr<-kuhr %>% select(ID, DATO, DIAGNOSER, chap)

# merge dates and health data
temp_df1<-left_join(temp_df1, kuhr, by="ID")

#keep only kuhr data from within the pre date and end date and remove
#diagnoses which are unrelated

temp_df1<-temp_df1 %>% filter((DATO>=pre3yr)&(DATO<=post2yr)) %>% 
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
examkuhr<-left_join(match_df, temp_df1, by="ID")

# remove unnecessary dfs
rm(temp_df, temp_df1)

######## create event for each pre and post period
examkuhr <- examkuhr %>%
  group_by(ID) %>% 
  mutate(eventpre3 = case_when((!is.na(chap))&(DATO>=pre3yr)&(DATO<pre2yr) ~ 1, TRUE~0)) %>% 
  mutate(eventpre2 = case_when((!is.na(chap))&(DATO>=pre2yr)&(DATO<pre1yr) ~ 1, TRUE~0)) %>% 
  mutate(eventpre1 = case_when((!is.na(chap))&(DATO>=pre1yr)&(DATO<time0) ~ 1, TRUE~0)) %>% 
  mutate(eventpost1 = case_when((!is.na(chap))&(DATO>=time0)&(DATO<post1yr) ~ 1, TRUE~0)) %>% 
  mutate(eventpost2 = case_when((!is.na(chap))&(DATO>=post1yr)&(DATO<=post2yr) ~ 1, TRUE~0)) %>% 
  ungroup()


## keep only one obs per ind
examkuhr <- examkuhr %>% 
  distinct(ID, .keep_all = T)

# make sure variables are factors
examkuhr$TREAT_SKR<-as.factor(examkuhr$TREAT_SKR)
examkuhr$STP<-as.factor(examkuhr$STP)
examkuhr$gender<-as.factor(examkuhr$gender)
examkuhr$examyear<-as.factor(examkuhr$examyear)
examkuhr$mor_edu<-as.factor(examkuhr$mor_edu)
examkuhr$far_edu<-as.factor(examkuhr$far_edu)

# run one regression for each period


#pre 3 years
regpre3 <-
  glm(
    eventpre3 ~ TREAT_SKR + gender + STP + examyear + mor_edu + far_edu + poly(gspstd, 3) + age + immback ,
    data = examkuhr,
    family = "binomial"
  )
summary(regpre3)
exppre3 <-
  summ(
    regpre3,
    confint = TRUE,
    ci.width = 0.95,
    robust = TRUE,
    exp = T
  )


#pre 2 years
regpre2 <-
  glm(
    eventpre2 ~ TREAT_SKR + gender + STP + examyear + mor_edu + far_edu + poly(gspstd, 3) + age + immback,
    data = examkuhr,
    family = "binomial"
  )
summary(regpre2)
exppre2 <-
  summ(
    regpre2,
    confint = TRUE,
    ci.width = 0.95,
    robust = TRUE,
    exp = T
  )


#exam year
regpre1 <-
  glm(
    eventpre1 ~ TREAT_SKR + gender + STP + examyear + mor_edu + far_edu + poly(gspstd, 3) + age + immback ,
    data = examkuhr,
    family = "binomial"
  )
summary(regpre1)
exppre1 <-
  summ(
    regpre1,
    confint = TRUE,
    ci.width = 0.95,
    robust = TRUE,
    exp = T
  )


# post 1 year
regpost1 <-
  glm(
    eventpost1 ~ TREAT_SKR + gender + STP + examyear + mor_edu + far_edu + poly(gspstd, 3) + age + immback ,
    data = examkuhr,
    family = "binomial"
  )
summary(regpost1)
exppost1 <-
  summ(
    regpost1,
    confint = TRUE,
    ci.width = 0.95,
    robust = TRUE,
    exp = T
  )


# post 2 years
regpost2 <-
  glm(
    eventpost2 ~ TREAT_SKR + gender + STP + examyear + mor_edu + far_edu + poly(gspstd, 3) + age + immback,
    data = examkuhr,
    family = "binomial"
  )
summary(regpost2)
exppost2 <-
  summ(
    regpost2,
    confint = TRUE,
    ci.width = 0.95,
    robust = TRUE,
    exp = T
  )


## Create regression plot with OR from each regression
a <- data.frame(exppre3[["coeftable"]])[2,]
b <- data.frame(exppre2[["coeftable"]])[2,]
c <- data.frame(exppre1[["coeftable"]])[2,]
d <- data.frame(exppost1[["coeftable"]])[2,]
e <- data.frame(exppost2[["coeftable"]])[2,]


rownames(a)<-c("-3")
rownames(b)<-c("-2")
rownames(c)<-c("-1")
rownames(d)<-c("0")
rownames(e)<-c("1")


plot<-bind_rows(a,b,c,d,e)

plot$time <-
  c(
    "-3",
    "-2",
    "-1",
    "0",
    "1"
  )

plot <- plot %>% 
  dplyr::rename("lower" = 'X2.5.') %>% 
  dplyr::rename("upper" = 'X97.5.') %>% 
  dplyr::rename("OR" = 'exp.Est..')

plot$time<-as.character(plot$time)
plot$time<-factor(plot$time, levels = unique(plot$time))

plot$time<-as.numeric(as.character(plot$time))

# plot
plot %>%
  ggplot(aes(
    x = time,
    y = OR,
    color = as.factor(time),
    group = 1
  )) +
  geom_line(color = c("#374E55FF")) +
  geom_point() +
  geom_linerange(aes(ymin = lower, ymax = upper)) +
  geom_hline(
    yintercept = 1,
    linetype = "dotted",
    alpha = 0.8,
    color = "#374E55FF"
  ) +
  theme_bw() +
  labs(x = "Years Since Exam", y = "Odds Ratio") +
  theme(
    legend.position = "none",
    panel.grid.minor = element_blank(),
    axis.title.x = element_text(margin = margin(t = 10)),
    axis.title.y = element_text(margin = margin(r = 10)),
    plot.caption = element_text(size = 8, margin = margin(t = 15)),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 12),
    plot.subtitle = element_text(size = 9, hjust = 0.5)
  ) +
  scale_color_manual(values = c(
    "#374E55FF",
    "#374E55FF",
    "#374E55FF",
    "#DF8F44FF",
    "#374E55FF"
  )) +
  scale_y_continuous(breaks = c(0.8, 1, 1.2, 1.4, 1.6),
                     limits = c(0.8, 1.65)) +
  scale_x_continuous(breaks = c(-3, -2, -1, 0, 1),
                     labels = c(-2, -1, "Exam Year", 1, 2))

########### pre and post trends IRR ################

# follow same procedure as main analysis but use number of visits 
temp_df <- match_df

# create exam date
temp_df$exammonth <- 06
temp_df$examday <- 30
temp_df$examyear <- as.numeric(as.character(temp_df$examyear))

temp_df <- temp_df %>%
  mutate(examdate = make_datetime(examyear, exammonth, examday))

# remove unnecessary variables
temp_df <- select(temp_df,-c("exammonth":"examday"))

temp_df$pre3yr<-ymd(temp_df$examdate) - days(1098)
temp_df$pre2yr<-ymd(temp_df$examdate) - days(732)
temp_df$pre1yr<-ymd(temp_df$examdate) - days(366)
temp_df$time0<-ymd(temp_df$examdate)
temp_df$post1yr<-ymd(temp_df$examdate) + days(366)
temp_df$post2yr<-ymd(temp_df$examdate) + days(732)


# temporary df with exam dates to merge with health data
temp_df1<-temp_df %>% select(ID, "examdate":"post2yr")

# select variables from health data
kuhr<-kuhr %>% select(ID, DATO, DIAGNOSER, chap)

# merge dates and health data
temp_df1<-left_join(temp_df1, kuhr, by="ID")

#keep only kuhr data from within the pre date and end date and remove
#diagnoses which are unrelated

temp_df1<-temp_df1 %>% filter((DATO>=pre3yr)&(DATO<=post2yr)) %>% 
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
examkuhr<-left_join(match_df, temp_df1, by="ID")

# remove unnecessary dfs
rm(temp_df, temp_df1)

######## create event for each pre and post period
examkuhr <- examkuhr %>%
  group_by(ID) %>% 
  mutate(eventpre3 = case_when((!is.na(chap))&(DATO>=pre3yr)&(DATO<pre2yr) ~ 1, TRUE~0)) %>% 
  mutate(eventpre2 = case_when((!is.na(chap))&(DATO>=pre2yr)&(DATO<pre1yr) ~ 1, TRUE~0)) %>% 
  mutate(eventpre1 = case_when((!is.na(chap))&(DATO>=pre1yr)&(DATO<time0) ~ 1, TRUE~0)) %>% 
  mutate(eventpost1 = case_when((!is.na(chap))&(DATO>=time0)&(DATO<post1yr) ~ 1, TRUE~0)) %>% 
  mutate(eventpost2 = case_when((!is.na(chap))&(DATO>=post1yr)&(DATO<=post2yr) ~ 1, TRUE~0)) %>% 
  ungroup()

### count number of visits in each period
countpre3 <-
  aggregate(examkuhr$eventpre3,
            by = list(ID = examkuhr$ID),
            FUN = sum)
countpre3$pre3num <- countpre3$x
countpre3$x <- NULL

countpre2 <-
  aggregate(examkuhr$eventpre2,
            by = list(ID = examkuhr$ID),
            FUN = sum)
countpre2$pre2num <- countpre2$x
countpre2$x <- NULL

countpre1 <-
  aggregate(examkuhr$eventpre1,
            by = list(ID = examkuhr$ID),
            FUN = sum)
countpre1$pre1num <- countpre1$x
countpre1$x <- NULL

countpost1 <-
  aggregate(examkuhr$eventpost1,
            by = list(ID = examkuhr$ID),
            FUN = sum)
countpost1$post1num <- countpost1$x
countpost1$x <- NULL

countpost2 <-
  aggregate(examkuhr$eventpost2,
            by = list(ID = examkuhr$ID),
            FUN = sum)
countpost2$post2num <- countpost2$x
countpost2$x <- NULL


# join to one dataframe
examkuhr<-left_join(match_df, countpre3, by="ID") %>% 
  left_join(countpre2, by="ID") %>% 
  left_join(countpre1, by="ID") %>% 
  left_join(countpost1, by="ID") %>% 
  left_join(countpost2, by="ID")

# make sure variables are factors
examkuhr$TREAT_SKR<-as.factor(examkuhr$TREAT_SKR)
examkuhr$STP<-as.factor(examkuhr$STP)
examkuhr$gender<-as.factor(examkuhr$gender)
examkuhr$examyear<-as.factor(examkuhr$examyear)
examkuhr$mor_edu<-as.factor(examkuhr$mor_edu)
examkuhr$far_edu<-as.factor(examkuhr$far_edu)

# run regression on each pre and post period
#pre 3 years
regpre3 <-
  glm.nb(
    pre3num ~ TREAT_SKR + gender + STP + examyear + mor_edu + far_edu + poly(gspstd, 3) + age + immback,
    data = examkuhr,
    link = log
  )
summary(regpre3)
exppre3 <-
  summ(
    regpre3,
    confint = TRUE,
    ci.width = 0.95,
    robust = TRUE,
    exp = T
  )


#pre 2 years
regpre2 <-
  glm.nb(
    pre2num ~ TREAT_SKR + gender + STP + examyear + mor_edu + far_edu + poly(gspstd, 3) + age + immback,
    data = eventstudy,
    link = log
  )
summary(regpre2)
exppre2 <-
  summ(
    regpre2,
    confint = TRUE,
    ci.width = 0.95,
    robust = TRUE,
    exp = T
  )


#pre 1 years
regpre1 <-
  glm.nb(
    pre1num ~ TREAT_SKR + gender + STP + examyear + mor_edu + far_edu + poly(gspstd, 3) + age + immback,
    data = examkuhr,
    link = log
  )
summary(regpre1)
exppre1 <-
  summ(
    regpre1,
    confint = TRUE,
    ci.width = 0.95,
    robust = TRUE,
    exp = T
  )


# post 1 year
regpost1 <-
  glm.nb(
    post1num ~ TREAT_SKR + gender + STP + examyear + mor_edu + far_edu + poly(gspstd, 3) + age + immback,
    data = examkuhr,
    link = log
  )
summary(regpost1)
exppost1 <-
  summ(
    regpost1,
    confint = TRUE,
    ci.width = 0.95,
    robust = TRUE,
    exp = T
  )


# post 2 years
regpost2 <-
  glm.nb(
    post2num ~ TREAT_SKR + gender + STP + examyear + mor_edu + far_edu + poly(gspstd, 3) + age + immback,
    data = examkuhr,
    link = log
  )
summary(regpost2)
exppost2 <-
  summ(
    regpost2,
    confint = TRUE,
    ci.width = 0.95,
    robust = TRUE,
    exp = T
  )

## Create regression plot with IRR from each regression
a <- data.frame(exppre3[["coeftable"]])[2,]
b <- data.frame(exppre2[["coeftable"]])[2,]
c <- data.frame(exppre1[["coeftable"]])[2,]
d <- data.frame(exppost1[["coeftable"]])[2,]
e <- data.frame(exppost2[["coeftable"]])[2,]

rownames(a)<-c("-3")
rownames(b)<-c("-2")
rownames(c)<-c("-1")
rownames(d)<-c("0")
rownames(e)<-c("1")


plot<-bind_rows(a,b,c,d,e)

plot$time <-
  c(
    "-3",
    "-2",
    "-1",
    "0",
    "1"
  )


plot <- plot %>% 
  dplyr::rename("lower" = 'X2.5.') %>% 
  dplyr::rename("upper" = 'X97.5.') %>% 
  dplyr::rename("IRR" = 'exp.Est..')

plot$time<-as.character(plot$time)
plot$time<-factor(plot$time, levels = unique(plot$time))
plot$time<-as.numeric(as.character(plot$time))

# plot
plot %>%
  ggplot(aes(
    x = time,
    y = IRR,
    color = as.factor(time),
    group = 1
  )) +
  geom_line(color = c("#374E55FF")) +
  geom_point() +
  geom_linerange(aes(ymin = lower, ymax = upper)) +
  geom_hline(
    yintercept = 1,
    linetype = "dotted",
    alpha = 0.8,
    color = "#374E55FF"
  ) +
  theme_bw() +
  labs(x = "Years Since Exam", y = "Incidence Rate Ratio") +
  theme(
    legend.position = "none",
    panel.grid.minor = element_blank(),
    axis.title.x = element_text(margin = margin(t = 10)),
    axis.title.y = element_text(margin = margin(r = 10)),
    plot.caption = element_text(size = 8, margin = margin(t = 15)),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 12),
    plot.subtitle = element_text(size = 9, hjust = 0.5)
  ) +
  scale_color_manual(values = c(
    "#374E55FF",
    "#374E55FF",
    "#374E55FF",
    "#DF8F44FF",
    "#DF8F44FF"
  )) +
  scale_y_continuous(breaks = c(0.8, 1, 1.2, 1.4, 1.6),
                     limits = c(0.8, 1.65)) +
  scale_x_continuous(
    breaks = c(-3, -2, -1, 0, 1),
    labels = c(-2, -1, "Exam Year", 1, 2),
    limits = c(-3, 1.5)
  )


