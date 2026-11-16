library(data.table)
library(ggplot2)
theme_set(theme_bw())
library(rms)

# RawData <-

RawData[, NONCARDIOVASCULARDISEASES := PreindexCOPD + PreindexMALIGN + PreindexDemencia]
RawData[, CARDIOVASCULARDISEASES := PreindexCHD + PreindexCVD + PreindexHF + PreindexAF + PreindexHT]
RawData[, METABOLICDISEASES := PreindexDM + PreindexHLP]
RawData[, MINORAMP := as.numeric(!is.na(PreIndexMINORAMPipszidate)) + as.numeric(!is.na(PREIndexMINORAMPkontradate))]
RawData$MINORAMP <- as.factor(ifelse(RawData$MINORAMP > 0, 1, 0))
RawData[, MAJORAMP := PreIndexMAJORAMP]
RawData[, VASCPROCEDURES := PreIndexlastOPENipszi + PreIndexlastOPENkontra + PreIndexlastENDOipszi +
          PreIndexlastENDOkontra]
RawData[, MEDICATION := (PreIndexStatin > 0) + (PreIndexAntithr > 0) + (PreIndexAnticoag > 0) +
          (PostIndexStatin > 0) + (PostIndexAntithr > 0) + (PostIndexAnticoag > 0)]
RawData[, PROCEDURE := PreindexWOUND + PreindexInfection]

RawData$SEX <- as.factor(RawData$SEX)

RawData$TimeToIndexEvent <- as.numeric(difftime(RawData$Indexendodate, RawData$PreIndexPADfirstdate, units = "days"))/365.241

RawData$status <- RawData$PostIndexDeath
RawData$time <- as.numeric(difftime(dplyr::if_else(RawData$status == 1, RawData$DEATHdate, RawData$Censordate),
                                    RawData$Indexendodate, units = "days"))/365.241

RawData$statusMajorAmp <- ifelse(RawData$PostIndexMAJORAMPipszi == 1, "Major amputation",
                                 ifelse(RawData$status == 1, "Death", "Censoring"))
RawData$statusMajorAmp <- factor(RawData$statusMajorAmp, levels = c("Censoring", "Major amputation", "Death"))
RawData$timeMajorAmp <- as.numeric(difftime(
  dplyr::if_else(RawData$statusMajorAmp == "Censoring", RawData$Censordate,
                 dplyr::if_else(RawData$statusMajorAmp == "Major amputation",
                                RawData$PostIndexMAJORAMPipszidate, RawData$DEATHdate)),
  RawData$Indexendodate, units = "days"))/365.241

stopifnot(sum(RawData$PostIndexMAJORAMPipszidate - RawData$PostIndexMINORAMPipszidate < 0, na.rm = TRUE) == 0)
RawData$statusMinorAmp <- ifelse(RawData$PostIndexMinorampipszi == 1, "Minor amputation",
                                 ifelse(RawData$PostIndexMAJORAMPipszi == 1, "Major amputation",
                                        ifelse(RawData$status == 1, "Death", "Censoring")))
RawData$statusMinorAmp <- factor(RawData$statusMinorAmp, levels = c("Censoring", "Minor amputation", "Major amputation", "Death"))
RawData$timeMinorAmp <- as.numeric(difftime(
  dplyr::if_else(RawData$statusMinorAmp == "Censoring", RawData$Censordate,
                 dplyr::if_else(RawData$statusMinorAmp == "Minor amputation",
                                RawData$PostIndexMINORAMPipszidate,
                                dplyr::if_else(RawData$statusMinorAmp == "Major amputation",
                                               RawData$PostIndexMAJORAMPipszidate, RawData$DEATHdate))),
  RawData$Indexendodate, units = "days"))/365.241

RawData2 <- RawData[
  , .(TAJ, time, status, timeMajorAmp, statusMajorAmp, timeMinorAmp, statusMinorAmp,
      AGEindex, SEX, NONCARDIOVASCULARDISEASES, CARDIOVASCULARDISEASES,
      METABOLICDISEASES, PreindexGFR, MINORAMP, MAJORAMP, #RENALDISEASES, MINORAMP, MAJORAMP,
      VASCPROCEDURES, MEDICATION, TimeToIndexEvent, PROCEDURE)]

summary(RawData2)

RawData2 <- na.omit(RawData2)

sum(RawData2$statusMajorAmp == "Major amputation")/sum(RawData2$timeMajorAmp)*1000

lapply(RawData2, function(x) table(x))

################################################################

PostIndexEventDates <- c(
  PostIndeMalign = "PostIndexMaligndate", PostIndexAMI = "PostIndexAMIdate", PostindexSTROKE = "PostindexStrokedate",
  PostIndexMinorampipszi = "PostIndexMINORAMPipszidate", PostIndexMinorampkontra = "PostIndexMINORAMPkontradate",
  PostIndexMAJORAMPipszi = "PostIndexMAJORAMPipszidate", PostIndexMAJORAMPkontra = "PostIndexMAJORAMPkontradate",
  POSTINDEXOPENipszi = "PostIndexOPENipszidate", POSTINDEXOPENkontra = "POSTINDEXOPENkontradate",
  POSTINDEXENDOipszi = "POSTIndexENDOipszidate", POSTINDEXENDOkontra = "POSTINDEXENDOkontradate")

EventTimes <- RawData[, lapply(.SD, function(x) as.numeric(difftime(x, Indexendodate, units = "days"))/365.241),
                      .SDcols = PostIndexEventDates]
setnames(EventTimes, names(PostIndexEventDates))
EventTimes[, `:=`(id = .I, time = RawData$time, status = RawData$status)]

RawData2CP <- survival::tmerge(
  as.data.frame(RawData[, .(id = .I, TAJ, AGEindex, SEX, NONCARDIOVASCULARDISEASES, CARDIOVASCULARDISEASES,
                            METABOLICDISEASES, PreindexGFR, MINORAMP, MAJORAMP, VASCPROCEDURES, MEDICATION,
                            TimeToIndexEvent, PROCEDURE, PreindexMALIGN, PreindexCHD, PreindexCVD)]),
  as.data.frame(EventTimes), id = id, status = event(time, status))
RawData2CP <- survival::tmerge(
  RawData2CP, as.data.frame(EventTimes), id = id,
  PostIndeMalign = tdc(PostIndeMalign), PostIndexAMI = tdc(PostIndexAMI), PostindexSTROKE = tdc(PostindexSTROKE),
  PostIndexMinorampipszi = tdc(PostIndexMinorampipszi), PostIndexMinorampkontra = tdc(PostIndexMinorampkontra),
  PostIndexMAJORAMPipszi = tdc(PostIndexMAJORAMPipszi), PostIndexMAJORAMPkontra = tdc(PostIndexMAJORAMPkontra),
  POSTINDEXOPENipszi = tdc(POSTINDEXOPENipszi), POSTINDEXOPENkontra = tdc(POSTINDEXOPENkontra),
  POSTINDEXENDOipszi = tdc(POSTINDEXENDOipszi), POSTINDEXENDOkontra = tdc(POSTINDEXENDOkontra))
RawData2CP <- as.data.table(RawData2CP)

for (v in names(PostIndexEventDates)) {
  RawData2CP[id %in% which(is.na(RawData[[v]]) | (RawData[[v]] == 1 & is.na(RawData[[PostIndexEventDates[v]]]))),
             (v) := NA]
}

RawData2CP[, `:=`(
  NONCARDIOVASCULARDISEASES = NONCARDIOVASCULARDISEASES + (PreindexMALIGN == 0 & PostIndeMalign == 1),
  CARDIOVASCULARDISEASES = CARDIOVASCULARDISEASES + (PreindexCHD == 0 & PostIndexAMI == 1) +
    (PreindexCVD == 0 & PostindexSTROKE == 1),
  MINORAMP = factor(as.numeric(MINORAMP == "1" | PostIndexMinorampipszi == 1 | PostIndexMinorampkontra == 1),
                    levels = c(0, 1)),
  MAJORAMP = MAJORAMP + PostIndexMAJORAMPipszi + PostIndexMAJORAMPkontra,
  VASCPROCEDURES = VASCPROCEDURES + POSTINDEXOPENipszi + POSTINDEXOPENkontra + POSTINDEXENDOipszi +
    POSTINDEXENDOkontra)]

RawData2CP <- na.omit(RawData2CP[
  , .(id, TAJ, tstart, tstop, status, AGEindex, SEX, NONCARDIOVASCULARDISEASES, CARDIOVASCULARDISEASES,
      METABOLICDISEASES, PreindexGFR, MINORAMP, MAJORAMP, VASCPROCEDURES, MEDICATION, TimeToIndexEvent, PROCEDURE)])

################################################################

dd <- datadist(RawData2CP)
options(datadist = "dd")

sum(RawData2CP$status)/15

fit <- cph(Surv(tstart, tstop, status) ~ rcs(AGEindex, 3) + SEX + scored(NONCARDIOVASCULARDISEASES) +
             rcs(CARDIOVASCULARDISEASES, 3) + scored(METABOLICDISEASES) + rcs(PreindexGFR, 3) +
             MINORAMP + scored(MAJORAMP) + rcs(VASCPROCEDURES, 3) + scored(MEDICATION) +
             rcs(TimeToIndexEvent, 3) + scored(PROCEDURE), data = RawData2CP, x = TRUE, y = TRUE)
fit <- robcov(fit, cluster = RawData2CP$TAJ)

length(coef(fit))

anova(fit)

fit <- cph(Surv(tstart, tstop, status) ~ rcs(AGEindex) + SEX + NONCARDIOVASCULARDISEASES +
             CARDIOVASCULARDISEASES + METABOLICDISEASES + PreindexGFR +
             MINORAMP + MAJORAMP + VASCPROCEDURES + MEDICATION +
             rcs(TimeToIndexEvent) + PROCEDURE , data = RawData2CP, x = TRUE, y = TRUE)
fit <- robcov(fit, cluster = RawData2CP$TAJ)

length(coef(fit))

anova(fit)

fit <- cph(Surv(tstart, tstop, status) ~ AGEindex + SEX + NONCARDIOVASCULARDISEASES +
             CARDIOVASCULARDISEASES + METABOLICDISEASES + PreindexGFR +
             MINORAMP + MAJORAMP + VASCPROCEDURES + MEDICATION +
             rcs(TimeToIndexEvent) + PROCEDURE, data = RawData2CP, x = TRUE, y = TRUE)
fit <- robcov(fit, cluster = RawData2CP$TAJ)

length(coef(fit))

anova(fit)

set.seed(1)
validate(fit, B = 1000, cluster = RawData2CP$TAJ)

survival::cox.zph(fit)

cairo_pdf("Modell_Mortalitas.pdf", width = 16, height = 9)
plot(Predict(fit, TimeToIndexEvent), anova = anova(fit), pval = TRUE, xlab = "Time from first detection of PAD to revascularization [year]")
plot(Predict(fit), anova = anova(fit), pval = TRUE)
dev.off()

tiff("SupplementaryFigureS1.tiff", width = 16, height = 9, units = "in", res = 300, compression = "lzw")
plot(Predict(fit, TimeToIndexEvent), anova = anova(fit), pval = TRUE, xlab = "Time from first detection of PAD to revascularization [year]", adj.subtitle = FALSE)
dev.off()

tab <- data.table(variable = names(coef(fit)), b = coef(fit), se = sqrt(diag(vcov(fit))))[
  !startsWith(variable, "TimeToIndexEvent"),
  .(variable, HR = exp(b), lower = exp(b - qnorm(0.975) * se), upper = exp(b + qnorm(0.975) * se),
    p = 2 * pnorm(-abs(b / se)))]

fwrite(tab, "Modell_Mortalitas.csv", dec = ",", sep =";", bom = TRUE)

################################################################

cif <- rbind(
  with(survival::survfit(survival::Surv(timeMinorAmp, ifelse(statusMinorAmp == "Minor amputation", 1, 0)) ~ 1,
                         data = RawData2, cluster = TAJ),
       data.table(time, outcome = "Minor amputation (KM)", estimate = 1 - surv, conf.low = 1- upper, conf.high = 1 - lower)),
  with(survival::survfit(survival::Surv(timeMinorAmp, as.factor(statusMinorAmp)) ~ 1,
                         data = RawData2, cluster = TAJ),
       data.table(time, outcome = "Minor amputation", estimate = pstate[,"Minor amputation"],
                  conf.low = lower[,"Minor amputation"], conf.high = upper[,"Minor amputation"])),
  with(survival::survfit(survival::Surv(timeMajorAmp, ifelse(statusMajorAmp == "Major amputation", 1, 0)) ~ 1,
                         data = RawData2, cluster = TAJ),
       data.table(time, outcome = "Major amputation (KM)", estimate = 1 - surv, conf.low = 1- upper, conf.high = 1 - lower)),
  with(survival::survfit(survival::Surv(timeMajorAmp, as.factor(statusMajorAmp)) ~ 1,
                         data = RawData2, cluster = TAJ),
       data.table(time, outcome = "Major amputation", estimate = pstate[,"Major amputation"],
                  conf.low = lower[,"Major amputation"], conf.high = upper[,"Major amputation"])),
  with(survival::survfit(survival::Surv(timeMajorAmp, ifelse(statusMajorAmp == "Censoring", 0, 1)) ~ 1,
                         data = RawData2, cluster = TAJ),
       data.table(time, outcome = "Major amputation or death", estimate = 1 - surv, conf.low = 1- upper, conf.high = 1 - lower)),
  with(survival::survfit(survival::Surv(timeMinorAmp, ifelse(statusMinorAmp == "Censoring", 0, 1)) ~ 1,
                         data = RawData2, cluster = TAJ),
       data.table(time, outcome = "Minor amputation or major amputation or death",
                  estimate = 1 - surv, conf.low = 1- upper, conf.high = 1 - lower)),
  with(survival::survfit(survival::Surv(time, status) ~ 1,
                         data = RawData2, cluster = TAJ),
       data.table(time, outcome = "Death", estimate = 1 - surv, conf.low = 1- upper, conf.high = 1 - lower))
)

fwrite(
  dcast(melt(cif, id.vars = c("time", "outcome"))[, (with(approx(time, value, c(30, 90, 6*30, 365, 3*365)/365), list(time = x * 365, value = y))),
                                                  .(outcome, variable)],
        outcome + time ~ variable, value.var = "value"),
  "Survival_Univariate.csv", dec = ",", sep =";", bom = TRUE)

ggplot(cif[outcome %in% c("Major amputation or death", "Minor amputation",
                          "Major amputation", "Death")],
       aes(x = time, y = estimate*100, ymin = conf.low*100, ymax = conf.high*100,
           color = outcome, fill = outcome, group = outcome)) + geom_line() +
  geom_ribbon(alpha = 0.2, color = NA) + theme(legend.position = "bottom") +
  labs(x = "Time [year]", y = "Cumulative incidence [%]", color = "", fill = "") +
  theme(text = element_text(family = "Garamond"))
ggsave("Survival_Univariate.pdf", width = 16/3*2, height = 9/3*2, device = cairo_pdf)
ggsave("Figure1.tiff", width = 16/3*2, height = 9/3*2, compression = "lzw")

################################################################

RawData2CPMajorAmp <- merge(
  RawData2CP[, !c("status", "MAJORAMP")],
  RawData[, .(id = .I, statusMajorAmp, timeMajorAmp = fifelse(timeMajorAmp == 0, 0.5/365.241, timeMajorAmp))],
  by = "id")[tstart < timeMajorAmp]
RawData2CPMajorAmp[, statusMajorAmp := factor(fifelse(tstop >= timeMajorAmp, as.character(statusMajorAmp), "Censoring"),
                                              levels = levels(statusMajorAmp))]
RawData2CPMajorAmp[, `:=`(tstop = pmin(tstop, timeMajorAmp), timeMajorAmp = NULL)]

fit2 <- cph(Surv(tstart, tstop, statusMajorAmp == "Major amputation") ~ rcs(AGEindex, 3) + SEX + scored(NONCARDIOVASCULARDISEASES) +
              rcs(CARDIOVASCULARDISEASES, 3) + scored(METABOLICDISEASES) + rcs(PreindexGFR, 3) +
              MINORAMP + rcs(VASCPROCEDURES, 3) + scored(MEDICATION) +
              rcs(TimeToIndexEvent, 3) + scored(PROCEDURE),
            data = RawData2CPMajorAmp, x = TRUE, y = TRUE)
fit2 <- robcov(fit2, cluster = RawData2CPMajorAmp$TAJ)

anova(fit2)

fit2 <- cph(Surv(tstart, tstop, statusMajorAmp == "Major amputation") ~ AGEindex + SEX + NONCARDIOVASCULARDISEASES +
              CARDIOVASCULARDISEASES + METABOLICDISEASES + PreindexGFR +
              MINORAMP + VASCPROCEDURES + scored(MEDICATION) +
              TimeToIndexEvent + PROCEDURE,
            data = RawData2CPMajorAmp, x = TRUE, y = TRUE)
fit2 <- robcov(fit2, cluster = RawData2CPMajorAmp$TAJ)

anova(fit2)

fit2 <- cph(Surv(tstart, tstop, statusMajorAmp == "Major amputation") ~ AGEindex + SEX + NONCARDIOVASCULARDISEASES +
              CARDIOVASCULARDISEASES + METABOLICDISEASES + PreindexGFR +
              MINORAMP + VASCPROCEDURES + MEDICATION +
              TimeToIndexEvent + PROCEDURE,
            data = RawData2CPMajorAmp, x = TRUE, y = TRUE)
fit2 <- robcov(fit2, cluster = RawData2CPMajorAmp$TAJ)

anova(fit2)

cairo_pdf("Modell_Majormap.pdf", width = 16, height = 9)
plot(Predict(fit2), anova = anova(fit2), pval = TRUE)
dev.off()

tab2 <- data.table(variable = names(coef(fit2)), b = coef(fit2), se = sqrt(diag(vcov(fit2))))[
  , .(variable, HR = exp(b), lower = exp(b - qnorm(0.975) * se), upper = exp(b + qnorm(0.975) * se),
      p = 2 * pnorm(-abs(b / se)))]

fwrite(tab2, "Modell_Majormap.csv", dec = ",", sep =";", bom = TRUE)

################################################################

fitdata <- rbind(
  cbind(tab, Type = "All-cause mortality"),
  cbind(tab2, Type = "Ipsilateral lower limb major amputation-free survival")
)

varnames <- c("Age", "Sex (female)", "Cardiovascular diseases", "Non-cardiovascular diseases",
              "Metabolic diseases", "GFR", "Lower limb minor amputations", "Lower limb major amputations",
              "Lower limb revascularisations", "Preventive medications", "Conditions of index procedure",
              "Time of PAD detection to index event")
fitdata <- merge(fitdata, data.table(
  variable = c("AGEindex", "SEX=2", "CARDIOVASCULARDISEASES", "NONCARDIOVASCULARDISEASES", "METABOLICDISEASES", "PreindexGFR",
               "MINORAMP=1", "MAJORAMP", "VASCPROCEDURES", "MEDICATION", "PROCEDURE", "TimeToIndexEvent"),
  varname = varnames), by = "variable")
fitdata[, varname := factor(varname, levels = varnames)]
fitdata <- fitdata[order(varname)]

ggplot(fitdata, aes(y = varname, x = HR, xmin = lower, xmax = upper, group = Type, color = Type)) +
  geom_point(position = position_dodge(width = 0.5)) +
  geom_errorbar(position = position_dodge(width = 0.5), width = 0.5) +
  geom_vline(xintercept = 1) + scale_y_discrete(limits = rev) +
  labs(x = "Hazard ratio", y = "") +
  theme(legend.position = "bottom", legend.title = element_blank(),
        text = element_text(family = "Garamond"))
ggsave("Modells.pdf", width = 16/3*2, height = 9/3*2, device = cairo_pdf)
ggsave("Figure2.tiff", width = 16/3*2, height = 9/3*2, compression = "lzw")

fitdata$formatted <- paste0(round(fitdata$HR, 3), " (", round(fitdata$lower, 3), " - ",
                            round(fitdata$upper, 3), ") p = ", round(fitdata$p, 4))
fwrite(dcast(fitdata, variable ~ Type, value.var = "formatted"), "Modells.csv",
       dec = ",", sep =";", bom = TRUE)