library(data.table)
library(survival)
library(ggplot2)
library(plotly)

colour_code = c("#f98b83", 
                "#a9ab00", 
                "#51d3a6",
                "#00b4fc", 
                "#f471ff", 
                "#808080")
colour_code2 = c("#f98b83", 
                "#a9ab00", 
                "red", 
                "#51d3a6",
                "#00b4fc", 
                "#f471ff", 
                "#808080")

predf = fread("./MEGA_bp_t2d_phenofile_for_survival_curves.txt.gz", na = c("NA","")) #Pheno file + partitioned PGS extracted from UKB
cov = fread("./MEGA_bp_t2d_phenofile.tsv.gz") #Pheno file + partitioned PGS extracted from UKB
#hypothesis: columns order is conserved
df = cbind(predf, cov[, c("age", "sex", "array", "PC1", "PC2", "PC3", "PC4", "PC5", "PC6")])

# Assume a start date for calculating the time (e.g., date of birth or study start date)
start_date <- as.Date(min(pmin(df$overall_hypertension_date, df$overall_T2D_date, na.rm = TRUE), na.rm = T))

# Calculate time to comorbidity and status THIS IS THE CRUCIAL STEP, A GUIDE IS HERE:
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC10957029/
# df$comorbidity_time <- as.Date(pmax(df$overall_hypertension_date, df$overall_T2D_date, na.rm = TRUE))
# df[is.na(comorbidity_time), comorbidity_time := as.Date(max(comorbidity_time))]
# df[is.infinite(comorbidity_time), comorbidity_time := 0]
df$comorbidity_status <- ifelse(!is.na(df$overall_hypertension_date) & !is.na(df$overall_T2D_date), 1, 0)
df[, comorbid_time := as.numeric(difftime(overall_hypertension_date, overall_T2D_date, units = "days"))/365.25]
df[, comorbid_time2 := abs(comorbid_time)]
#Second, you need to have a follow-up time duration variable which is the time between the start of the study and the first event, or the time between the start of the study and the end of the study for participants who did not have an event. In other words, this is the duration that a participant was “at risk” for the event.
df[is.na(comorbid_time), comorbid_time2 := max(df$comorbid_time2, na.rm = T)]

cumulative_hazard <- data.frame()
concord <- data.frame()

#for top10_cluster2 (hard clustering, see makesuperpheno)
cat("\nReady for assoc per pgs groups")
for (i in unique(df$top10_cluster2)) {

  subdf = df[top10_cluster2 == i, ]

  cox_model <- coxph(Surv(comorbid_time2, comorbidity_status) ~ 1 + age + sex + array + PC1 + PC2 + PC3 + PC4 + PC5 + PC6, data = subdf)
  survival_object <- survfit(cox_model, newdata = subdf)

  # Extract cumulative hazard and CI
  cum_hazard <- summary(survival_object, times = unique(subdf$comorbid_time2))
  cum <- data.frame(time = cum_hazard$time,
                    hazard = cum_hazard$cumhaz,
                    lower = -log(cum_hazard$lower),
                    upper = -log(cum_hazard$upper),
                    cluster = rep(i, length(cum_hazard$time)))

  # Combine results
  cumulative_hazard <- rbind(cumulative_hazard, cum)
  
  # Combine concordance
  con <- concordance(cox_model, newdata = subdf)
  concord <- rbind(concord, data.frame(cluster = i,
                              c_index = con$concordance, 
                              c_index_var = con$var))
  
  cat("\n\tGroup", i, "done - out of", length(unique(df$top10_cluster2)), "test")

}

p = ggplot(data = cumulative_hazard, aes(x = time, y = hazard, color = factor(cluster))) +
  geom_step(size = 1) +
  labs(x = "Follow-up (years)", y = "Cumulative Hazard",
       color = paste0("Being in the top\n", 33, "% of PGS\ndistribution")) +
  xlim(0, 25) +
  theme_minimal()
ggsave(p, filename = "./pgs_after_clustering/wcovar_cox_haz_handmade_top10_cluster2_0.67.png", height = 8, width = 12, dpi = 400)

pint = ggplotly(p)
htmlwidgets::saveWidget(pint, file = "./pgs_after_clustering/wcovar_cox_haz_handmade_top10_cluster3_0.67.html")

pp = ggplot(data = cumulative_hazard, aes(x = time, y = hazard, color = factor(cluster))) +
  geom_step(size = 1) +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = factor(cluster)), alpha = 0.1) +
  labs(x = "Follow-up (years)", y = "Cumulative Hazard", 
       color = paste0("Being in the top\n", 33, "% of PGS\ndistribution"),
       fill = paste0("Being in the top\n", 33, "% of PGS\ndistribution")) +
  scale_colour_manual(values = colour_code2, labels = group_labels2) +
  scale_fill_manual(values = colour_code2, labels = group_labels2) +
  xlim(0, 25) +
  theme_minimal()
ggsave(pp, filename = "./pgs_after_clustering/wcovar_cox_haz_handmade_top10_cluster2_0.67_wCI.png", height = 8, width = 12, dpi = 400)
