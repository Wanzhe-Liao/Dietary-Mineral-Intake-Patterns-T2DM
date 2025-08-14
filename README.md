# Dietary-Mineral-Intake-Patterns-T2DM
## Project Overview

This is a comprehensive NHANES-based epidemiological study analyzing dietary mineral intake patterns in 7,909 Type 2 Diabetes patients. The project uses hierarchical clustering on 9 mineral proportions to identify distinct dietary patterns and examines their associations with mortality, cardiovascular disease, and cancer outcomes.

**Key Study Design:**

- Population: NHANES 1999-2018 T2DM patients (≥18 years)
- Primary Analysis: 4-cluster mineral intake patterns (hierarchical clustering with Ward.D2)
- Outcomes: All-cause mortality, cause-specific mortality, life expectancy
- Follow-up: Complete follow-up through 2019 (median ~8 years)

## Core Data Structure

**Primary Dataset:** `DM_ALL.RData` contains 7,909 patients with:

- Demographic/clinical variables (SEQN_new, RIDAGEYR, RIAGENDR, etc.)
- Mineral intake data (DRXTCALC through DRXTSELE) - 9 minerals total
- Survival data (MORTSTAT, PERMTH_INT, VNUCOD_LEADING)
- Cluster assignment (clusters: 1-4)
- Survey weights (WTDRD1) - **ALWAYS use for population-representative analyses**

## Analysis Pipeline Architecture

### 1. Data Preprocessing (`NHANES.R`)

- Multi-source data merging (demographics, dietary, mortality, questionnaire)
- T2DM case definition: physician diagnosis OR medication use OR HbA1c ≥6.5% OR FBG ≥126 mg/dL
- Missing data imputation using missForest for categorical variables (see lines 398-404)
- Variable recoding and standardization

### 2. Clustering Analysis (`hierarchical clustering.R`)

**Optimal clustering determination:**

- PCA dimensionality reduction (first 3 PCs explain 85% variance)
- Multiple validation methods: elbow, BIC (Mclust), gap statistic
- UMAP/t-SNE visualization for cluster validation

Key clustering workflow:

```r
# Mineral proportion calculation and scaling
scaled_data <- element_pro[,Ccol]  # Ccol contains 9 mineral columns
dist_matrix <- dist(scaled_data, method = "euclidean")
hc_model <- hclust(dist_matrix, method = "ward.D2")
clusters <- cutree(hc_model, k = 4)  # Optimal k=4 determined by multiple criteria
```

### 3. Primary Survival Analyses

**Overall Survival (`survival_analysis.R`):**

- Model 1: Basic adjustment (age, sex, race)
- Model 2: Full adjustment (+ BMI, BP, HbA1c, comorbidities, medications)
- Generates both full follow-up and 6-year restricted analyses
- Each cluster vs. all others comparison approach

**Cause-Specific Mortality (`cause_specific_mortality.R`):**

- Fine-Gray competing risks models for 10 major death causes
- Automated plot generation for all cluster-cause combinations
- Model 1/Model 2 adjustment schemes parallel to overall survival

**Life Expectancy (`life_expectancy_analysis.R`):**

- Age-as-timescale survival models
- Conditional life expectancy calculations from age 45-100
- Between-cluster differences in years of life lost/gained

### 4. Alternative Intake Pattern Analyses

**Tertile-based Analysis (`mineral_intake_survival.R`):**

- High/Medium/Low total mineral intake groups
- Parallel analytical framework to cluster-based approach
- Cross-validation of clustering results

**Correlation Analysis (`mineral_correlation_analysis.R`):**

- Spearman correlations between all 9 minerals
- Log-transformation for improved linearity
- Comprehensive visualization suite (heatmaps, scatterplot matrices)

### 5. Baseline Characteristics (`DM/baseline_table_gtsummary.R`)

Uses gtsummary package with proper survey weight handling for publication-ready tables.

### 6. Mineral Intake Range Analysis (`comprehensive_mineral_ranges_comparison.R`)

**Optimized Clinical Target Ranges for Cluster 3 (Optimal Pattern):**

- **Three Range Types**: 95% percentile (2.5%-97.5%), IQR (25%-75%), Optimized (IQR+RDA/UL integrated)
- **Clinical Priority System**: HIGH (Ca, Mg, Zn, K), MEDIUM (P, Fe, Se), CAUTION (Na), LOW (Cu)
- **Safety Integration**: Lower bounds ≥80% RDA, Upper bounds ≤90% UL where established
- **T2DM-Specific Optimization**: Emphasizes minerals critical for glucose control and CVD protection

**Key Scripts:**

- `high_mineral_cluster3_mineral_ranges.R` - Complete analysis with all statistical measures
- `optimized_mineral_ranges_display.R` - Streamlined clinical target display
- `comprehensive_mineral_ranges_comparison.R` - Three-method comparison table

## Important Analysis Conventions

### Survey Weights

**CRITICAL:** Always use `weights = WTDRD1` in survival models and `survey` package functions for population representativeness.

### Clustering Interpretation

- Cluster 1: Low overall mineral intake pattern
- Cluster 2: Moderate balanced intake
- Cluster 3: High balanced intake (reference/protective pattern)
- Cluster 4: High sodium, low other minerals

### Statistical Models

**Cox Proportional Hazards Standard Format:**

```r
# Basic template for survival analysis
coxph(Surv(PERMTH_INT, MORTSTAT) ~ cluster_factor + covariates, 
      data = data, weights = WTDRD1)
```

**Cause-specific mortality:** Use competing risks (Fine-Gray models) when analyzing specific death causes.

## File Organization

- `/DM/` - Main analysis scripts and results
- `/survival_plots/` - Kaplan-Meier curves and forest plots
- `/cause_specific_plots/` - Death cause analysis visualizations
- `/subgroup_plots/` - Stratified analysis results
- `/life_expectancy_results/` - Life expectancy comparison outputs
- `/mineral_profile/` - Cluster 3 mineral range analyses and clinical guidance tables

## Key Variables Reference

**Minerals (proportion-based):** DRXTCALC (Calcium), DRXTPHOS (Phosphorus), DRXTMAGN (Magnesium), DRXTIRON (Iron), DRXTZINC (Zinc), DRXTCOPP (Copper), DRXTSODI (Sodium), DRXTPOTA (Potassium), DRXTSELE (Selenium)

**Primary Outcomes:** MORTSTAT (0=alive, 1=dead), PERMTH_INT (follow-up months)

**CVD Composite:** MCQ160B+MCQ160C+MCQ160D+MCQ160E+MCQ160F (heart failure, CHD, angina, MI, stroke)

**Key Covariates:** LBXGH (HbA1c), VNAVEBPXSY/VNLBAVEBPXDI (BP), alcohol_week, hyp_med, hbc_med

## Standard Analysis Commands

**Load main dataset:**

```r
load("DM_ALL.RData")  # Primary dataset with 7,909 T2DM patients
```

**Complete analysis pipeline (run in sequence):**

```r
# 1. Data preprocessing and clustering (if needed)
source("NHANES.R")  # Raw data processing and T2DM case definition
source("hierarchical clustering.R")  # Generate 4-cluster solution

# 2. Primary analyses
source("survival_analysis.R")  # Overall survival analysis (Model 1 & 2)
source("cause_specific_mortality.R")  # Death cause-specific analysis
source("life_expectancy_analysis.R")  # Life expectancy calculations

# 3. Baseline characteristics and correlations
source("DM/baseline_table_gtsummary.R")  # Publication-ready Table 1
source("mineral_correlation_analysis.R")  # Mineral intake correlations

# 4. Alternative intake pattern analyses
source("mineral_intake_survival.R")  # Tertile-based mineral intake analysis
source("mineral_intake_subgroup_analysis.R")  # Subgroup analyses
source("mineral_intake_cause_specific_mortality.R")  # Cause-specific by intake levels

# 5. Clinical mineral range analysis
source("high_mineral_cluster3_mineral_ranges.R")  # Complete Cluster 3 range analysis
source("comprehensive_mineral_ranges_comparison.R")  # Display all three range types
```

**Quick validation checks:**

```r
# Verify data integrity
table(DM_ALL$clusters)  # Should show 4 clusters with reasonable distributions
summary(DM_ALL$MORTSTAT)  # Death events
range(DM_ALL$PERMTH_INT, na.rm=TRUE)  # Follow-up time range
```

## Specialized Analysis Features

### Subgroup Analyses

All primary analyses include automatic stratification by:

- Age (<60 vs ≥60 years)
- Sex (male vs female) 
- BMI categories (<25, 25-29.9, ≥30)
- HbA1c control (<7.0% vs ≥7.0%)
- Medication patterns (no meds, OAD only, insulin ±OAD)
- CVD status and hypertension medication use

### Quality Control Checks

Scripts include built-in validation:

```r
# Example from cause_specific_mortality.R (lines 29-40)
required_vars <- c("PERMTH_INT", "MORTSTAT", "VNUCOD_LEADING", "clusters", "WTDRD1")
missing_vars <- required_vars[!required_vars %in% names(DM_ALL)]
death_count <- sum(DM_ALL$MORTSTAT == 1, na.rm = TRUE)
if (death_count < 50) warning("Insufficient death events for reliable analysis")
```

### Parallel Analysis Approaches

The codebase implements both:

1. **Cluster-based analyses**: Using 4-group hierarchical clustering patterns
2. **Tertile-based analyses**: Using total mineral intake tertiles (High/Medium/Low)
   - Cross-validation approach ensuring robustness of findings
   - Alternative grouping method for sensitivity analyses

## Output Standards

- All plots use 300 DPI for publication quality
- Forest plots include 95% CIs and p-values with proper positioning
- Results tables saved as both CSV and RData formats
- Survival curves show both Model 1 and Model 2 adjusted versions
- Always generate both full follow-up and 6-year restricted analyses
- Automated directory creation for organized output storage

## Critical Notes

- **Survey weights mandatory**: Population weights (WTDRD1) essential for external validity to US adults with diabetes
- **Complete-case analysis**: No imputation for survival outcomes (only for baseline covariates)
- **Competing risks**: All cause-specific mortality analyses use Fine-Gray models
- **Model consistency**: Model 1 (basic) and Model 2 (full) adjustment schemes used consistently across all analyses
- **Methodological rigor**: This is defensive epidemiological research following established standards

## Clinical Translation: Optimized Mineral Target Ranges

**Based on Cluster 3 Analysis (n=1,849, 23.4% of T2DM cohort):**

**HIGH Priority Interventions (Focus areas for T2DM management):**

- **Calcium**: 800-1170 mg/day (bone health, CVD protection)
- **Magnesium**: 256-315 mg/day (insulin sensitivity, glucose control)
- **Zinc**: 6.4-13.4 mg/day (immune function, glycemic control)  
- **Potassium**: 2800-3414 mg/day (BP control, CVD protection)

**CAUTION (Strict monitoring required):**

- **Sodium**: 1703-2070 mg/day (limit <2300mg per guidelines)

**MEDIUM Priority (Monitor and maintain adequate intake):**

- **Phosphorus**: 829-1513 mg/day (bone formation)
- **Iron**: 8.2-18.2 mg/day (avoid excess due to oxidative stress)
- **Selenium**: 59.6-113 mcg/day (antioxidant function)

**Range Calculation Method:**

- **Optimized ranges** = IQR (25th-75th percentile) integrated with RDA/UL safety guidelines
- **Lower bounds**: max(IQR_25th, 80% of RDA)
- **Upper bounds**: min(IQR_75th, 90% of UL where established)
- **More clinically relevant than wide 95% percentile ranges**

**Clinical Application:**
Use these evidence-based targets for individual T2DM patient nutrition counseling, focusing on HIGH priority minerals first while ensuring sodium control.
