# 📊 Logistic Regression Estimators under Multicollinearity

A research-oriented statistical modeling project that investigates the performance of several shrinkage and biased estimators for multicollinear logistic regression models using simulation studies and real-world datasets.

---

## 🚀 Features

- Implementation of advanced logistic regression estimators
- Multicollinearity diagnostics and analysis
- Comparative model evaluation
- Simulation studies for estimator performance
- Real-world medical dataset applications
- Classification performance analysis

---

## 📚 Implemented Estimators

- Maximum Likelihood Estimator (MLE)
- Ridge Logistic Estimator
- Liu Logistic Estimator
- Almost Unbiased Ridge Estimator (AURE)
- Almost Unbiased Liu Estimator (AULE)
- Adjusted Logistic Liu Estimator (ALLE)
- Kibria–Lukman Logistic Estimator (LKLE)
- Logistic Dorugade Estimator (LDE)
- Logistic Modified Ridge-Type Estimator (LMRT)
- Logistic Stein Ridge Estimator (LSRE)
- Logistic Stein Kibria Estimator (LSKE)

---

## 🧮 Logistic Regression Model

The project studies the logistic regression model:

\[
P(Y=1|X)=\frac{e^{X\beta}}{1+e^{X\beta}}
\]

where shrinkage estimators are used to reduce instability caused by multicollinearity among predictors.

---

## 📊 Evaluation Metrics

Models are compared using:

- Accuracy
- Sensitivity
- Specificity
- Precision
- F1-Score
- Balanced Accuracy
- ROC Curve
- AUC
- Log-Loss
- Brier Score

---

## 🔬 Dataset Applications

- Diabetes Dataset
- Prostate Cancer Dataset
- Simulated Multicollinear Data

---

## 🛠️ Technologies Used

- **R Programming**
- Logistic Regression
- Matrix Algebra
- Statistical Modeling
- Multicollinearity Analysis

### Packages

```r
MASS
pROC
glmnet
readxl
caret
```

---

## 📂 Project Structure

```bash
├── data/
│   ├── diabetes.csv
│   ├── prostate_cancer_data.xlsx
│
├── scripts/
│   ├── estimators.R
│   ├── simulation_study.R
│   ├── model_comparison.R
│   ├── evaluation_metrics.R
│
├── results/
│   ├── accuracy_tables.csv
│   ├── roc_curves.png
│   ├── coefficient_comparison.csv
│
├── README.md
```

---

## 📈 Key Findings

- Ridge-based estimators improved prediction stability under multicollinearity.
- LMRT, Ridge, and LDE achieved higher classification accuracy compared to MLE.
- LSRE produced the highest AUC among evaluated estimators.
- Shrinkage methods reduced variance while maintaining strong predictive performance.

---

## ▶️ How to Run

```r
# Load dataset
data <- read.csv("data/diabetes.csv")

# Source estimator functions
source("scripts/estimators.R")

# Run model comparison
source("scripts/model_comparison.R")
```

---

## 👨‍💻 Author

### Aditya Pardeshi
M.Sc. Statistics

Interested in:
- Statistical Modeling
- Data Science
- Machine Learning
- Biostatistics
- Regression Analysis

---

## ⭐ Future Improvements

- Cross-validation based tuning parameter selection
- Bayesian shrinkage estimators
- R Shiny interactive dashboard
- Automated estimator ranking system
- High-dimensional logistic regression extensions

---

## 📜 License

This project is intended for academic and research purposes.
