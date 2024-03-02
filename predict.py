import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import tensorflow as tf
import warnings

from sklearn.exceptions import ConvergenceWarning
from sklearn.linear_model import LinearRegression, Lasso, Ridge
from sklearn.metrics import mean_squared_error, r2_score
from sklearn.model_selection import train_test_split, GridSearchCV

# Ignore ConvergenceWarning
warnings.simplefilter('ignore', ConvergenceWarning)

# Number of snps to use
N = [20,40,60,80,100]

# Data files
snp_file = 'data/emmax_out/EMMAX.0_5_FT10.top'
geno_file = 'data/filtered.coded_call_method_54.tair9.FT10.csv'
pheno_file = 'data/filtered.FT10.txt'


def read_data(geno_file, pheno_file, snp_file):
    # Read in data
    colnames = ['Chromosome', 'Positions', 'pvalue', 'AdjustedR2', 'coefficient', 'Sd_Err', 'MAF_count']
    snps = pd.read_csv(snp_file, skiprows=2, names=colnames)
    geno = pd.read_csv(geno_file)
    pheno = pd.read_csv(pheno_file, index_col=0)

    # Sort snps by p-value
    snps = snps.sort_values(by='pvalue')
    # snps.head(100).to_csv('data/ssnps_jawamix5.csv')

    return geno, pheno, snps


def create_train_test(geno, pheno, snps, n):
    # Keep only top n snps
    top_n_snps = snps.head(n)

    # Filter genotypes for top_n_snps
    geno_filtered = pd.merge(geno, top_n_snps, how='inner', on=['Chromosome', 'Positions'])
    geno_filtered = geno_filtered[geno.columns] # Drop extra columns from top_n_snps
    geno_filtered = geno_filtered.drop(labels=['Chromosome','Positions'], axis=1) # Drop chr and pos columns

    # Separate into training and test
    geno_filtered = geno_filtered.astype(float)
    geno_filtered.index = geno_filtered.index.astype('int64')
    X_train, X_test, y_train, y_test = train_test_split(geno_filtered.T, pheno, test_size=0.2, random_state=42)
    return X_train, X_test, y_train, y_test


def plot(y_test, y_pred, model_name):
    plt.figure(figsize=(10, 6))
    plt.scatter(y_test, y_pred, alpha=0.5)
    plt.xlabel('Actual phenotype')
    plt.ylabel('Predicted phenotype')
    plt.title(f'Actual vs. predicted phenotypes using {model_name}')
    plt.plot([y_test.min(), y_test.max()], [y_test.min(), y_test.max()], 'k--', lw=2)
    plt.show()


def run_model(model_func, model_name, alphas=None):
    geno, pheno, snps = read_data(geno_file, pheno_file, snp_file)
    for n in N:
        X_train, X_test, y_train, y_test = create_train_test(geno,pheno,snps,n)
        if alphas is not None:
            model = GridSearchCV(estimator=model_func(), param_grid=dict(alpha=alphas), cv=5, scoring='neg_mean_squared_error')
            model.fit(X_train, y_train)
            best_alpha = model.best_estimator_.alpha
            print(f"Best alpha for {model_name} with {n} SNPs: {best_alpha}")
            model = model_func(alpha=best_alpha)
        else:
            model = model_func()
        model.fit(X_train, y_train)
        y_pred = model.predict(X_test)
        mse = mean_squared_error(y_test, y_pred)
        r2 = r2_score(y_test, y_pred)
        print(f"MSE using {model_name} with {n} SNPs: {mse}")
        print(f"R^2 using {model_name} with {n} SNPs: {r2}")
        plot(y_test, y_pred, f"{model_name} ({n} SNPs)")


# Linear Regression
run_model(LinearRegression, "Linear Regression")

# Lasso
alphas = np.logspace(-4, 2, num=50)
run_model(Lasso, "Lasso", alphas)

# Ridge
run_model(Ridge, "Ridge", alphas)










# Implementation using tensorflow
#
# # Create the model
# model = tf.keras.Sequential([
#     tf.keras.layers.Dense(1)
# ])
#
# # Compile the model
# model.compile(loss=tf.keras.losses.mse,
#                 optimizer=tf.keras.optimizers.legacy.SGD(learning_rate=0.01),
#                 metrics=["mse"])
#
# # Fit the model
# fit = model.fit(X_train, y_train, epochs=100)
#
# # Plot the model fit
# pd.DataFrame(fit.history).plot()
# plt.xlabel('epochs')
# plt.ylabel('loss')
# plt.show()
#
# # Evaluate the model
# model.evaluate(X_test, y_test)