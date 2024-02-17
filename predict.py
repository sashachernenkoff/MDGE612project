import matplotlib.pyplot as plt
import pandas as pd
import tensorflow as tf

from sklearn.model_selection import train_test_split
from sklearn.linear_model import LinearRegression, Lasso, Ridge
from sklearn.metrics import mean_squared_error, r2_score


# Number of snps to use
N = 20


def plot(y_test, y_pred):
    plt.figure(figsize=(10, 6))
    plt.scatter(y_test, y_pred, alpha=0.5)
    plt.xlabel('Actual phenotype')
    plt.ylabel('Predicted phenotype')
    plt.title('Actual vs. Predicted Phenotypes')
    plt.plot([y_test.min(), y_test.max()], [y_test.min(), y_test.max()], 'k--', lw=2)
    plt.show()


# Read in data
colnames = ['Chromosome', 'Positions', 'pvalue', 'AdjustedR2', 'coefficient', 'Sd_Err', 'MAF_count']
snps = pd.read_csv('data/emmax_out/EMMAX.0_5_FT10.top', skiprows=2, names=colnames)
geno = pd.read_csv('data/filtered.coded_call_method_54.tair9.FT10.csv')
pheno = pd.read_csv('data/filtered.FT10.txt', index_col=0)

# Sort snps by p-value and keep only top N snps
snps = snps.sort_values(by='pvalue', axis=0).head(N)

# Get snps from geno
snps = pd.merge(geno, snps, how='inner', on=['Chromosome', 'Positions'])
snps = snps[geno.columns]

# Separate into training and test
geno = snps.T.iloc[2:]
geno = geno.astype(float)
geno.index = geno.index.astype('int64')
X_train, X_test, y_train, y_test = train_test_split(geno, pheno, test_size=0.2, random_state=42)

# Linear regression model
model = LinearRegression()
model.fit(X_train, y_train)

y_pred = model.predict(X_test)

mse = mean_squared_error(y_test, y_pred)
print(f"MSE using linear regression: {mse}")
r2 = r2_score(y_test, y_pred)
print(f"R^2 using linear regression: {r2}")

plot(y_test,y_pred)

# Lasso (L1) regularization
model = Lasso(alpha=1)
model.fit(X_train, y_train)

y_pred = model.predict(X_test)

mse = mean_squared_error(y_test, y_pred)
print(f"MSE using lasso regularization: {mse}")
r2 = r2_score(y_test, y_pred)
print(f"R^2 using lasso regularization: {r2}")

plot(y_test,y_pred)

# Ridge (L2) regularization
model = Ridge(alpha=1)
model.fit(X_train, y_train)

y_pred = model.predict(X_test)

mse = mean_squared_error(y_test, y_pred)
print(f"MSE using ridge regularization: {mse}")
r2 = r2_score(y_test, y_pred)
print(f"R^2 using ridge regularization: {r2}")

plot(y_test,y_pred)




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