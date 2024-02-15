import numpy as np
import pandas as pd
import torch

from sklearn.model_selection import train_test_split
from torch.utils.data import Dataset, DataLoader

# Define the linear model
class LinearRegression(torch.nn.Module):
    def __init__(self, input_size):
        super(LinearRegression, self).__init__()
        self.linear = torch.nn.Linear(input_size, 1)

    def forward(self, x):
        return self.linear(x)


# Define the Dataset class
class GenoPhenoDataset(Dataset):
    def __init__(self, genotypes, phenotypes):
        self.genotypes = genotypes
        self.phenotypes = phenotypes

    def __len__(self):
        return len(self.phenotypes)

    def __getitem__(self, idx):
        return self.genotypes[idx], self.phenotypes[idx]


# Read in data and trim to ssnps
geno = pd.read_csv('data/coded_call_method_54.tair9.FT10.csv')
geno.columns = [int(col) if col.isdigit() else col for col in geno.columns] # convert ids to ints
pheno = pd.read_csv('FT10.txt', sep='\t')
ssnps = pd.read_csv('data/ssnps.csv', index_col=0)

# Impute missing phenotypes
mean = pheno['5_FT10'].mean()
pheno.fillna(mean, inplace=True)

# Filter by ssnps (from map.py)
geno_ssnps = geno.loc[ssnps.index]

# Filter and sort pheno to match ids in geno
geno_ids = geno_ssnps.iloc[:,2:].columns.tolist()
pheno = pheno[pheno.iloc[:,0].isin(geno_ids)]
pheno['sort_order'] = pheno.iloc[:, 0].apply(lambda x: geno_ids.index(x))
pheno = pheno.sort_values('sort_order')
pheno.drop('sort_order', axis=1, inplace=True)
pheno = pheno.reset_index(drop=True)

# Separate into training and test
geno_T = geno_ssnps.T.iloc[2:]
X_train, X_test, y_train, y_test = train_test_split(geno_T, pheno['5_FT10'], test_size=0.2, random_state=42)

# Convert sets to tensors
X_train_tensor = torch.tensor(X_train.values, dtype=torch.float32)
y_train_tensor = torch.tensor(y_train.values, dtype=torch.float32)
X_test_tensor = torch.tensor(X_test.values, dtype=torch.float32)
y_test_tensor = torch.tensor(y_test.values, dtype=torch.float32)

# Parameters
input_size = X_train_tensor.shape[1]
learning_rate = 0.01
batch_size = 1
epochs = 100

# Dataset and DataLoader
dataset = GenoPhenoDataset(X_train_tensor, y_train_tensor)
dataloader = DataLoader(dataset, batch_size=batch_size, shuffle=True)

# Model, Loss, Optimizer
model = LinearRegression(input_size)
criterion = torch.nn.MSELoss()
optimizer = torch.optim.Adam(model.parameters(), lr=learning_rate)

# Training Loop

for epoch in range(epochs):
    for i, (inputs, targets) in enumerate(dataloader):
        # Forward pass
        outputs = model(inputs)
        loss = criterion(outputs, targets)

        # Backward pass and optimization
        optimizer.zero_grad()
        loss.backward()
        torch.nn.utils.clip_grad_norm_(model.parameters(), max_norm=1.0)
        optimizer.step()

    # Print loss after every epoch
    print(f'Epoch [{epoch + 1}/{epochs}], Loss: {loss.item():.4f}')


# Make predictions
model.eval()
with torch.no_grad():
    predictions = model(X_test_tensor)

test_loss = criterion(predictions, y_test_tensor.unsqueeze(1))
print(f'Test Mean Squared Error: {test_loss.item():.4f}')