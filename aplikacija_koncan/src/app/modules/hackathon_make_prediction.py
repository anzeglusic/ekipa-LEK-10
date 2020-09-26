from sklearn.neighbors import KNeighborsRegressor
from sklearn.model_selection import train_test_split
from sklearn import preprocessing

import pandas as pd

class ML_model:
	def __init__(self, path, n_neighbors=20, weights="distance"):
		self.model = self.create_model(path, n_neighbors=n_neighbors, weights=weights)

	def create_model(self, path, n_neighbors, weights):
		'''
		Function that creates the model and returns it.
		'''
		model = KNeighborsRegressor(n_neighbors=n_neighbors, weights=weights) # creates the model
		# trains / fits the data
		df = pd.read_csv(path) # importing our dataset for training, last column is the BA->target to predict
		#LEON
		df.drop(columns= ['Unnamed: 0', 'Name', 'Biološka uporabnost (%)', 'Updated SMILES'], inplace=True)
		

		X_before_norm = df.iloc[:, :-1] 
		#X = self.normalize(X_before_norm) # Normalize the features
		X = X_before_norm
		y = df.iloc[:, -1]
		
		self.X_train, self.X_test, self.y_train, self.y_test = train_test_split(X, y, test_size=0.2) # make train/test split, it is randomized
		model.fit(self.X_train, self.y_train) # the actual fitting of the model
		
		return model

	def normalize(self, features_dataset):
		'''
		For better model performance, all the features must be individually scaled (between -1 and 1) or
		normalized (between 0 and 1).
		'''
		min_max_scaler = preprocessing.MinMaxScaler()
		X = min_max_scaler.fit_transform(features_dataset)
		X = pd.DataFrame(X, columns=features_dataset.columns)
		return X

	def make_prediction(self, features):
		BA = self.model.predict([features])[0]
		return BA

if __name__ == "__main__":
	df = pd.read_csv("clean_10.csv") # importing our dataset for training, last column is the BA->target to predict
	X = df.iloc[:, :-1]
	model = ML_model("clean_10.csv")
	BA = model.make_prediction(X.iloc[0, :])
	print(f"New drug BA: {BA:2.2f}%")