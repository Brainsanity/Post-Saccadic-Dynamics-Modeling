# fit a final xgboost model on the housing dataset and make a prediction
import numpy as np
from pandas import read_csv
import xgboost as xgb
from sklearn.multioutput import MultiOutputRegressor
import matplotlib.pyplot as plt

# load the dataset
url = 'C:/0 Codes/SmoothPursuit/Matlab/saccades.txt' #'https://raw.githubusercontent.com/jbrownlee/Datasets/master/housing.csv'
outFile = 'C:/0 Codes/SmoothPursuit/Matlab/prediction.txt'
dataframe = read_csv(url, header=None)
data = dataframe.values
# split dataset into input and output columns
nFeatures = np.int(data.shape[1]/2)
x, y = data[:, :nFeatures], data[:, nFeatures:]
print(x.shape)
print(y.shape)

N = np.int(np.round(x.shape[0] * 0.9))
trainx = x[0:N,:]
trainy = y[0:N,:]
testx = x[N:,:]
testy = y[N:,:]
# # define model
# model = xgb.XGBRegressor(n_estimators=1000, max_depth=7, eta=0.1, subsample=0.7, colsample_bytree=0.8)
# # fit model
# model.fit(x, y)
# # define new data
# row = [0.00632,18.00,2.310,0,0.5380,6.5750,65.20,4.0900,1,296.0,15.30,396.90,4.98]
# new_data = np.asarray([row])
# # make a prediction
# yhat = model.predict(x[1,:]) #(new_data)

# # summarize prediction
# print('Ground true: %.1f', y[1,:])
# print('Predicted: %.1f' % yhat)


## Use sklearn.multioutput.MultiOutputRegressor as a wrapper of xgb.XGBRegressor
# fitting
reg = MultiOutputRegressor(xgb.XGBRegressor(objective='reg:pseudohubererror', eval_metric='rmse', verbose=True)).fit(trainx, trainy)
trainyHat = reg.predict(trainx)
testyHat = reg.predict(testx)

np.savetxt(outFile, reg.predict(x), fmt=' %f', delimiter=',')

# predicting
print(np.mean((trainyHat - trainy).flatten()**2, axis=0))
print(np.mean((testyHat - testy).flatten()**2, axis=0))

plt.subplot(2,2,1)
for i in np.arange(0, trainx.shape[0]):
	# plt.plot(np.arange(0,100), trainx[i,:], 'b')
	plt.plot(np.arange(0,100), trainy[i,:], 'r')
	plt.plot(np.arange(0,100), trainyHat[i,:], 'g')

plt.subplot(2,2,3)
for i in np.arange(0, trainx.shape[0]):
	# plt.plot(np.arange(0,100), trainx[i,:], 'b')
	# plt.plot(np.arange(0,100), trainy[i,:], 'r')
	plt.plot(np.arange(0,100), trainyHat[i,:] - trainy[i,:], 'k')


plt.subplot(2,2,2)
for i in np.arange(0, testx.shape[0]):
	# plt.plot(np.arange(0,100), testx[i,:], 'b')
	plt.plot(np.arange(0,100), testy[i,:], 'r')
	plt.plot(np.arange(0,100), testyHat[i,:], 'g')

plt.subplot(2,2,4)
for i in np.arange(0, testx.shape[0]):
	# plt.plot(np.arange(0,100), testx[i,:], 'b')
	# plt.plot(np.arange(0,100), testy[i,:], 'r')
	plt.plot(np.arange(0,100), testyHat[i,:] - testy[i,:], 'k')

plt.show()