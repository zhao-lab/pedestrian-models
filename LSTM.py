
import numpy as np
import matplotlib.pyplot as plt
from pandas import read_csv
import math
from keras.models import Sequential
from keras.layers import Dense,TimeDistributed,Masking, Dropout
from keras.layers import LSTM
from sklearn.preprocessing import MinMaxScaler
from sklearn.metrics import mean_squared_error
import scipy.io
mat = scipy.io.loadmat('discrete_1134.mat')
batch_size=1134
train_size=int(round(batch_size*0.7))
valid_size=int(round(batch_size*0.2))
test_size=batch_size-train_size-valid_size
max_time=150
prediction_step=1
dis_range=8    # 0~8m/s
dis_binsize=0.25    #every 0.25 per bin 
N_bin=round(dis_range/dis_binsize)
k=1
U_org=mat['U']
v_org=mat['v']
U=np.zeros((batch_size,max_time,3))
v_con=np.zeros((batch_size,max_time,1))   #continuous velocity
v=np.zeros((batch_size,max_time,N_bin))
print(type(U_org[0,k-1]))

for k in range(0,batch_size):
    a=U_org[0,k]
    U[k,0:a.shape[0],:]=a[:,0:3]
    U[k,a.shape[0]-prediction_step:,:]=0.     #make sure to mask out the U such that the non-zero U and v are equal in length
   # U[k-1,t_s:61,:]=np.zeros((61-t_s,4))
    a=v_org[0,k-1]
    v_con[k,0:a.shape[0],:]=a
    #v_con[k,a.shape[0]:,:]=-np.ones((max_time-a.shape[0],1))
    
for k in range(0,batch_size):
    a=U_org[0,k]
    for j in range(0,a.shape[0]):
        bin_index=round(v_con[k,j,0]/dis_binsize)
        if bin_index>N_bin-1:
            bin_index=N_bin-1
        
        if bin_index<0:
            bin_index=0
        
        v[k,j,int(bin_index)]=1
        

U_train=U[0:train_size,0:max_time-prediction_step,0:3]
v_train=v[0:train_size,prediction_step:max_time,:]
U_validation=U[train_size:train_size+valid_size,0:max_time-prediction_step,0:3]
v_validation=v[train_size:train_size+valid_size,prediction_step:max_time,:]
U_test=U[train_size+valid_size:,0:max_time-prediction_step,0:3]

np.random.seed(7)

model = Sequential()
model.add(Masking(mask_value=0.,batch_input_shape=(train_size,max_time-prediction_step, 3)))
model.add(LSTM(
        batch_input_shape=(train_size,max_time-prediction_step,3),   #236 is the batch size(number of independent samples, 62 is the maximum time length, 4 is the input data dimension)
        output_dim=25,   #hidden layer unit number  25
        return_sequences=True,
        stateful=False,
        activation='relu'
        ))
model.add(Dropout(0.2))
model.add(LSTM(
            35,    #35
            return_sequences=True,
            activation='relu'))
model.add(Dropout(0.2))
#model.add(LSTM(10,return_sequences=True))   #the second RNN layer
model.add(TimeDistributed(Dense(N_bin, activation='softmax')))   # is the output dim
model.compile(loss='categorical_crossentropy', optimizer='adam',metrics=['accuracy'])
history = model.fit(U_train, v_train, epochs=200, verbose=2, validation_data=(U_validation, v_validation), shuffle=False)


#U_test=U[1018,0:max_time-prediction_step,0:3].reshape([1,max_time-prediction_step,3])
#v_test=v[1018,prediction_step:max_time,:].reshape([1,max_time-prediction_step,N_bin])
testPredict=model.predict(U_test)
v_test=v_con[train_size+valid_size:,0:max_time-prediction_step,:]
scipy.io.savemat('LSTM_mask_prediction',{'predicted_v':testPredict,'true_v':v_test})
#==============================================================================
# v_con_test=np.zeros((max_time-prediction_step))
# predict_con=np.zeros((max_time-prediction_step))
# for k in range(0,max_time-prediction_step):
#     v_con_test[k]=np.argmax(v_test[0,k,:])*dis_binsize
#     predict_con[k]=np.argmax(testPredict[0,k,:])*dis_binsize
#     
# plt.plot(predict_con)
# plt.plot(v_con_test)
#==============================================================================
#plt.show()

