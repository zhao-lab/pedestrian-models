
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
max_time=150
prediction_step=1
dis_range=8    # 0~8m/s
dis_binsize=0.25    #every 0.25 per bin 
past_input=3    #use past three step input
N_bin=round(dis_range/dis_binsize)
k=1
U_org=mat['U']
v_org=mat['v']
#U=np.zeros((batch_size,max_time,3))
#v_con=np.zeros((batch_size,max_time,1))   #continuous velocity
#v=np.zeros((batch_size,max_time,N_bin))
print(type(U_org[0,k-1]))

#get training data
U=np.zeros((1,past_input*3))
v_con=np.zeros((1,1))
for k in range(0,batch_size):
    a=U_org[0,k]
    b=a.shape[0]
    for j in range(0,b-past_input-prediction_step):
        U=np.concatenate((U, np.reshape(a[j:j+past_input,0:3],(1,past_input*3))))
        
    a=v_org[0,k]
    for j in range(0,b-past_input-prediction_step):
        v_con=np.concatenate((v_con, np.reshape(a[j+past_input+prediction_step,0],(1,1))))
        
U=U[1:,:]   #discard the first element that is zero 
v_con=v_con[1:,:]    
v=np.zeros((v_con.shape[0],N_bin))
    
for k in range(0,v_con.shape[0]):
    bin_index=round(v_con[k,0]/dis_binsize)
    if bin_index>N_bin-1:
        bin_index=N_bin-1
        
    if bin_index<0:
        bin_index=0
        
    v[k,int(bin_index)]=1
        
total_size=v_con.shape[0]
train_size=int(round(total_size*0.7))
valid_size=int(round(total_size*0.2))
test_size=total_size-train_size-valid_size

U_train=U[0:train_size,:]
v_train=v[0:train_size,:]
U_validation=U[train_size:train_size+valid_size,:]
v_validation=v[train_size:train_size+valid_size,:]
U_test=U[train_size+valid_size:,:]
v_test=v_con[train_size+valid_size:,:]

scipy.io.savemat('rawUv_1steps',{'U_train':U_train,'v_train':v_train,'U_validation':U_validation,'v_validation':v_validation,'U_test':U_test,'v_test':v_test})


#==============================================================================
# np.random.seed(7)
# 
# 
# model = Sequential()
# #model.add(Masking(mask_value=0,batch_input_shape=(train_size, max_time-prediction_step, 3)))
# model.add(Dense(
#         batch_input_shape=(train_size,3*past_input),   #236 is the batch size(number of independent samples, 62 is the maximum time length, 4 is the input data dimension)
#         output_dim=175,   #hidden layer unit number
#         activation='relu'
#         ))
# model.add(Dropout(0.2))
# model.add(Dense(
#         output_dim=175,   #hidden layer unit number
#         activation='relu'
#         ))
# model.add(Dropout(0.2))
# #model.add(LSTM(10,return_sequences=True))   #the second RNN layer
# model.add(Dense(N_bin, activation='softmax'))   # is the output dim
# model.compile(loss='categorical_crossentropy', optimizer='adam',metrics=['accuracy'])
# history = model.fit(U_train, v_train, epochs=500, verbose=2, validation_data=(U_validation, v_validation), shuffle=False)
# #print(history.history['loss'][-1])
# #loss_record[i,j]=history.history['loss'][-1]
# 
# #U_test=U[1018,0:max_time-prediction_step,0:3].reshape([1,max_time-prediction_step,3])
# #v_test=v[1018,prediction_step:max_time,:].reshape([1,max_time-prediction_step,N_bin])
# 
# testPredict=model.predict(U_test)
# scipy.io.savemat('FF_test_prediction_1step',{'predicted_v':testPredict,'true_v':v_test})
#==============================================================================
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

