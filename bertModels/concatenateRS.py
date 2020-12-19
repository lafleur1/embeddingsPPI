#trying difference of two representations and then single layer and neuron 
#bert embedding models 
#Bert embeddings are just (1,782) for a protein.
#using random sample dataset 

from math import ceil

#from __future__ import print_function
import keras
from keras.models import Sequential, Model, load_model
from keras.layers import Dense, Dropout, Activation, Flatten, Input, Lambda, Subtract
from keras.layers import Conv1D, MaxPooling1D, Conv2D, MaxPooling2D, LSTM, GRU, Masking, TimeDistributed, BatchNormalization
from keras.layers import Concatenate, Reshape, Bidirectional, Softmax, LocallyConnected1D
from keras.callbacks import ModelCheckpoint, EarlyStopping
from keras import regularizers
from keras import backend as K
import keras.losses
from keras.utils import plot_model
import tensorflow as tf

import pandas as pd

import math

import os
import pickle
import numpy as np

import scipy.sparse as sp
import scipy.io as spio

import matplotlib.pyplot as plt
import matplotlib.cm as cm


from contextlib import redirect_stdout

import pickle
from sklearn.utils import shuffle

#set up generators - must open the correct npy files and return them from the bert_folder 

def simple_concatenation_model(layer1Size, layer2Size):
    #Inputs
    res_1 = Input(shape=(768,))
    res_2 = Input(shape=(768,))

    #Outputs
    true_interacts = Input(shape=(1,))
    
    #layers 
    layer_dense_pair_1 = Dense(layer1Size, activation='relu')
    layer_dense_pair_2 = Dense(layer2Size, activation = 'relu')
    pred_layer = Dense(1, activation='sigmoid', kernel_initializer='zeros')
    
    #network
    dense_out_pair1 = layer_dense_pair_1(Concatenate(axis=-1)([res_1, res_2])) #concatenate two sequences 
    dense_out_pair2 = layer_dense_pair_2(dense_out_pair1)
    pred_interacts = pred_layer(dense_out_pair2)

    interaction_model = Model(
        inputs=[
            res_1,
            res_2
        ],
        outputs=pred_interacts
    )
    interaction_loss = Lambda(sigmoid_nll, output_shape = (1,))([true_interacts, pred_interacts])

    loss_model = Model(
        [
            res_1,
            res_2,
            true_interacts
        ],
        interaction_loss
    )
    return interaction_model, loss_model

def setupDataframesForGenerators(data):
    #break up by train/valid/test splits already outlined in DF 
    train = data[data['Dataset'] == 'Train'].copy().reset_index(drop = True)
    valid = data[data['Dataset'] == 'Valid'].copy().reset_index(drop = True)
    test = data[data['Dataset'] == 'Test'].copy().reset_index(drop = True)
    return train, valid, test

def GeneratorFromDF(dataframe, folder, batch_size = 32, shuffle_data = True):
    #keras sequence generators is being a pain- making own custom generators     
    while True:
        dataframe = dataframe.sample(frac=1).reset_index(drop = True) #shuffle the rows of the dataframe
        for offset in range(0, dataframe.shape[0], batch_size):
            dfSlice = dataframe.iloc[offset:offset+batch_size]
            seq1 = dfSlice['1'].to_list()
            seq2 = dfSlice['2'].to_list()
            labels = dfSlice.Interacts.to_list()
            x_1 = []
            x_2 = []
            x_3 = []
            Y = []
            for i in range(0, len(seq1)):
                #tupleToLoad = tuples[i].split(",")
                loadForward = np.load(folder + str(seq1[i]) + ".npy").flatten()
                loadedBackwards = np.load(folder + str(seq2[i]) + ".npy").flatten()
                x_1.append(loadForward)
                x_2.append(loadedBackwards)
                x_3.append(labels[i])
            x_1 = np.array(x_1)
            x_2 = np.array(x_2)
            x_3 = np.array(x_3)
            yield [x_1, x_2, x_3], x_3
            

def sigmoid_nll(inputs) :
    y_true, y_pred = inputs
    y_pred = K.clip(y_pred, K.epsilon(), 1.0 - K.epsilon())
    return K.mean(-y_true * K.log(y_pred) - (1.0 - y_true) * K.log(1.0 - y_pred), axis=-1)


def mean_squared_error(inputs) :
    y_true, y_pred = inputs
    return K.mean(K.square(y_pred - y_true), axis=-1)

def mean_absolute_error(inputs) :
    y_true, y_pred = inputs
    return K.mean(K.abs(y_pred - y_true), axis=-1)

def logcosh(inputs) :
    y_true, y_pred = inputs
    def _logcosh(x):
        return x + K.softplus(-2. * x) - K.log(2.)
    return K.mean(_logcosh(y_pred - y_true), axis=-1)


def simple_concatenation_model(layer1Size, layer2Size):
    #Inputs
    res_1 = Input(shape=(768,))
    res_2 = Input(shape=(768,))

    #Outputs
    true_interacts = Input(shape=(1,))
    
    #layers 
    layer_dense_pair_1 = Dense(layer1Size, activation='relu')
    layer_dense_pair_2 = Dense(layer2Size, activation = 'relu')
    pred_layer = Dense(1, activation='sigmoid', kernel_initializer='zeros')
    
    #network
    dense_out_pair1 = layer_dense_pair_1(Concatenate(axis=-1)([res_1, res_2])) #concatenate two sequences 
    dense_out_pair2 = layer_dense_pair_2(dense_out_pair1)
    pred_interacts = pred_layer(dense_out_pair2)

    interaction_model = Model(
        inputs=[
            res_1,
            res_2
        ],
        outputs=pred_interacts
    )
    interaction_loss = Lambda(sigmoid_nll, output_shape = (1,))([true_interacts, pred_interacts])

    loss_model = Model(
        [
            res_1,
            res_2,
            true_interacts
        ],
        interaction_loss
    )
    return interaction_model, loss_model


#set up keras callbacks, adam, etc. 
def trainEvalAndSaveModel(model_name, size1, size2,lr = 0.001, epochs = 1):
    opt = keras.optimizers.Adam(lr, beta_1=0.9, beta_2=0.999)
    callbacks =[EarlyStopping(monitor='val_loss', min_delta=0.001, patience=3, verbose=0, mode='auto')]
    iModel, loss_model = simple_concatenation_model(size1,size2)
    loss_model.compile(loss=lambda true, pred: pred, optimizer=opt, metrics = ['accuracy'])
    trainSteps = ceil(train.shape[0] / 32)
    validationSteps = ceil(valid.shape[0] / 32)
    print ("trianing steps: ", ceil(train.shape[0] / 32))
    print ("validation steps: ", ceil(valid.shape[0] / 32))

    history = loss_model.fit_generator(generator=trainGenerator, steps_per_epoch = trainSteps, 
                        validation_data=validationGenerator,validation_steps = validationSteps,
                        epochs=epochs,
                        use_multiprocessing=False,
                        workers=1, callbacks=callbacks)
    
    
    # summarize history for accuracy
    plt.plot(history.history['acc'])
    plt.plot(history.history['val_acc'])
    plt.title('model accuracy')
    plt.ylabel('accuracy')
    plt.xlabel('epoch')
    plt.legend(['train', 'validation'], loc='upper left')
    plt.show()
    # summarize history for loss
    plt.plot(history.history['loss'])
    plt.plot(history.history['val_loss'])
    plt.title('model loss')
    plt.ylabel('loss')
    plt.xlabel('epoch')
    plt.legend(['train', 'validation'], loc='upper left')
    plt.show()
    #save model to folder
    # Save model and weights
    
    save_dir = 'saved_models'
    
    if not os.path.isdir(save_dir):
        os.makedirs(save_dir)

    model_path = os.path.join(save_dir, model_name + '.h5')
    history_path =  os.path.join(save_dir, model_name + "_history.csv")

    iModel.save(model_path)
    print('Saved scrambler model at %s ' % (model_path))
    
    #save history to csv 
    historyDF = pd.DataFrame(history.history)
    historyDF.to_csv(history_path)


#set up datsets
totalSet = pd.read_csv("randomBreakup.csv")
train, valid, test = setupDataframesForGenerators(totalSet)
print ('train size: ', train.shape)
print ("valid size: ", valid.shape)
print ("test size: ", test.shape)
print (totalSet.columns) #Index(['Unnamed: 0', '1', '2', 'Interacts', 'Type', 'Dataset'], dtype='object')

#set up generators 
trainGenerator = GeneratorFromDF(train, "../bert_base_embeddings/")
validationGenerator = GeneratorFromDF(valid,  "../bert_base_embeddings/")

combos = [("concat_rs_300_100_lr_001", 300, 100, 0.001)]

for combo in combos:
    trainEvalAndSaveModel(combo[0], combo[1], combo[2],  combo[3])