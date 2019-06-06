import numpy as np 
import os
import numpy as np
from keras.models import *
from keras.layers import *
from keras.optimizers import *
from keras.callbacks import ModelCheckpoint, LearningRateScheduler
from keras import backend as keras
from keras.activations import softmax

num_classes=2

def shakenet(pretrained_weights = None,input_size = (6002,1)):
    inputs = Input(input_size)
    conv1=Conv1D(8, 10, activation = 'relu', padding = 'valid', kernel_initializer = 'he_normal')(inputs)
    conv1=Conv1D(8, 10, activation = 'relu', padding = 'valid', kernel_initializer = 'he_normal')(conv1)
    pool1=MaxPooling1D(pool_size=16)(conv1)
    conv2=Conv1D(8, 10, activation = 'relu', padding = 'valid', kernel_initializer = 'he_normal')(pool1)
    conv2=Conv1D(8, 10, activation = 'relu', padding = 'valid', kernel_initializer = 'he_normal')(conv2)
    pool2=MaxPooling1D(pool_size=16)(conv2)
    lstm3=LSTM(10)(pool2)
    dense3=Dense(num_classes,activation='softmax')(lstm3)

    model = Model(input = inputs, output = dense3)

    model.compile(optimizer = Adam(lr = 1e-3), loss = 'categorical_crossentropy', metrics = ['accuracy'])

    if(pretrained_weights):
    	model.load_weights(pretrained_weights)

    return model


