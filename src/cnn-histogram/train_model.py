import numpy as np
import pandas as pd
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import RobustScaler
import argparse
from keras.models import Sequential
from keras.layers import Dense, Dropout, Activation, Flatten
from keras.layers import Conv2D, MaxPooling2D
from keras.optimizers import SGD
from keras.callbacks import EarlyStopping, ModelCheckpoint, TensorBoard
from keras.utils import to_categorical
from pathlib import Path
from datetime import datetime

parser = argparse.ArgumentParser()

parser.add_argument('-save_dir',
                    dest='save_dir',
                    help='Path to save the results')

args = parser.parse_args()

def load_data(idx_samples, data_path, gene_conn_list):
    n_samples = 4
    tot_xdata = np.zeros((int((len(idx_samples)) * len(gene_conn_list)/n_samples), 32, 32, 1))
    tot_ydata = pd.DataFrame(columns=['interaction'])
    for i in idx_samples:
        file_name = f'Nxdata_tf_{str(i)}.npy'
        path = data_path / 'NEPDF_data' / file_name
        xdata = np.load(path)
        idx_low = int(i * len(gene_conn_list)/n_samples)
        idx_high = int((i+1) * len(gene_conn_list)/n_samples)
        ydata = gene_conn_list.iloc[idx_low:idx_high]['interaction'].to_frame()
        tot_xdata[idx_low:idx_high, :, :, :] = xdata
        tot_ydata = tot_ydata.append(ydata, ignore_index=True)
    print('Number of NaNs on y: %d' %tot_ydata['interaction'].isnull().sum())
    return tot_xdata, tot_ydata

def create_NN(X_train, num_classes):
    model = Sequential()
    model.add(Conv2D(32, (3, 3), padding='same',
                     input_shape=X_train.shape[1:]))
    model.add(Activation('relu'))
    model.add(Conv2D(32, (3, 3)))
    model.add(Activation('relu'))
    model.add(MaxPooling2D(pool_size=(2, 2)))
    model.add(Dropout(0.25))

    model.add(Conv2D(64, (3, 3), padding='same'))
    model.add(Activation('relu'))
    model.add(Conv2D(64, (3, 3)))
    model.add(Activation('relu'))
    model.add(MaxPooling2D(pool_size=(2, 2)))
    model.add(Dropout(0.25))

    model.add(Conv2D(128, (3, 3), padding='same'))
    model.add(Activation('relu'))
    model.add(Conv2D(128, (3, 3)))
    model.add(Activation('relu'))
    model.add(MaxPooling2D(pool_size=(2, 2)))
    model.add(Dropout(0.25))

    model.add(Flatten())
    model.add(Dense(512))
    model.add(Activation('relu'))
    model.add(Dropout(0.5))
    
    model.add(Dense(num_classes))
    model.add(Activation('softmax'))
    sgd = SGD(lr=0.01, decay=1e-6, momentum=0.9, nesterov=True)
    model.compile(optimizer=sgd,loss='categorical_crossentropy',metrics=['accuracy'])
    return model
    
            
if __name__ == '__main__':
    
    save_dir = Path(args.save_dir)
    # NN parameters
    batch_size = 1024 # mini batch for training
    num_classes = 3   # categories of labels (0: no iterataction; 1: A regulates B; 2: B regulates A)
    epochs = 100     # iterations of trainning, with GPU 1080, 600 for KEGG and Reactome, 200 for tasks for GTRD
    

    # Load the tresholded informration
    gene_conn_list = pd.read_csv('/output/TransformedData_DerivedData/CNN/histogram/thr_networks/corr_coeff_thr_CREBBP.csv')
    index_samples = [0, 1, 2, 3]
    xdata, ydata = load_data(index_samples, save_dir, gene_conn_list)

    # Split train test
    indices = ydata.index.to_list()
    X_train, X_test, y_train, y_test, idx_train, idx_test  = train_test_split(xdata, ydata, indices, test_size=0.33,
            random_state=42)
    print(X_train.shape, 'x_train samples')
    print(X_test.shape, 'x_test samples')
    
    if num_classes > 2:
        y_train = to_categorical(y_train, num_classes)
        y_test = to_categorical(y_test, num_classes)
    
    print(y_train.shape, 'y_train samples')
    print(y_test.shape, 'y_test samples')
    
    model = create_NN(X_train, num_classes)
    
    early_stopping = EarlyStopping(monitor='val_acc', patience=50, verbose=0, mode='auto')
    weights_path = save_dir / 'weights.hdf5'
    checkpoint1 = ModelCheckpoint(filepath=str(save_dir / 'weights-{epoch:02d}-{val_loss:.2f}.hdf5'), monitor='val_loss',
                                  verbose=1, save_best_only=False, save_weights_only=False, mode='auto', period=1)
    checkpoint2 = ModelCheckpoint(filepath=weights_path, monitor='val_acc', verbose=1,
                                  save_best_only=True, mode='auto', period=1)
    logdir = "logs/" + datetime.now().strftime("%Y%m%d-%H%M%S")
    tensorboard_callback = TensorBoard(log_dir=logdir)
    callbacks_list = [checkpoint2, early_stopping, tensorboard_callback]

    history = model.fit(X_train, y_train,
                  batch_size=batch_size,
                  epochs=epochs,validation_split=0.2,
                  shuffle=True, callbacks=callbacks_list)
    
    # Save model and weights
    model_path = save_dir / 'training'
    model_name = model_path / f'keras_cnn_trained_model_{epochs}.h5'
    model.save(model_name)
    print(f'Saved trained model at {model_path}')
    # Score trained model.
    scores = model.evaluate(X_test, y_test, verbose=1)
    print('Test loss:', scores[0])
    print('Test accuracy:', scores[1])
    y_predict = model.predict(X_test)
    
    # Save files
    indexes = {'idx_train': idx_train, 'idx_test': idx_test}
    np.save(model_path / 'train_test_idx.npy', indexes)
    np.save(model_path / 'end_y_test.npy', y_test)
    np.save(model_path / 'end_y_predict.npy', y_predict)



