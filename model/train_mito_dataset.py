import tensorflow as tf
from tensorflow.keras.callbacks import EarlyStopping
from tensorflow.keras.layers import Dense, Dropout, Flatten, LSTM, Conv1D, MaxPooling1D, Attention, Concatenate, Input, Layer, Softmax, Multiply
from tensorflow.keras.optimizers import SGD 
from keras.metrics import categorical_accuracy
import numpy as np
from Bio import SeqIO
from tensorflow.keras.callbacks import EarlyStopping
early_stopping = EarlyStopping(monitor='val_loss', patience=2)
np.random.seed(1337)

train_path = './data/human_1-1500.fasta'
valid_path = './data/human_1501-2000.fasta'
bases = "ACGT"
base_idx = dict((c, i) for i, c in enumerate(bases))
idx_base = dict((i, c) for i, c in enumerate(bases))

maxlen = 16  # 16-th Markov
step = 1
batch_size = 64
epochs = 20
input_dim = len(bases)
epoch_gap = 500

def read_fasta(path):
    records = list(SeqIO.parse(path, "fasta"))
    text = ""
    for record in records:
        text += str(record.seq)
    return text

def read_data(path):
    text = read_fasta(path)
    for i in range(0, len(text) - maxlen, step):
        sentence = text[i: i + maxlen]
        next_char = text[i + maxlen]
        yield sentence, next_char

def onehotencode(sentences, next_bases):
    x = np.zeros((batch_size, maxlen, len(bases)), dtype=np.float32)
    y = np.zeros((batch_size, len(bases)), dtype=np.float32)
    for i, sentence in enumerate(sentences):
        for t, char in enumerate(sentence):
            x[i, t, base_idx[char]] = 1
        y[i, base_idx[next_bases[i]]] = 1
    return x, y

def get_batch(stream):
    sentences = []
    next_bases = []
    for sentence, next_char in stream:
        sentences.append(sentence)
        next_bases.append(next_char)
        if len(sentences) == batch_size:
            data_tuple = onehotencode(sentences,next_bases) 
            yield data_tuple
            sentences = []
            next_bases = []

def model_CNN_LSTM_ATT():
    print('Build a CNN_LSTM_ATT model')
    inputs = Input(shape=(maxlen, len(bases)))
    conv_layer = Conv1D(filters=1024, 
                        kernel_size=4,
                        trainable=True, 
                        padding='valid', 
                        activation='relu', 
                        strides=1)(inputs)
    maxpool_layer = MaxPooling1D(pool_size=3)(conv_layer)
    dropout1 = Dropout(0.3)(maxpool_layer)

    lstm_layer = LSTM(256, return_sequences=True)(dropout1)

    attention_logits = Dense(1, activation='tanh')(lstm_layer)
    attention_weights = Softmax(axis=1)(attention_logits)
    attention = Multiply()([lstm_layer, attention_weights])
   
    merged = Concatenate(axis=-1)([lstm_layer, attention])

    dropout2 = Dropout(0.3)(merged)
    flatten = Flatten()(dropout2)
    dense1 = Dense(1024, activation='relu')(flatten)

    output = Dense(len(bases), activation='softmax')(dense1)

    model = tf.keras.Model(inputs=inputs, outputs=output)

    optimizer = tf.keras.optimizers.RMSprop(learning_rate=0.001)
    model.compile(loss='categorical_crossentropy', optimizer=optimizer)

    return model

def test_on_validset(epoch):
    valid_loss_value = []
    valid_acc_value = []
    for i, batch in enumerate(get_batch(read_data(valid_path))):
        _input = batch[0]
        _labels = batch[1]
        x = model.test_on_batch(_input,_labels)
        if(i%epoch_gap==0):
            valid_loss_value.append(x)
            acc_tensor = categorical_accuracy(_labels, model.predict_on_batch(_input)) 
            valid_acc_tmp = tf.reduce_mean(acc_tensor).numpy()
            valid_acc_value.append(valid_acc_tmp)

    average_loss = np.mean(valid_loss_value)
    average_acc = np.mean(valid_acc_value)
    print(epoch, '\t', average_loss, '\t', average_acc)
    if average_loss < 0.05:
        print('loss is less than 0.05.')
        # model.stop_training = True

    return valid_loss_value, valid_acc_value


model = model_CNN_LSTM_ATT()
model.build(input_shape=(batch_size, maxlen, len(bases)))
model.summary()

train_loss = [] 
train_acc = []
valid_loss = []
valid_acc = []
for epoch in range(epochs):
    print("this is train epoch: ", epoch)
    for i, batch in enumerate(get_batch(read_data(train_path))):
        _input = batch[0]
        _labels = batch[1]
        x = model.train_on_batch(_input,_labels) # loss
        if(i%epoch_gap==0): 
            train_loss.append(x)
            acc_tensor = categorical_accuracy(_labels, model.predict_on_batch(_input))
            train_accuracy = tf.reduce_mean(acc_tensor).numpy()
            train_acc.append(train_accuracy)
            print(epoch, '\t', x, '\t', train_accuracy)


    valid_loss_value, valid_acc_value = test_on_validset(epoch)
    valid_loss.extend(valid_loss_value)
    valid_acc.extend(valid_acc_value)
    if model.stop_training:
        break

saved_model_path = "./model/CNN_LSTM_ATT_model_{}".format(maxlen)
tf.saved_model.save(model, saved_model_path)
print("Saved model successfully")





