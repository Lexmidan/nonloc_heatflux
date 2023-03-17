# ML-class
## TensorFlow 9/21/22

NN template
1. define the model
  # Create sequential NN object
  `model = tf.keras.Sequential()`

  # Add input layer: size coresponds to the number of x components
  `model.add(tf.keras.layers.InputLayer(input_shape=(features.shape[1],), name='Input_layer'))`

  # Add hidden layers
  `model.add(tf.keras.layers.Dense(64, activation='relu', name='Dense_layer_1'))`
  `model.add(tf.keras.layers.Dense(32, name='Dense_layer_2'))`
  `model.add(tf.keras.layers.Dense(16, name='Dense_layer_3'))`
  `model.add(tf.keras.layers.Dense(8, name='Dense_layer_4'))`

  # Add output layer: size corresponds to the number of y components
  `model.add(tf.keras.layers.Dense(target.shape[1], name='Output_layer'))`

  In case of needing more customized NN try to use (see `APIs.ipynb`)
  `https://www.tensorflow.org/guide/keras/functional`

2. compile the model
  `model.compile(tf.optimizers.RMSprop(0.0001), loss='binary_crossentropy', metrics=['accuracy'])`
  `model.compile(tf.optimizers.RMSprop(0.001), loss='mse')`

3. fit the model (learning from x->y)
  `model.fit(x=features.to_numpy(), y=target.to_numpy(),epochs=100, callbacks=[tensorboard_callback] , batch_size=32, validation_split=0.2)`

4. Evaluate the model (should be data not used in the learning process)
  `loss, accuracy = model.evaluate(features.to_numpy(), target.to_numpy())`


The three most common loss functions are:

`binary_crossentropy` for binary classification

`sparse_categorical_crossentropy` for multi-class classification

`mse` (mean squared error) for regression
