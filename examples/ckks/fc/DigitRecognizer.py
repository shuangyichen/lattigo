''' DIGIT RECOGNIZER '''
import os
os.environ['TF_CPP_MIN_LOG_LEVEL'] = '2'
# Imports
# import cv2
from tensorflow import keras
import numpy as np

# Image Processing & Prediction
def preprocess_model():

    # Initialization of Variables
    frame_size = 240
    x, y, w, h = 0, 0, frame_size, frame_size

    # Keras Model Initialization
    model = keras.models.load_model('/home/ubuntu/test0822/lattigo/examples/ckks/fc/lr.h5')

    for idx,i in enumerate(model.layers):
        # weight = i.get_weights()
        # print(i.name)
        # name = str(idx)+".npy"
        # np.save(name , weight)
        if "conv" in i.name:
            print(i.name)
            weight = i.get_weights()
            name0 = str(idx)+"0.csv"
            name1 = str(idx)+"1.csv"
            # print(type(weight[0]))
            print(weight[0].shape)
            filter = weight[0].reshape((weight[0].shape[3],weight[0].shape[2],weight[0].shape[1],weight[0].shape[0]))
            print(filter.shape)
            np.savetxt(name0,filter.reshape(-1), delimiter=',')
            np.savetxt(name1, weight[1], delimiter=',')
            # /.save(name1 , weight[1])
        elif "dense" in i.name:
            print(i.name)
            weight = i.get_weights()
            name0 = str(idx)+"0.csv"
            name1 = str(idx)+"1.csv"
            np.savetxt(name0, weight[0].reshape(-1), delimiter=',')
            np.savetxt(name1, weight[1], delimiter=',')
            # np.save(name0 , weight[0])
            # np.save(name1 , weight[1])


    # weights = np.load("0.npy",allow_pickle=True)
    # print(weights)
    # weight = model.get_weights()
    # np.save('weight.npy' , weight)# , fmt='%s',delimiter=',')

    # Initialize the Camera
    # camera = cv2.VideoCapture(0)

    # while True:

    #     # Grab the image
    #     return_value, image = camera.read()

    #     # Dummy Variable
    #     original_image = image

    #     # Highlight the readable area
    #     cv2.rectangle(image, (x, y), (x + w, y + h), (255, 0, 0), 2)

    #     # Image Processing
    #     image = cv2.resize(image[y:y + h, x:x + w], (28, 28))   # Downsize to 28x28 pixels
    #     im_gray = cv2.cvtColor(image, cv2.COLOR_BGR2GRAY)   # Convert to GRAYSCALE
    #     (thresh, im_bw) = cv2.threshold(im_gray, 128, 255, cv2.THRESH_BINARY | cv2.THRESH_OTSU) # Sharpen the image

    #     # Image ready for prediction
    #     x_pred = im_bw.reshape(28, 28, 1)   # Reshape to 28x28x1 (size*size*color_channel)
    #     batch = np.array([x_pred])  # Create a batch (batch_size * size * size * color_channel)

    #     output = str(model.predict_classes(batch, verbose=0))

    #     # Add the ouput to the visible area of the camera
    #     cv2.putText(original_image, "Predicted Value is " + output, (10, 300), cv2.FONT_HERSHEY_SIMPLEX, 0.7, (255, 0, 0), 2)
    #     cv2.putText(original_image, "Press \'E\' to Quit", (10, 340), cv2.FONT_HERSHEY_DUPLEX, 0.7, (0, 255, 0), 1)
    #     cv2.imshow('image', original_image)

    #     # Press 'E' to quit
    #     if cv2.waitKey(1) & 0xFF == ord('e'):
    #         break

    # # Release the camera
    # camera.release()
    # cv2.destroyAllWindows()


def main():
    preprocess_model()

if __name__ == '__main__':
    main()
