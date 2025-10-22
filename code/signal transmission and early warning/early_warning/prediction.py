import tensorflow as tf
import numpy as np
from tensorflow.keras.models import load_model
from .alarm import check_alarm_level

def load_prediction_model(model_path):
    try:
        model = load_model(model_path)  # Load the model
        print("Model loaded successfully.")
        return model
    except Exception as e:
        print(f"Error loading model: {e}")
        return None

def predict_acceleration_data(model, data):
    try:
        # Reshape data to match the model's input shape
        data = np.reshape(data, (1, -1))  # Reshape to (1, n_features)
        prediction = model.predict(data)
        
        # Check the prediction's alarm level
        alarm_level, peak_value = check_alarm_level(prediction)
        print(f"Prediction Peak Value: {peak_value} | Alarm Level: {alarm_level}")
        
        return prediction, alarm_level
    except Exception as e:
        print(f"Error during prediction: {e}")
        return None, None
