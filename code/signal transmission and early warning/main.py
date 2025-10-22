


import time
from signal_decode.read_serial import read_acceleration_data
from upload.mysql_connector import insert_to_mysql
from early_warning.prediction import load_prediction_model, predict_acceleration_data
from config.config import COM_PORT, BAUDRATE

def main():
    # Set up serial communication
    ser = serial.Serial(COM_PORT, baudrate=BAUDRATE, timeout=1)

    # Load the pre-trained model
    model_path = 'models/Saved DEC model 1028.hdf5'  # Replace with your model path
    model = load_prediction_model(model_path)

    if model is None:
        print("Model loading failed. Exiting.")
        return

    print("Starting to read data from serial port...")

    while True:
        # Continuously read acceleration data from the serial port
        acceleration_data, data_alarm_level = read_acceleration_data(ser)

        if acceleration_data is not None:
            # Insert the data into MySQL
            insert_to_mysql(acceleration_data)

            # Make a prediction using the model
            prediction, prediction_alarm_level = predict_acceleration_data(model, acceleration_data)

            # Display alarm levels for data and prediction
            print(f"Data Alarm Level: {data_alarm_level}")
            print(f"Prediction Alarm Level: {prediction_alarm_level}")

        # Control the frequency of data reading
        time.sleep(1)  # Read data every second

if __name__ == "__main__":
    main()
