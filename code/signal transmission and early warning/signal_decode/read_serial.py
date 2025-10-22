import serial
import numpy as np
from .alarm import check_alarm_level

def read_acceleration_data(ser):
    try:
        # Read a line of data from the serial port (assuming it's a comma-separated float sequence)
        line = ser.readline().decode('utf-8').strip()
        # Convert the string data into a list of floats
        data = list(map(float, line.split(',')))
        
        # Check the alarm level for the current data
        alarm_level, peak_value = check_alarm_level(data)
        print(f"Data Peak Value: {peak_value} | Alarm Level: {alarm_level}")
        
        return data, alarm_level
    except Exception as e:
        print(f"Error reading data from serial port: {e}")
        return None, None
