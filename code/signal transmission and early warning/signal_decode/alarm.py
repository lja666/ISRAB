import numpy as np

def check_alarm_level(data):
    """
    Check the peak value of the acceleration data and return an alarm level.
    """
    peak_value = np.max(np.abs(data))  # Get the peak value (absolute value)
    
    if peak_value > 300:
        return "Severe Alarm", peak_value
    elif peak_value > 100:
        return "High Alarm", peak_value
    elif peak_value > 50:
        return "Moderate Alarm", peak_value
    else:
        return "Normal", peak_value
