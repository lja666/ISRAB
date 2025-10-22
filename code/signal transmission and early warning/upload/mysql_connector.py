import pymysql
from datetime import datetime

def insert_to_mysql(data):
    try:
        # Connect to MySQL database
        connection = pymysql.connect(
            host='localhost',      # Database host address
            user='user',  # Username
            password='******',  # Password
            database='shm.qd_dataset_2625'  # Database name
        )

        with connection.cursor() as cursor:
            timestamp = datetime.now()  # Get the current timestamp
            # Convert acceleration data into a string (comma-separated)
            acceleration_values = ','.join(map(str, data))

            # SQL query to insert the data
            sql = "INSERT INTO acceleration_data (timestamp, acceleration_values) VALUES (%s, %s)"
            cursor.execute(sql, (timestamp, acceleration_values))
        
        connection.commit()
    except Exception as e:
        print(f"Error inserting data into MySQL: {e}")
    finally:
        connection.close()
