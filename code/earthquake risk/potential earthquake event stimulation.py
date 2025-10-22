# -*- coding: utf-8 -*-
"""
Created on Tue Oct 21 23:35:35 2025

@author: Administrator
"""

import geopandas as gpd
import numpy as np
import pandas as pd
import requests
import json
from shapely.geometry import Point, LineString, box
from geopy.distance import geodesic
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import random
from datetime import datetime, timedelta
import warnings
import os
warnings.filterwarnings('ignore')

#%% 
print("=== earthquake risk analysis ===")
print("run ...\n")

# ========== Initialize Variables ========== 
faults = None
fault_data = []
gr_parameters = {}
simulation_results = []
earthquake_data = None

# ========== Set Parameters ========== 
# Example parameters (should be modified based on actual situation)
shapefile_path = "\shapefile\gem_active_faults.shp"  # Path to your shapefile
bbox_coords = [126.0, 31.0, 136.0, 44.0]  # Example bounding box [minx, miny, maxx, maxy]
# 44.211231°N 147.007985°E
print(f"Analysis area: {bbox_coords}")
print(f"Shapefile path: {shapefile_path}")
print("Note: If the shapefile is missing, the program will use simulated data\n")

# ========== Step 1: Extract Faults from Shapefile within the Given Range ========== 
print("Step 1: Extracting faults from shapefile within the specified range...")

try:
    # Read shapefile
    gdf = gpd.read_file(shapefile_path)
    gdf.set_crs("EPSG:4326", allow_override=True, inplace=True)

    # Ensure the coordinate system is WGS84
    if gdf.crs != 'EPSG:4326':
        gdf = gdf.to_crs('EPSG:4326')
    
    # Create bounding box
    bbox = box(bbox_coords[0], bbox_coords[1], bbox_coords[2], bbox_coords[3])
    
    # Extract features within the bounding box
    faults = gdf[gdf.geometry.intersects(bbox)].copy()
    faults = faults.reset_index(drop=True)
    
    print(f"Found {len(faults)} faults within the specified range")
    
except FileNotFoundError:
    print("Shapefile not found, creating simulated fault data...")
    
    # Create simulated fault data
    np.random.seed(42)
    n_faults = 5
    minx, miny, maxx, maxy = bbox_coords
    
    for i in range(n_faults):
        # Generate random fault endpoints
        start_lon = np.random.uniform(minx, maxx)
        start_lat = np.random.uniform(miny, maxy)
        
        # Generate fault direction and length
        angle = np.random.uniform(0, 2*np.pi)
        length_deg = np.random.uniform(0.05, 0.3)  # Length in degrees
        
        end_lon = start_lon + length_deg * np.cos(angle)
        end_lat = start_lat + length_deg * np.sin(angle)
        
        # Ensure endpoint is within the bounding box
        end_lon = np.clip(end_lon, minx, maxx)
        end_lat = np.clip(end_lat, miny, maxy)
        
        # Create LineString geometry
        geom = LineString([(start_lon, start_lat), (end_lon, end_lat)])
        
        fault_info = {
            'fault_id': f"demo_fault_{i}",
            'start_lon': start_lon,
            'start_lat': start_lat,
            'end_lon': end_lon,
            'end_lat': end_lat,
            'geometry': geom
        }
        fault_data.append(fault_info)
    
    print(f"Created {len(fault_data)} simulated faults")

#%% 
# Assuming gdf is the GeoDataFrame that has been read and CRS set
fig, ax = plt.subplots(figsize=(10, 10), dpi=300)  # High resolution DPI

# Plot map
gdf.plot(ax=ax, color='lightblue', edgecolor='black', linewidth=2)

# Set title
plt.title("All faults", fontsize=16)    

plt.show()

#%% 
fig, ax = plt.subplots(figsize=(10, 10), dpi=300)  # High resolution DPI

# Plot map
faults.plot(ax=ax, color='lightblue', edgecolor='black', linewidth=2)

# Set title
plt.title("Selected faults", fontsize=16)   

#%% 
# ========== Step 2: Calculate Fault Endpoint Coordinates ========== 
print("\nStep 2: Calculating fault endpoint coordinates...")

if faults is not None:
    for idx, fault in faults.iterrows():
        geom = fault.geometry
        
        if geom.geom_type == 'LineString':
            coords = list(geom.coords)
            start_point = coords[0]  # (lon, lat)
            end_point = coords[-1]   # (lon, lat)
            
            fault_info = {
                'fault_id': idx,
                'start_lon': start_point[0],
                'start_lat': start_point[1],
                'end_lon': end_point[0],
                'end_lat': end_point[1],
                'geometry': geom
            }
            fault_data.append(fault_info)
        
        elif geom.geom_type == 'MultiLineString':
            # Handle MultiLineString
            for i, line in enumerate(geom.geoms):
                coords = list(line.coords)
                start_point = coords[0]
                end_point = coords[-1]
                
                fault_info = {
                    'fault_id': f"{idx}_{i}",
                    'start_lon': start_point[0],
                    'start_lat': start_point[1],
                    'end_lon': end_point[0],
                    'end_lat': end_point[1],
                    'geometry': line
                }
                fault_data.append(fault_info)

print(f"Calculated endpoint coordinates for {len(fault_data)} fault segments")

#%% 
# ========== Step 3: Calculate Fault Length ========== 
print("\nStep 3: Calculating fault length...")

for fault in fault_data:
    start_coord = (fault['start_lat'], fault['start_lon'])
    end_coord = (fault['end_lat'], fault['end_lon'])
    
    # Use geodesic to calculate great-circle distance (km)
    distance = geodesic(start_coord, end_coord).kilometers
    fault['length_km'] = distance
    
    # Also calculate the actual length along the segment (more accurate)
    geom = fault['geometry']
    if geom.geom_type == 'LineString':
        coords = list(geom.coords)
        total_length = 0
        for i in range(len(coords) - 1):
            p1 = (coords[i][1], coords[i][0])  # (lat, lon)
            p2 = (coords[i+1][1], coords[i+1][0])
            total_length += geodesic(p1, p2).kilometers
        fault['actual_length_km'] = total_length
    else:
        fault['actual_length_km'] = distance

lengths = [f['actual_length_km'] for f in fault_data]
print(f"Fault length range: {min(lengths):.2f} - {max(lengths):.2f} km")
print(f"Average length: {np.mean(lengths):.2f} km")

#%% 
# ========== Step 4: Estimate Maximum Possible Magnitude ========== 
print("\nStep 4: Estimating maximum possible magnitude...")

for fault in fault_data:
    length = fault['actual_length_km']
    
    # Wells & Coppersmith (1994) strike-slip fault relationship
    max_magnitude = 5.16 + 1.12 * np.log10(length)
    fault['max_magnitude'] = max_magnitude
    
    # Set reasonable magnitude range
    fault['min_magnitude'] = max(5.0, max_magnitude - 3.0)

magnitudes = [f['max_magnitude'] for f in fault_data]
print(f"Estimated maximum magnitude range: {min(magnitudes):.2f} - {max(magnitudes):.2f}")

#%% 
# ========== Step 5: Retrieve Historical Earthquake Data from USGS ========== 
print("\nStep 5: Retrieving historical earthquake data from USGS...")

# Define URL
url = "https://earthquake.usgs.gov/fdsnws/event/1/query.geojson?starttime=1900-07-12%2000:00:00&endtime=2025-07-19%2023:59:59&maxlatitude=47.093&minlatitude=30.44&maxlongitude=152.93&minlongitude=127.441&minmagnitude=5&orderby=time"

# Use requests.get() to fetch data
response = requests.get(url)

# Ensure the request was successful
if response.status_code == 200:
    # Get JSON data
    data = response.json()

    # Extract earthquake event data
    features = data['features']

    # Create a dictionary with earthquake information
    earthquake_data = []
    for feature in features:
        properties = feature['properties']
        geometry = feature['geometry']

        earthquake_info = {
            'id': feature['id'],
            'mag': properties['mag'],
            'place': properties['place'],
            'time': properties['time'],
            'updated': properties['updated'],
            'url': properties['url'],
            'detail': properties['detail'],
            'mmi': properties['mmi'],
            'alert': properties['alert'],
            'status': properties['status'],
            'tsunami': properties['tsunami'],
            'coordinates': geometry['coordinates'],
        }
        earthquake_data.append(earthquake_info)

    # Convert data to DataFrame
    earthquake_data = pd.DataFrame(earthquake_data)

    # Display first few rows
    print(earthquake_data.head())

else:
    print("Request failed, status code:", response.status_code)

#%% 
# ========== Step 6: Calculate GR Law Parameters ========== 
print("\nStep 6: Calculating GR law parameters...")

if earthquake_data is None or len(earthquake_data) == 0:
    print("No earthquake data, using typical a and b values")
    gr_parameters = {'a': 4.5, 'b': 1.0}
else:
    # Calculate observation time span (years)
    times = pd.to_datetime(earthquake_data['time'], unit='ms')
    time_span_years = (times.max() - times.min()).days / 365.25
    
    # Set magnitude range
    min_mag = earthquake_data['mag'].min()
    max_mag = earthquake_data['mag'].max()
    if pd.isna(min_mag) or pd.isna(max_mag):
        min_mag, max_mag = 5.0, 9.0
    else:
        min_mag, max_mag = float(min_mag), float(max_mag)
    magnitude_bins = np.arange(np.floor(min_mag * 2) / 2, 
                             np.ceil(max_mag * 2) / 2 + 0.1, 0.1)
    
    # Calculate cumulative counts
    magnitudes = []
    cumulative_counts = []
    annual_rates = []
    
    for mag in magnitude_bins:
        count = len(earthquake_data[earthquake_data['mag'] >= mag])
        annual_rate = count / time_span_years
        
        if count > 0:  # Only include magnitudes with earthquakes
            magnitudes.append(mag)
            cumulative_counts.append(count)
            annual_rates.append(annual_rate)
    
    if len(magnitudes) < 3:
        print("Not enough data points, using typical values")
        gr_parameters = {'a': 4.5, 'b': 1.0}
    else:
        # Linear regression fitting log10(N) = a - b*M
        magnitudes = np.array(magnitudes)
        log_rates = np.log10(np.array(annual_rates))
        
        # Filter out infinite values
        valid_indices = np.isfinite(log_rates)
        magnitudes_valid = magnitudes[valid_indices]
        log_rates_valid = log_rates[valid_indices]
        
        if len(magnitudes_valid) < 2:
            gr_parameters = {'a': 4.5, 'b': 1.0}
        else:
            # Linear fitting
            coeffs = np.polyfit(magnitudes_valid, log_rates_valid, 1)
            b_value = -coeffs[0]  # Negative of the slope
            a_value = coeffs[1]   # Intercept
            
            # Ensure parameters are within reasonable range
            b_value = max(0.5, min(2.0, b_value))
            a_value = max(0, min(10, a_value))
            
            gr_parameters = {
                'a': a_value,
                'b': b_value,
                'time_span_years': time_span_years,
                'n_earthquakes': len(earthquake_data)
            }
            
            plt.figure(figsize=(7,5),dpi=300)
            # Plot scatter points
            plt.scatter(magnitudes_valid, log_rates_valid, color='blue', label='Raw Data')
            # Plot fitted curve
            mag_fit = np.linspace(magnitudes_valid.min(), magnitudes_valid.max(), 100)
            log_rate_fit = a_value - b_value * mag_fit
            plt.plot(mag_fit, log_rate_fit, color='red', label='Fit Curve')
            plt.xlabel('M')
            plt.ylabel('log10(annual_rate)')
            plt.title('Gutenberg-Richter Para Fit')
            plt.legend()
            plt.grid(True, linestyle='--', alpha=0.5)
            plt.tight_layout()
            plt.show()

print(f"GR parameters: a = {gr_parameters.get('a', 'N/A'):.3f}, b = {gr_parameters.get('b', 'N/A'):.3f}")
if 'n_earthquakes' in gr_parameters:
    print(f"Based on {gr_parameters['n_earthquakes']} earthquake records, time span {gr_parameters.get('time_span_years', 0):.1f} years")
