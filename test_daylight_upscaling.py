#!/usr/bin/env python
"""
Test script to reproduce the daylight upscaling issue with arrays of datetime objects.

This script isolates the problem that occurs when daylight_ET_from_instantaneous_LE
receives arrays of datetime objects instead of single datetime values.

Run this script to see the TypeError that needs to be fixed in the daylight-evapotranspiration package.
"""

import numpy as np
import pandas as pd
from datetime import datetime
from rasters import MultiPoint, WGS84

# Import the function that's causing issues
from daylight_evapotranspiration import daylight_ET_from_instantaneous_LE

# Create test data that mimics what STIC-JPL passes
print("Creating test data...")

# Sample data for 3 locations/times
n_samples = 3

# Create sample flux data
LE_instantaneous_Wm2 = np.array([150.5, 200.3, 180.7])
Rn_instantaneous_Wm2 = np.array([400.2, 450.8, 420.5])
G_instantaneous_Wm2 = np.array([50.1, 60.2, 55.3])

# Create sample geometry (3 points)
lons = np.array([-120.5, -119.8, -121.2])
lats = np.array([37.5, 38.2, 36.8])
geometry = MultiPoint(x=lons, y=lats, crs=WGS84)

# Create sample times as a list of datetime objects
time_UTC_list = [
    datetime(2023, 6, 15, 18, 30, 0),  # ~10:30 AM local time in California
    datetime(2023, 6, 15, 19, 0, 0),   # ~11:00 AM local time
    datetime(2023, 6, 15, 19, 30, 0),  # ~11:30 AM local time
]

# Calculate day of year for reference
day_of_year = np.array([166, 166, 166])

print(f"\nTest configuration:")
print(f"  Number of samples: {n_samples}")
print(f"  LE range: {LE_instantaneous_Wm2.min():.1f} - {LE_instantaneous_Wm2.max():.1f} W/m²")
print(f"  Rn range: {Rn_instantaneous_Wm2.min():.1f} - {Rn_instantaneous_Wm2.max():.1f} W/m²")
print(f"  G range: {G_instantaneous_Wm2.min():.1f} - {G_instantaneous_Wm2.max():.1f} W/m²")
print(f"  Longitude range: {lons.min():.2f} - {lons.max():.2f}")
print(f"  Latitude range: {lats.min():.2f} - {lats.max():.2f}")
print(f"  Times: {time_UTC_list[0]} to {time_UTC_list[-1]}")

# Test 1: Single datetime (this should work)
print("\n" + "="*70)
print("TEST 1: Single datetime value (expected to work)")
print("="*70)
try:
    result_single = daylight_ET_from_instantaneous_LE(
        LE_instantaneous_Wm2=LE_instantaneous_Wm2[0],
        Rn_instantaneous_Wm2=Rn_instantaneous_Wm2[0],
        G_instantaneous_Wm2=G_instantaneous_Wm2[0],
        day_of_year=day_of_year[0],
        time_UTC=time_UTC_list[0],
        geometry=geometry
    )
    print("✓ SUCCESS: Single datetime processed correctly")
    print(f"  Result keys: {list(result_single.keys())}")
    for key, value in result_single.items():
        if isinstance(value, (int, float, np.number)):
            print(f"    {key}: {value:.3f}")
        elif hasattr(value, 'shape'):
            print(f"    {key}: array with shape {value.shape}, mean={np.nanmean(value):.3f}")
except Exception as e:
    print(f"✗ FAILED: {type(e).__name__}: {e}")

# Test 2: List of datetime objects (this will fail with current implementation)
print("\n" + "="*70)
print("TEST 2: List of datetime objects (expected to fail)")
print("="*70)
print("This is the problematic case that needs to be fixed...")
try:
    result_list = daylight_ET_from_instantaneous_LE(
        LE_instantaneous_Wm2=LE_instantaneous_Wm2,
        Rn_instantaneous_Wm2=Rn_instantaneous_Wm2,
        G_instantaneous_Wm2=G_instantaneous_Wm2,
        day_of_year=day_of_year,
        time_UTC=time_UTC_list,
        geometry=geometry
    )
    print("✓ SUCCESS: List of datetimes processed correctly")
    print(f"  Result keys: {list(result_list.keys())}")
    for key, value in result_list.items():
        if hasattr(value, 'shape'):
            print(f"    {key}: array with shape {value.shape}, mean={np.nanmean(value):.3f}")
        else:
            print(f"    {key}: {value}")
except Exception as e:
    print(f"✗ FAILED: {type(e).__name__}: {e}")
    print(f"\n  Full error traceback:")
    import traceback
    traceback.print_exc()

# Test 3: Array of datetime objects (also likely to fail)
print("\n" + "="*70)
print("TEST 3: NumPy array of datetime objects")
print("="*70)
try:
    time_UTC_array = np.array(time_UTC_list)
    result_array = daylight_ET_from_instantaneous_LE(
        LE_instantaneous_Wm2=LE_instantaneous_Wm2,
        Rn_instantaneous_Wm2=Rn_instantaneous_Wm2,
        G_instantaneous_Wm2=G_instantaneous_Wm2,
        day_of_year=day_of_year,
        time_UTC=time_UTC_array,
        geometry=geometry
    )
    print("✓ SUCCESS: NumPy array of datetimes processed correctly")
    print(f"  Result keys: {list(result_array.keys())}")
    for key, value in result_array.items():
        if hasattr(value, 'shape'):
            print(f"    {key}: array with shape {value.shape}, mean={np.nanmean(value):.3f}")
        else:
            print(f"    {key}: {value}")
except Exception as e:
    print(f"✗ FAILED: {type(e).__name__}: {e}")

print("\n" + "="*70)
print("SUMMARY")
print("="*70)
print("""
The issue occurs when passing lists or arrays of datetime objects to 
daylight_ET_from_instantaneous_LE. The function needs to be modified to:

1. Detect when time_UTC is a list or array of datetime objects
2. Loop through each datetime and corresponding geometry point
3. Call calculate_solar_hour_of_day and calculate_solar_day_of_year 
   for each element individually
4. Combine the results into arrays

The fix should be implemented in:
  daylight-evapotranspiration/daylight_evapotranspiration/daylight_evapotranspiration.py
  
Around lines 289-303 where hour_of_day and day_of_year are calculated.
""")
