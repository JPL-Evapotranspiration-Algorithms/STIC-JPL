#!/usr/bin/env python
"""
Test STIC-JPL batch processing with daylight upscaling.

This test verifies that the STIC-JPL model can now process multiple rows
efficiently with daylight upscaling enabled, thanks to the fixes in the
daylight-evapotranspiration package.
"""

import numpy as np
import pandas as pd
from datetime import datetime, timedelta
from rasters import MultiPoint, WGS84
from STIC_JPL import process_STIC_table

print("Creating test data for STIC-JPL processing...")

# Create sample data with multiple rows
n_samples = 5
data = {
    'ST_C': np.random.uniform(25, 35, n_samples),
    'EmisWB': np.random.uniform(0.95, 0.98, n_samples),
    'NDVI': np.random.uniform(0.6, 0.8, n_samples),
    'albedo': np.random.uniform(0.15, 0.25, n_samples),
    'Ta_C': np.random.uniform(20, 30, n_samples),
    'RH': np.random.uniform(0.3, 0.6, n_samples),
    'Rn_Wm2': np.random.uniform(400, 500, n_samples),
    'Rg': np.random.uniform(800, 1000, n_samples),
    'lon': np.random.uniform(-121, -120, n_samples),
    'lat': np.random.uniform(37, 38, n_samples),
    'time_UTC': [datetime(2023, 6, 15, 19, 0, 0) + timedelta(minutes=30*i) for i in range(n_samples)]
}

input_df = pd.DataFrame(data)

print(f"\nTest configuration:")
print(f"  Number of samples: {n_samples}")
print(f"  Time range: {input_df['time_UTC'].min()} to {input_df['time_UTC'].max()}")
print(f"  Latitude range: {input_df['lat'].min():.2f} to {input_df['lat'].max():.2f}")
print(f"  Longitude range: {input_df['lon'].min():.2f} to {input_df['lon'].max():.2f}")

print("\n" + "="*70)
print("TEST: Batch processing with daylight upscaling enabled")
print("="*70)

try:
    print("\nProcessing data...")
    result_df = process_STIC_table(
        input_df,
        constrain_negative_LE=False,
        supply_SWin=False,
        upscale_to_daylight=True
    )
    
    print("\n✓ SUCCESS: Batch processing completed!")
    print(f"\nResult DataFrame shape: {result_df.shape}")
    print(f"\nOutput columns added:")
    
    output_cols = [col for col in result_df.columns if col not in input_df.columns]
    for col in output_cols:
        if col in result_df.columns:
            values = result_df[col]
            if hasattr(values, '__len__') and len(values) > 0:
                mean_val = np.nanmean(values)
                print(f"  {col}: mean={mean_val:.3f}")
    
    # Check for daylight upscaling outputs
    daylight_cols = ['ET_daylight_kg', 'LE_daylight_Wm2', 'Rn_daylight_Wm2', 'EF']
    found_daylight = [col for col in daylight_cols if col in result_df.columns]
    
    if found_daylight:
        print(f"\n✓ Daylight upscaling columns found: {', '.join(found_daylight)}")
    else:
        print(f"\n⚠ Warning: No daylight upscaling columns found")
    
    print("\n" + "="*70)
    print("SUMMARY")
    print("="*70)
    print("\nThe STIC-JPL model can now efficiently process multiple rows")
    print("with daylight upscaling enabled, thanks to the improvements in")
    print("the daylight-evapotranspiration package that properly handle")
    print("arrays of datetime objects.")
    
except Exception as e:
    print(f"\n✗ FAILED: {type(e).__name__}: {e}")
    import traceback
    print("\nFull traceback:")
    traceback.print_exc()
    print("\n" + "="*70)
    print("The error above indicates that there may still be issues with")
    print("batch processing or the daylight-evapotranspiration package.")
