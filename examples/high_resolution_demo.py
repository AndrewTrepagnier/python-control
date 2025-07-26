#!/usr/bin/env python3
"""
High-Resolution Root Locus GUI Demo.

This example demonstrates the precomputed gain table with 10,000 resolution points
and Catmull-Rom spline interpolation for ultra-smooth green dot movement.
"""

import control as ct
import numpy as np

def main():
    """Demonstrate the high-resolution root locus GUI."""
    
    print("High-Resolution Root Locus GUI Demo")
    print("=" * 40)
    
    try:
        # Create a complex system with multiple asymptotes
        s = ct.tf('s')
        sys = (s**2 + 2*s + 2) / (s**4 + 5*s**3 + 8*s**2 + 6*s + 2)
        
        print(f"System: {sys}")
        print()
        print("Features:")
        print("- Precomputed gain table with 10,000 resolution points")
        print("- Catmull-Rom spline interpolation for ultra-smooth movement")
        print("- Log-spaced gains for better resolution at lower gains")
        print("- Green dot should slide like butter along the curves!")
        print()
        
        # Create the high-resolution GUI
        gui = ct.root_locus_gui(sys, title="High-Resolution Root Locus")
        
        # Show info about the gain table
        if gui.gain_table is not None:
            print(f"Gain table created with {len(gui.gain_table['gains'])} points")
            print(f"Gain range: {gui.gain_table['gains'][0]:.2e} to {gui.gain_table['gains'][-1]:.2e}")
            print(f"Number of curves: {len(gui.gain_table['curves'])}")
            for i, curve in enumerate(gui.gain_table['curves']):
                print(f"  Curve {i}: {len(curve['points'])} points")
        else:
            print("Gain table creation failed")
        
        print()
        print("Displaying plot...")
        print("Move your mouse over the root locus for ultra-smooth green dot movement!")
        
        gui.show()
        
        print("\nDemo completed!")
        
    except Exception as e:
        print(f"Error during demo: {e}")
        import traceback
        traceback.print_exc()

if __name__ == "__main__":
    main() 