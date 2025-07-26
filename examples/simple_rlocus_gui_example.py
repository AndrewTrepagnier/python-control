#!/usr/bin/env python3
"""
Simple example demonstrating the Root Locus GUI.

This example shows how to create an interactive root locus plot
with hover functionality.
"""

import numpy as np
import control as ct

def main():
    """Run a simple example of the root locus GUI."""
    
    print("Root Locus GUI - Simple Example")
    print("=" * 40)
    
    try:
        s = ct.tf('s')
        sys = 1 / (s**2 + 2*s + 1)
        
        print(f"System: {sys}")
        print("Creating interactive root locus plot...")
        
        gui = ct.root_locus_gui(sys, 
                               title="Simple Root Locus Example",
                               show_grid_lines=True)
        
        print("Displaying plot...")
        print("Hover over the root locus curves to see gain, damping, and frequency information.")
        gui.show()
        
        print("\nExample completed successfully!")
        
    except ImportError as e:
        print(f"Error: {e}")
        print("\nTo use this example, make sure matplotlib is installed:")
        print("pip install matplotlib")
    except Exception as e:
        print(f"Error during example: {e}")

if __name__ == "__main__":
    main() 