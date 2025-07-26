#!/usr/bin/env python3
"""
Complex Root Locus GUI Example - Beautiful System with Multiple Asymptotes.

This example demonstrates the interactive root locus GUI with a complex system
that has multiple asymptotes and curves, showcasing the smooth green dot
cursor marker functionality.
"""

import control as ct
import numpy as np

def main():
    """Demonstrate the beautiful complex root locus GUI."""
    
    print("Complex Root Locus GUI - Beautiful System Demo")
    print("=" * 50)
    
    try:
        # Create a beautiful complex system with multiple asymptotes
        s = ct.tf('s')
        sys = (s**2 + 2*s + 2) / (s**4 + 5*s**3 + 8*s**2 + 6*s + 2)
        
        print("System created:")
        print(f"Numerator: {s**2 + 2*s + 2}")
        print(f"Denominator: {s**4 + 5*s**3 + 8*s**2 + 6*s + 2}")
        print()
        print("Features to explore:")
        print("- Multiple asymptotes and curves")
        print("- Smooth green dot cursor marker")
        print("- Real-time gain, damping, and frequency display")
        print("- Works beyond ±1 bounds")
        print("- Hover anywhere near the curves!")
        print()
        
        # Create the interactive GUI
        gui = ct.root_locus_gui(sys, title="Beautiful Complex Root Locus")
        
        # Show the plot
        gui.show()
        
        print("\nDemo completed! The green dot should slide smoothly along all curves.")
        
    except ImportError as e:
        print(f"Error: {e}")
        print("\nTo use this example, make sure matplotlib is installed:")
        print("pip install matplotlib")
    except Exception as e:
        print(f"Error during demo: {e}")

if __name__ == "__main__":
    main() 