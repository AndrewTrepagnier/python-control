# Interactive Plotting Tools

This module provides interactive plotting capabilities for the Python Control Systems Library using matplotlib.

## Root Locus GUI

The `root_locus_gui` function creates an interactive root locus plot with hover functionality, similar to MATLAB's root locus GUI.

### Features

- **Hover Information**: Hover over the root locus to see gain, damping ratio, and frequency
- **Original Plot Style**: Uses the same visual style as the original matplotlib root locus plots
- **Interactive Info Box**: Small info box in the corner shows real-time information
- **Poles and Zeros**: Visual display of open-loop poles and zeros
- **Customizable**: Various options for display and interaction

### Basic Usage

```python
import control as ct

# Create a system
s = ct.tf('s')
sys = 1 / (s**2 + 2*s + 1)

# Create interactive root locus plot
gui = ct.root_locus_gui(sys)
gui.show()
```

### Advanced Usage

```python
# Customize the plot
gui = ct.root_locus_gui(
    sys,
    title="My Root Locus",
    show_grid_lines=True,
    damping_lines=True,
    frequency_lines=True
)
gui.show()
```

### Parameters

- `sys`: LTI system (SISO only)
- `gains`: Custom gain range (optional)
- `xlim`, `ylim`: Axis limits (optional)
- `grid`: Show s-plane grid (default: True)
- `show_poles_zeros`: Show poles and zeros (default: True)
- `show_grid_lines`: Show grid lines (default: True)
- `damping_lines`: Show damping ratio lines (default: True)
- `frequency_lines`: Show frequency lines (default: True)
- `title`: Plot title

### Hover Information

When you hover over the root locus, you can see:

- **Gain**: The current gain value
- **Pole**: The pole location in the s-plane
- **Damping**: Damping ratio (for complex poles)
- **Frequency**: Natural frequency (for complex poles)

### Installation

The interactive tools require matplotlib:

```bash
pip install matplotlib
```

### Examples

See the `examples/` directory for more detailed examples:

- `simple_rlocus_gui_example.py`: Basic usage

### Comparison with MATLAB

This GUI provides similar functionality to MATLAB's root locus tool:

| Feature | MATLAB | Python Control |
|---------|--------|----------------|
| Hover information | ✓ | ✓ |
| Grid lines | ✓ | ✓ |
| Poles/zeros display | ✓ | ✓ |
| Custom gain ranges | ✓ | ✓ |
| Desktop application | ✓ | ✓ |
| Jupyter integration | ✗ | ✓ |

### Comparison with Existing Functionality

The python-control library already has some interactive features:

- **Original click functionality**: `ct.pole_zero_plot(rl_data, interactive=True)` allows clicking to see gain
- **This GUI adds**: Hover-based interaction (more intuitive) with real-time info box

### Troubleshooting

If you get an ImportError, make sure matplotlib is installed:

```bash
pip install matplotlib
```

For Jupyter notebooks, you may need to enable matplotlib rendering:

```python
%matplotlib inline
``` 