"""
Interactive Root Locus GUI using Matplotlib.

This module provides an interactive root locus plot using matplotlib that allows
users to hover over the root locus to see how gain changes, similar to
MATLAB's root locus GUI.

Author: [Your Name]
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.widgets import TextBox
import warnings
from typing import Optional, Union, List, Tuple

from ..rlocus import root_locus_map, root_locus_plot
from ..pzmap import _find_root_locus_gain, _create_root_locus_label
from ..lti import LTI
from ..config import _get_param


class RootLocusGUI:
    """Interactive root locus GUI using matplotlib."""
    
    def __init__(self, sys: LTI, 
                 gains: Optional[np.ndarray] = None,
                 xlim: Optional[Tuple[float, float]] = None,
                 ylim: Optional[Tuple[float, float]] = None,
                 grid: bool = True,
                 show_poles_zeros: bool = True,
                 show_grid_lines: bool = True,
                 damping_lines: bool = True,
                 frequency_lines: bool = True,
                 title: Optional[str] = None,
                 **kwargs):
        """
        Initialize the interactive root locus GUI.
        
        Parameters
        ----------
        sys : LTI
            Linear time-invariant system (SISO only)
        gains : array_like, optional
            Gains to use in computing the root locus. If not given, gains are
            automatically chosen to include the main features.
        xlim : tuple, optional
            Limits for the x-axis (real part)
        ylim : tuple, optional
            Limits for the y-axis (imaginary part)
        grid : bool, optional
            If True, show the s-plane grid with damping and frequency lines
        show_poles_zeros : bool, optional
            If True, show the open-loop poles and zeros
        show_grid_lines : bool, optional
            If True, show the grid lines for damping and frequency
        damping_lines : bool, optional
            If True, show lines of constant damping ratio
        frequency_lines : bool, optional
            If True, show lines of constant frequency
        title : str, optional
            Title for the plot
        **kwargs
            Additional arguments passed to root_locus_map
        """
        
        if not sys.issiso():
            raise ValueError("System must be single-input single-output (SISO)")
        
        self.sys = sys
        self.gains = gains
        self.xlim = xlim
        self.ylim = ylim
        self.grid = grid
        self.show_poles_zeros = show_poles_zeros
        self.show_grid_lines = show_grid_lines
        self.damping_lines = damping_lines
        self.frequency_lines = frequency_lines
        self.title = title
        self.kwargs = kwargs
        
        # Get root locus data with wider limits if not specified
        if xlim is None and ylim is None:
            # Use wider default limits to allow the green dot to follow beyond ±1
            xlim = (-5, 2)  # Wider x range
            ylim = (-3, 3)  # Wider y range to go beyond ±1
        
        self.rl_data = root_locus_map(sys, gains=gains, xlim=xlim, ylim=ylim, **kwargs)
        
        # Initialize GUI elements
        self.fig = None
        self.ax = None
        self.info_text = None
        self.locus_lines = []
        self.cursor_marker = None  # Green dot that follows the cursor
        
        # Create the plot using the original root locus plotting
        self._create_plot()
        self._setup_interactivity()
    
    def _create_plot(self):
        """Create the root locus plot using the original plotting function."""
        
        # Use the original root locus plotting function
        if self.title is None:
            if self.rl_data.sysname:
                title = f"Root Locus: {self.rl_data.sysname}"
            else:
                title = "Root Locus"
        else:
            title = self.title
        
        # Create the plot using the original function
        self.cplt = root_locus_plot(self.rl_data, grid=self.grid, title=title)
        
        # Get the figure and axis
        self.fig = self.cplt.figure
        self.ax = self.cplt.axes[0, 0]  # Get the main axis
        
        # Store the locus lines for hover detection
        if hasattr(self.cplt, 'lines') and len(self.cplt.lines) > 0:
            # The locus lines are typically in the third column (index 2)
            if len(self.cplt.lines.shape) > 1 and self.cplt.lines.shape[1] > 2:
                self.locus_lines = self.cplt.lines[0, 2]  # First system, locus lines
        
        # Create info box
        self._create_info_box()
        
        # Create cursor marker (initially hidden)
        self._create_cursor_marker()
    
    def _create_info_box(self):
        """Create the information display box."""
        
        # Create text for the info box in the upper left corner
        self.info_text = self.ax.text(
            0.02, 0.98, "Hover over root locus\nto see gain information",
            transform=self.ax.transAxes,
            fontsize=10,
            verticalalignment='top',
            bbox=dict(
                boxstyle="round,pad=0.3", 
                facecolor='lightblue', 
                alpha=0.9,
                edgecolor='black',
                linewidth=1
            )
        )
    
    def _create_cursor_marker(self):
        """Create the cursor marker (green dot) that follows the mouse."""
        
        # Create a small green dot marker that will follow the cursor
        self.cursor_marker, = self.ax.plot(
            [], [], 'go',  # Green circle marker
            markersize=8,
            markeredgecolor='darkgreen',
            markeredgewidth=1.5,
            markerfacecolor='lime',
            alpha=0.8,
            zorder=10  # Ensure it's on top of other elements
        )
        
        # Initially hide the marker
        self.cursor_marker.set_visible(False)
    
    def _setup_interactivity(self):
        """Set up mouse event handlers."""
        
        # Connect mouse motion event
        self.fig.canvas.mpl_connect('motion_notify_event', self._on_mouse_move)
        
        # Connect mouse leave event
        self.fig.canvas.mpl_connect('axes_leave_event', self._on_mouse_leave)
    
    def _on_mouse_move(self, event):
        """Handle mouse movement events."""
        
        if event.inaxes != self.ax:
            self._hide_info_box()
            self._hide_cursor_marker()
            return
        
        # Find the closest point on the root locus
        closest_point, closest_gain = self._find_closest_point(event.xdata, event.ydata)
        
        if closest_point is not None:
            self._update_info_box(closest_point, closest_gain)
            self._update_cursor_marker(closest_point)
        else:
            self._hide_info_box()
            self._hide_cursor_marker()
    
    def _on_mouse_leave(self, event):
        """Handle mouse leave events."""
        self._hide_info_box()
        self._hide_cursor_marker()
    
    def _find_closest_point(self, x, y):
        """Find the closest point on the root locus to the given coordinates."""
        
        if self.rl_data.loci is None:
            return None, None
        
        min_distance = float('inf')
        closest_point = None
        closest_gain = None
        closest_indices = None
        
        # Search through all locus points
        for i, gain in enumerate(self.rl_data.gains):
            for j, locus in enumerate(self.rl_data.loci[i, :]):
                s = locus
                distance = np.sqrt((s.real - x)**2 + (s.imag - y)**2)
                
                if distance < min_distance:
                    min_distance = distance
                    closest_point = s
                    closest_gain = gain
                    closest_indices = (i, j)
        
        # Use a very large threshold to make the green dot always responsive
        if min_distance < 10.0:  # Very large threshold for maximum responsiveness
            # If we found a close point, try to interpolate for smoother movement
            if closest_indices is not None:
                interpolated_point, interpolated_gain = self._interpolate_point(x, y, closest_indices)
                if interpolated_point is not None:
                    return interpolated_point, interpolated_gain
            
            return closest_point, closest_gain
        
        return None, None
    
    def _interpolate_point(self, x, y, closest_indices):
        """Interpolate between nearby points for smoother movement."""
        
        i, j = closest_indices
        
        # Get neighboring points for interpolation
        neighbors = []
        gains = []
        
        # Check points around the closest point
        for di in [-1, 0, 1]:
            for dj in [-1, 0, 1]:
                ni, nj = i + di, j + dj
                if (0 <= ni < len(self.rl_data.gains) and 
                    0 <= nj < self.rl_data.loci.shape[1]):
                    neighbor = self.rl_data.loci[ni, nj]
                    if neighbor is not None and not np.isnan(neighbor):
                        neighbors.append(neighbor)
                        gains.append(self.rl_data.gains[ni])
        
        if len(neighbors) < 2:
            return None, None
        
        # Find the two closest neighbors to the mouse position
        distances = [np.sqrt((n.real - x)**2 + (n.imag - y)**2) for n in neighbors]
        sorted_indices = np.argsort(distances)
        
        # Get the two closest points
        p1 = neighbors[sorted_indices[0]]
        p2 = neighbors[sorted_indices[1]]
        g1 = gains[sorted_indices[0]]
        g2 = gains[sorted_indices[1]]
        
        # Calculate interpolation weight based on distance
        d1 = distances[sorted_indices[0]]
        d2 = distances[sorted_indices[1]]
        
        if d1 + d2 == 0:
            return p1, g1
        
        # Weighted interpolation
        w1 = d2 / (d1 + d2)
        w2 = d1 / (d1 + d2)
        
        # Interpolate the complex point
        interpolated_point = w1 * p1 + w2 * p2
        
        # Interpolate the gain
        interpolated_gain = w1 * g1 + w2 * g2
        
        return interpolated_point, interpolated_gain
    
    def _update_info_box(self, s, gain):
        """Update the information box with current point data."""
        
        if s is None or gain is None:
            return
        
        # Calculate damping ratio and frequency
        if s.imag != 0:
            wn = abs(s)
            zeta = -s.real / wn
            info_text = f"Gain: {gain:.3f}\n"
            info_text += f"Pole: {s:.3f}\n"
            info_text += f"Damping: {zeta:.3f}\n"
            info_text += f"Frequency: {wn:.3f} rad/s"
        else:
            info_text = f"Gain: {gain:.3f}\n"
            info_text += f"Pole: {s:.3f}\n"
            info_text += "Real pole"
        
        # Update the text
        self.info_text.set_text(info_text)
        
        # Make sure the text is visible
        self.info_text.set_visible(True)
        
        # Redraw
        self.fig.canvas.draw_idle()
    
    def _hide_info_box(self):
        """Hide the information box."""
        
        self.info_text.set_visible(False)
        self.fig.canvas.draw_idle()
    
    def _update_cursor_marker(self, s):
        """Update the cursor marker position."""
        
        if s is None:
            self._hide_cursor_marker()
            return
        
        # Update the marker position to the closest point on the root locus
        self.cursor_marker.set_data([s.real], [s.imag])
        self.cursor_marker.set_visible(True)
        
        # Redraw
        self.fig.canvas.draw_idle()
    
    def _hide_cursor_marker(self):
        """Hide the cursor marker."""
        
        self.cursor_marker.set_visible(False)
        self.fig.canvas.draw_idle()
    
    def show(self):
        """Show the interactive plot."""
        plt.show()
    
    def save(self, filename, **kwargs):
        """Save the plot to a file."""
        self.fig.savefig(filename, **kwargs)


def root_locus_gui(sys: LTI, **kwargs) -> RootLocusGUI:
    """
    Create an interactive root locus GUI using matplotlib.
    
    Parameters
    ----------
    sys : LTI
        Linear time-invariant system (SISO only)
    **kwargs
        Additional arguments passed to RootLocusGUI
        
    Returns
    -------
    RootLocusGUI
        Interactive root locus GUI object
        
    Examples
    --------
    >>> import control as ct
    >>> import numpy as np
    >>> 
    >>> # Create a simple system
    >>> s = ct.tf('s')
    >>> sys = (s + 1) / (s**2 + 2*s + 1)
    >>> 
    >>> # Create interactive root locus GUI
    >>> gui = ct.root_locus_gui(sys)
    >>> gui.show()
    """
    
    return RootLocusGUI(sys, **kwargs)


# Convenience function for quick plotting
def rlocus_gui(sys: LTI, **kwargs) -> RootLocusGUI:
    """
    Convenience function for creating root locus GUI.
    
    This is a shorthand for root_locus_gui().
    
    Parameters
    ----------
    sys : LTI
        Linear time-invariant system (SISO only)
    **kwargs
        Additional arguments passed to root_locus_gui
        
    Returns
    -------
    RootLocusGUI
        Interactive root locus GUI object
    """
    return root_locus_gui(sys, **kwargs)


# Keep the advanced function for future implementation
def root_locus_gui_advanced(sys: LTI, **kwargs):
    """
    Advanced root locus GUI with additional features.
    
    This version includes:
    - Multiple subplots (root locus + step response)
    - Real-time gain adjustment
    - System information panel
    
    Parameters
    ----------
    sys : LTI
        Linear time-invariant system (SISO only)
    **kwargs
        Additional arguments passed to root_locus_gui
        
    Returns
    -------
    RootLocusGUI
        Interactive root locus GUI object
    """
    
    # For now, just return the basic GUI
    # TODO: Implement advanced features
    return root_locus_gui(sys, **kwargs)
