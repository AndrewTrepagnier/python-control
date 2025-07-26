"""
Interactive Root Locus GUI using Matplotlib.

This module provides an interactive root locus plot that allows users to hover
over the root locus to see gain, damping, and frequency information.
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
        
        # Set default limits if not specified
        if xlim is None and ylim is None:
            xlim = (-5, 2)
            ylim = (-3, 3)
        
        self.rl_data = root_locus_map(sys, gains=gains, xlim=xlim, ylim=ylim, **kwargs)
        
        # Initialize GUI elements
        self.fig = None
        self.ax = None
        self.info_text = None
        self.locus_lines = []
        self.cursor_marker = None
        
        # Precomputed high-resolution gain table
        self.gain_table = None
        self.gain_resolution = 10000  # High resolution for smooth interpolation
        
        self._create_plot()
        self._setup_interactivity()
        self._create_gain_table()
    
    def _create_gain_table(self):
        """Create a high-resolution precomputed table of gains and corresponding points."""
        
        if self.rl_data.loci is None or len(self.rl_data.gains) == 0:
            return
        
        # Create high-resolution gain array
        min_gain = np.min(self.rl_data.gains)
        max_gain = np.max(self.rl_data.gains)
        
        # Handle edge cases where min_gain might be zero or very small
        if min_gain <= 0:
            min_gain = 1e-6  # Small positive value
        
        if max_gain <= min_gain:
            max_gain = min_gain * 10  # Ensure we have a range
        
        # Use log spacing for better resolution at lower gains
        self.gain_table = {
            'gains': np.logspace(np.log10(min_gain), np.log10(max_gain), self.gain_resolution),
            'curves': []  # Store each locus as a separate curve
        }
        
        # Extract each locus as a separate curve for smooth interpolation
        num_loci = self.rl_data.loci.shape[1]
        
        for locus_idx in range(num_loci):
            curve_points = []
            curve_gains = []
            
            # Extract valid points for this locus
            for gain_idx, gain in enumerate(self.rl_data.gains):
                point = self.rl_data.loci[gain_idx, locus_idx]
                if point is not None and not np.isnan(point):
                    curve_points.append(point)
                    curve_gains.append(gain)
            
            if len(curve_points) > 3:  # Need at least 4 points for Catmull-Rom
                self.gain_table['curves'].append({
                    'points': np.array(curve_points),
                    'gains': np.array(curve_gains),
                    'lengths': self._compute_curve_lengths(curve_points)
                })
    
    def _compute_curve_lengths(self, points):
        """Compute cumulative arc lengths along the curve."""
        if len(points) < 2:
            return [0.0]
        
        lengths = [0.0]
        for i in range(1, len(points)):
            segment_length = abs(points[i] - points[i-1])
            lengths.append(lengths[-1] + segment_length)
        
        return np.array(lengths)
    
    def _find_closest_point_high_res(self, x, y):
        """Find the closest point using curve-following interpolation."""
        
        if self.gain_table is None or len(self.gain_table['curves']) == 0:
            return self._find_closest_point(x, y)
        
        target_point = complex(x, y)
        min_distance = float('inf')
        best_interpolated_point = None
        best_interpolated_gain = None
        
        # Check each curve
        for curve in self.gain_table['curves']:
            points = curve['points']
            gains = curve['gains']
            
            # Find the closest point on this curve
            distances = np.abs(points - target_point)
            closest_idx = np.argmin(distances)
            min_curve_distance = distances[closest_idx]
            
            if min_curve_distance < min_distance:
                min_distance = min_curve_distance
                
                # If we're close enough to this curve, interpolate along it
                if min_curve_distance < 10.0 and len(points) >= 4:
                    # Find the best interpolation point along the curve
                    interpolated_point, interpolated_gain = self._interpolate_along_curve(
                        target_point, points, gains, closest_idx
                    )
                    
                    if interpolated_point is not None:
                        best_interpolated_point = interpolated_point
                        best_interpolated_gain = interpolated_gain
        
        return best_interpolated_point, best_interpolated_gain
    
    def _interpolate_along_curve(self, target_point, points, gains, closest_idx):
        """Interpolate along a curve using Catmull-Rom splines."""
        
        if len(points) < 4:
            return points[closest_idx], gains[closest_idx]
        
        # Find the best segment for interpolation
        best_t = 0.0
        best_distance = float('inf')
        
        # Try interpolation in different segments around the closest point
        for start_idx in range(max(0, closest_idx - 2), min(len(points) - 3, closest_idx + 1)):
            if start_idx + 3 >= len(points):
                continue
            
            # Get 4 consecutive points for Catmull-Rom
            p0, p1, p2, p3 = points[start_idx:start_idx + 4]
            g0, g1, g2, g3 = gains[start_idx:start_idx + 4]
            
            # Try different interpolation parameters
            for t in np.linspace(0, 1, 50):  # 50 samples per segment
                # Interpolate the point
                interpolated_point = self._catmull_rom_interpolate(t, p0, p1, p2, p3)
                
                # Interpolate the gain
                interpolated_gain = self._catmull_rom_interpolate(t, g0, g1, g2, g3)
                
                # Check distance to target
                distance = abs(interpolated_point - target_point)
                
                if distance < best_distance:
                    best_distance = distance
                    best_t = t
                    best_point = interpolated_point
                    best_gain = interpolated_gain
        
        if best_distance < 10.0:
            return best_point, best_gain
        
        return points[closest_idx], gains[closest_idx]
    
    def _catmull_rom_interpolate(self, t, y0, y1, y2, y3):
        """Catmull-Rom spline interpolation between four points."""
        
        t2 = t * t
        t3 = t2 * t
        
        # Catmull-Rom coefficients
        p0 = -0.5 * t3 + t2 - 0.5 * t
        p1 = 1.5 * t3 - 2.5 * t2 + 1.0
        p2 = -1.5 * t3 + 2.0 * t2 + 0.5 * t
        p3 = 0.5 * t3 - 0.5 * t2
        
        return y0 * p0 + y1 * p1 + y2 * p2 + y3 * p3
    
    def _create_plot(self):
        """Create the root locus plot."""
        
        if self.title is None:
            if self.rl_data.sysname:
                title = f"Root Locus: {self.rl_data.sysname}"
            else:
                title = "Root Locus"
        else:
            title = self.title
        
        self.cplt = root_locus_plot(self.rl_data, grid=self.grid, title=title)
        
        self.fig = self.cplt.figure
        self.ax = self.cplt.axes[0, 0]
        
        if hasattr(self.cplt, 'lines') and len(self.cplt.lines) > 0:
            if len(self.cplt.lines.shape) > 1 and self.cplt.lines.shape[1] > 2:
                self.locus_lines = self.cplt.lines[0, 2]
        
        self._create_info_box()
        self._create_cursor_marker()
    
    def _create_info_box(self):
        """Create the information display box."""
        
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
        """Create the cursor marker."""
        
        self.cursor_marker, = self.ax.plot(
            [], [], 'go',
            markersize=8,
            markeredgecolor='darkgreen',
            markeredgewidth=1.5,
            markerfacecolor='lime',
            alpha=0.8,
            zorder=10
        )
        
        self.cursor_marker.set_visible(False)
    
    def _setup_interactivity(self):
        """Set up mouse event handlers."""
        
        self.fig.canvas.mpl_connect('motion_notify_event', self._on_mouse_move)
        self.fig.canvas.mpl_connect('axes_leave_event', self._on_mouse_leave)
    
    def _on_mouse_move(self, event):
        """Handle mouse movement events."""
        
        if event.inaxes != self.ax:
            self._hide_info_box()
            self._hide_cursor_marker()
            return
        
        closest_point, closest_gain = self._find_closest_point_high_res(event.xdata, event.ydata)
        
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
        
        for i, gain in enumerate(self.rl_data.gains):
            for j, locus in enumerate(self.rl_data.loci[i, :]):
                s = locus
                distance = np.sqrt((s.real - x)**2 + (s.imag - y)**2)
                
                if distance < min_distance:
                    min_distance = distance
                    closest_point = s
                    closest_gain = gain
                    closest_indices = (i, j)
        
        if min_distance < 10.0:
            if closest_indices is not None:
                interpolated_point, interpolated_gain = self._interpolate_point(x, y, closest_indices)
                if interpolated_point is not None:
                    return interpolated_point, interpolated_gain
            
            return closest_point, closest_gain
        
        return None, None
    
    def _interpolate_point(self, x, y, closest_indices):
        """Interpolate between nearby points for smoother movement."""
        
        i, j = closest_indices
        
        neighbors = []
        gains = []
        
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
        
        distances = [np.sqrt((n.real - x)**2 + (n.imag - y)**2) for n in neighbors]
        sorted_indices = np.argsort(distances)
        
        p1 = neighbors[sorted_indices[0]]
        p2 = neighbors[sorted_indices[1]]
        g1 = gains[sorted_indices[0]]
        g2 = gains[sorted_indices[1]]
        
        d1 = distances[sorted_indices[0]]
        d2 = distances[sorted_indices[1]]
        
        if d1 + d2 == 0:
            return p1, g1
        
        w1 = d2 / (d1 + d2)
        w2 = d1 / (d1 + d2)
        
        interpolated_point = w1 * p1 + w2 * p2
        interpolated_gain = w1 * g1 + w2 * g2
        
        return interpolated_point, interpolated_gain
    
    def _update_info_box(self, s, gain):
        """Update the information box with current point data."""
        
        if s is None or gain is None:
            return
        
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
        
        self.info_text.set_text(info_text)
        self.info_text.set_visible(True)
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
        
        self.cursor_marker.set_data([s.real], [s.imag])
        self.cursor_marker.set_visible(True)
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
    Create an interactive root locus GUI.
    
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
    """
    
    return RootLocusGUI(sys, **kwargs)


def rlocus_gui(sys: LTI, **kwargs) -> RootLocusGUI:
    """
    Convenience function for creating root locus GUI.
    
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
