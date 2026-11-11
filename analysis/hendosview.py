"""
HENDOS - HErmite N-body Dynamics Of Systems
The viewer: look at a run.

    python hendosview.py [data directory] [runid]

        data directory   where the run's files are (default: data)
        runid            which run (default: the first one found there;
                         the 'Next run' button steps through the others)

Three views, chosen at the top left:
    3D view            the stars of one snapshot. Drag with the mouse to
                       rotate. Slider or Play to move through the run.
    Snapshot analysis  a quantity worked out from every snapshot, against
                       time (the analyses are in hendosanalysis.py)
    Run-time files     the files the code wrote while it ran,
                       <runid>_ana_<name>.dat

Needs numpy and matplotlib (scipy as well for the Q parameter). The window
is made with matplotlib's own buttons, so there is nothing else to install.
Movies need ffmpeg (mp4); without it an animated GIF is written instead.

The drawing is done by two functions that do not need the window, so they
can be used in your own scripts:
    draw_particles(ax, snap, ...)     one snapshot on a 3D axis
    plot_curves(fig, rect, ...)       curves against time
"""
import os
import sys

import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib import animation
from matplotlib.colors import LinearSegmentedColormap, LogNorm, Normalize
from matplotlib.figure import Figure
from matplotlib.backends.backend_agg import FigureCanvasAgg
from matplotlib.lines import Line2D
from matplotlib.widgets import Button, CheckButtons, RadioButtons, Slider

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import hendosanalysis as ha                     # noqa: E402
from hendosio import Run, find_runs             # noqa: E402

# ----------------------------------------------------------------------
# colours. One set for a light background and one for a dark background.
#   series   up to three hues for separate things (checked to be told
#            apart by colour-blind readers)
#   ramp     one hue, light to dark, for 'how much'
# ----------------------------------------------------------------------
BLUES = ['#cde2fb', '#b7d3f6', '#9ec5f4', '#86b6ef', '#6da7ec', '#5598e7',
         '#3987e5', '#2a78d6', '#256abf', '#1c5cab', '#184f95', '#104281', '#0d366b']
LIGHT = dict(surface='#fcfcfb', panel='#f1f0ec', text='#0b0b0b', text2='#52514e',
             grid='#dddcd6', neutral='#8b8a84',
             series=['#2a78d6', '#eb6834', '#1baf7a'],
             ramp=BLUES[3:],             # low = light, but still clear on white
             ordinal=BLUES[3:])
DARK = dict(surface='#1a1a19', panel='#272726', text='#f0efec', text2='#c3c2b7',
            grid='#3a3a37', neutral='#8f8e87',
            series=['#3987e5', '#d95926', '#199e70'],
            ramp=BLUES[:11][::-1],       # low = dark, but still clear on black
            ordinal=BLUES[:11][::-1])

# how each particle type is drawn: (which colour, marker shape).
# Shape as well as colour, so types can be told apart without colour.
# TO ADD A TYPE: add a line here.
TYPE_STYLE = {
    'star':         ('neutral', 'o'),
    'white_dwarf':  ('series0', 'D'),
    'neutron_star': ('series1', 's'),
    'black_hole':   ('text',    'X'),
    'planet':       ('series2', '^'),
}
MASS_BINS = [0.5, 2.0, 8.0]                         # Msun: edges of the mass bins
COLOUR_MODES = ['mass', 'mass bins', 'type', 'speed', 'bound / unbound', 'timestep']


def theme_colour(theme, name):
    if name.startswith('series'):
        return theme['series'][int(name[-1])]
    return theme[name]


def ramp_cmap(theme):
    return LinearSegmentedColormap.from_list('hendos', theme['ramp'])


def ordinal_colours(theme, k):
    """k colours in order along the one-hue ramp."""
    steps = theme['ordinal']
    if k == 1:
        return [steps[len(steps) // 2]]
    return [steps[int(round(i * (len(steps) - 1) / (k - 1)))] for i in range(k)]


# ----------------------------------------------------------------------
# drawing one snapshot
# ----------------------------------------------------------------------
def draw_particles(ax, snap, colour='mass', theme=LIGHT, size_by_mass=True,
                   shape_by_type=False, ring_massive=False, centre=True,
                   halfwidth=None, point_size=None):
    """Draw the stars of one snapshot on a 3D axis.

    colour         one of COLOUR_MODES
    size_by_mass   bigger points for more massive stars
    shape_by_type  a different marker for each particle type (always done
                   when colouring by type)
    ring_massive   ring the ten most massive stars
    centre         put the centre of mass at the origin
    halfwidth      half the side of the box shown (pc); None = fit the stars
    point_size     area of a typical point; None = chosen from the number

    Returns what the key needs: ('bar', mappable, label) for a colour bar,
    or ('legend', handles, title) for a legend.
    """
    n = snap.n
    com, cov = ha.centre_of_mass(snap)
    pos = snap.pos - com if centre else snap.pos
    s0 = point_size if point_size is not None else float(np.clip(6000.0 / n, 3.0, 36.0))
    if size_by_mass:
        sizes = s0 * np.clip((snap.m / np.median(snap.m)) ** (2.0 / 3.0), 0.4, 30.0)
    else:
        sizes = np.full(n, s0)

    # ---- what colour is each star?
    values = None          # numbers for a colour bar
    colours = None         # or a colour for each star
    handles = []
    if colour == 'type':
        shape_by_type = True
        colours = np.empty(n, dtype=object)
        for t in np.unique(snap.type):
            cname, marker = TYPE_STYLE.get(t, ('neutral', 'P'))
            colours[snap.type == t] = theme_colour(theme, cname)
            handles.append(Line2D([], [], ls='', marker=marker, ms=7,
                                  color=theme_colour(theme, cname),
                                  label='%s (%d)' % (t.replace('_', ' '),
                                                     (snap.type == t).sum())))
        key = ('legend', handles, 'type')
    elif colour == 'mass bins':
        edges = MASS_BINS
        which = np.digitize(snap.m, edges)
        cols = ordinal_colours(theme, len(edges) + 1)
        colours = np.array(cols, dtype=object)[which]
        labels = (['< %g' % edges[0]]
                  + ['%g - %g' % (edges[i], edges[i + 1]) for i in range(len(edges) - 1)]
                  + ['> %g' % edges[-1]])
        for i, lab in enumerate(labels):
            handles.append(Line2D([], [], ls='', marker='o', ms=7, color=cols[i],
                                  label='%s (%d)' % (lab, (which == i).sum())))
        key = ('legend', handles, 'mass (Msun)')
    elif colour == 'bound / unbound':
        e = ha.specific_energy(snap)
        if e is None:
            colours = np.full(n, theme['neutral'], dtype=object)
            key = ('legend', [Line2D([], [], ls='', marker='o', ms=7,
                                     color=theme['neutral'],
                                     label='no potential in this file')], '')
        else:
            colours = np.where(e < 0, theme['series'][0], theme['series'][1]).astype(object)
            key = ('legend',
                   [Line2D([], [], ls='', marker='o', ms=7, color=theme['series'][0],
                           label='bound (%d)' % (e < 0).sum()),
                    Line2D([], [], ls='', marker='o', ms=7, color=theme['series'][1],
                           label='unbound (%d)' % (e >= 0).sum())],
                   'energy')
    else:
        if colour == 'speed':
            values, label, log = ha.speed(snap), 'speed (km/s)', False
        elif colour == 'timestep' and snap.dt is not None:
            values, label, log = snap.dt, 'timestep (yr)', True
        else:
            values, label, log = snap.m, 'mass (Msun)', True
        lo, hi = float(np.min(values)), float(np.max(values))
        if hi <= lo:
            hi = lo * 1.01 + 1e-30
        norm = LogNorm(lo, hi) if (log and lo > 0) else Normalize(lo, hi)
        key = ('bar', None, label)

    # ---- stars outside a fixed box are not drawn (3D axes do not clip)
    mid = np.zeros(3) if centre else com
    if halfwidth is None:
        halfwidth = 1.05 * np.percentile(np.sqrt(((pos - mid) ** 2).sum(1)), 98)
        inside = np.ones(n, bool)
    else:
        inside = (np.abs(pos - mid) <= halfwidth).all(1)

    # ---- draw, one marker shape at a time
    groups = [(None, 'o')]
    if shape_by_type:
        groups = [(t, TYPE_STYLE.get(t, ('neutral', 'P'))[1]) for t in np.unique(snap.type)]
    mappable = None
    for t, marker in groups:
        pick = inside if t is None else (inside & (snap.type == t))
        if not pick.any():
            continue
        # special types are few and small in mass: keep them visible
        sz = sizes[pick] if t in (None, 'star') else np.maximum(sizes[pick], 2.5 * s0)
        kw = dict(s=sz, marker=marker, depthshade=False,
                  edgecolors=theme['surface'], linewidths=0.3)
        if values is not None:
            mappable = ax.scatter(pos[pick, 0], pos[pick, 1], pos[pick, 2],
                                  c=values[pick], cmap=ramp_cmap(theme), norm=norm, **kw)
        else:
            ax.scatter(pos[pick, 0], pos[pick, 1], pos[pick, 2],
                       c=list(colours[pick]), **kw)
    if values is not None:
        key = ('bar', mappable, key[2])
        if shape_by_type:
            # shapes need a key as well as the colour bar
            key = key + ([Line2D([], [], ls='', marker=TYPE_STYLE.get(t, ('', 'P'))[1],
                                 ms=7, color=theme['text2'], label=t.replace('_', ' '))
                          for t in np.unique(snap.type)],)

    if ring_massive:
        top = np.argsort(snap.m)[::-1][:10]
        top = top[inside[top]]
        ax.scatter(pos[top, 0], pos[top, 1], pos[top, 2], s=np.maximum(sizes[top], s0) * 5,
                   facecolors='none', edgecolors=theme['text'], linewidths=1.1,
                   depthshade=False)

    # ---- the box
    ax.set_xlim(mid[0] - halfwidth, mid[0] + halfwidth)
    ax.set_ylim(mid[1] - halfwidth, mid[1] + halfwidth)
    ax.set_zlim(mid[2] - halfwidth, mid[2] + halfwidth)
    ax.set_box_aspect((1, 1, 1))
    style_3d(ax, theme)
    return key


def style_3d(ax, theme):
    """Quiet axes: the stars are the picture."""
    ax.set_facecolor(theme['surface'])
    for axis, name in ((ax.xaxis, 'x (pc)'), (ax.yaxis, 'y (pc)'), (ax.zaxis, 'z (pc)')):
        axis.set_pane_color(theme['surface'])
        axis.pane.set_edgecolor(theme['grid'])
        axis._axinfo['grid']['color'] = theme['grid']
        axis._axinfo['grid']['linewidth'] = 0.5
        axis.line.set_color(theme['grid'])
        axis.label.set_color(theme['text2'])
        axis.set_label_text(name)
        axis.set_tick_params(colors=theme['text2'], labelsize=8)


def draw_key(fig, cax, key, theme, anchor=(0.86, 0.62)):
    """The colour bar or legend that goes with draw_particles. Returns the legend (or None)."""
    legend = None
    cax.cla()
    if key[0] == 'bar':
        cax.set_visible(True)
        bar = fig.colorbar(key[1], cax=cax)
        bar.set_label(key[2], color=theme['text2'])
        bar.ax.tick_params(colors=theme['text2'], labelsize=8, which='both')
        bar.outline.set_edgecolor(theme['grid'])
        if len(key) > 3:
            legend = fig.legend(handles=key[3], title='type', loc='upper left',
                                bbox_to_anchor=(anchor[0], 0.26), frameon=False,
                                labelcolor=theme['text'], fontsize=9)
    else:
        cax.set_visible(False)
        legend = fig.legend(handles=key[1], title=key[2], loc='center left',
                            bbox_to_anchor=anchor, frameon=False,
                            labelcolor=theme['text'], fontsize=9)
    if legend is not None:
        legend.get_title().set_color(theme['text2'])
    return legend


# ----------------------------------------------------------------------
# curves against time
# ----------------------------------------------------------------------
def plot_curves(fig, rect, groups, theme=LIGHT, logy=False, tmark=None, title=''):
    """Plot curves against time, one panel for each group, stacked.

    rect     [left, bottom, width, height] of the whole plot in the figure
    groups   list of (y label, t, {label: values}). Curves with different
             units go in different groups, so each panel has one y axis.
    logy     logarithmic y axes (where the values are all positive)
    tmark    a time to mark with a vertical line
    Returns the list of axes made.
    """
    axes = []
    ng = max(len(groups), 1)
    gap = 0.05 if ng > 1 else 0.0
    height = (rect[3] - gap * (ng - 1)) / ng
    for g, (ylabel, t, curves) in enumerate(groups):
        bottom = rect[1] + (ng - 1 - g) * (height + gap)
        ax = fig.add_axes([rect[0], bottom, rect[2], height])
        axes.append(ax)
        ax.set_facecolor(theme['surface'])
        labels = list(curves)
        cols = (theme['series'][:len(labels)] if len(labels) <= 3
                else ordinal_colours(theme, len(labels)))
        for lab, col in zip(labels, cols):
            y = curves[lab]
            ax.plot(t, y, color=col, lw=2, label=lab, solid_capstyle='round')
            ok = np.isfinite(y)
            if 1 < len(labels) <= 4 and ok.any():
                # name the curve at its end, in text colour
                ax.annotate(' ' + lab, (t[ok][-1], y[ok][-1]), color=theme['text2'],
                            fontsize=9, va='center', annotation_clip=False)
        allv = np.concatenate([np.asarray(v, float) for v in curves.values()]) if curves else np.array([])
        allv = allv[np.isfinite(allv)]
        if logy and allv.size and (allv > 0).all():
            ax.set_yscale('log')
        if not allv.size:
            ax.text(0.5, 0.5, 'nothing to plot', transform=ax.transAxes,
                    ha='center', color=theme['text2'])
        if tmark is not None:
            ax.axvline(tmark, color=theme['text2'], lw=0.8, ls=(0, (4, 3)))
        ax.set_ylabel(ylabel, color=theme['text2'])
        ax.grid(True, color=theme['grid'], lw=0.6)
        ax.set_axisbelow(True)
        for side in ('top', 'right'):
            ax.spines[side].set_visible(False)
        for side in ('left', 'bottom'):
            ax.spines[side].set_color(theme['grid'])
        ax.tick_params(colors=theme['text2'], labelsize=9, which='both')
        if len(labels) > 1:
            leg = ax.legend(frameon=False, fontsize=9, labelcolor=theme['text'],
                            ncol=min(len(labels), 5), loc='lower left',
                            bbox_to_anchor=(0.0, 1.01), borderaxespad=0)
            leg.set_zorder(5)
        if g == ng - 1:
            ax.set_xlabel('time (Myr)', color=theme['text2'])
        else:
            ax.tick_params(labelbottom=False)
    if axes and title:
        axes[0].set_title(title, color=theme['text'], fontsize=12, loc='left', pad=44)
    return axes


# ----------------------------------------------------------------------
# the window
# ----------------------------------------------------------------------
VIEWS = ['3D view', 'Snapshot analysis', 'Run-time files']
OPTIONS = ['size by mass', 'shape by type', 'ring 10 most massive', 'fixed scale',
           'centre on c.o.m.', 'rotate when playing', 'dark background', 'log y axis (plots)']
OPTION_DEFAULTS = [True, False, False, True, True, False, False, False]

MAIN3D = [0.20, 0.09, 0.64, 0.88]
MAIN2D = [0.26, 0.15, 0.56, 0.72]
KEYBAR = [0.885, 0.36, 0.013, 0.40]
X0, W = 0.012, 0.166


class Viewer:
    def __init__(self, datadir='data', runid=None):
        self.datadir = datadir
        self.runs = find_runs(datadir)
        if not self.runs:
            raise SystemExit('No snapshots (<runid>_snap_<time>Myr.dat) found in %s' % datadir)
        self.runids = list(self.runs)
        if runid is None:
            runid = self.runids[0]
        if runid not in self.runs:
            raise SystemExit('No run %s in %s. Runs there: %s'
                             % (runid, datadir, ', '.join(self.runids)))
        self.view = VIEWS[0]
        self.colour = COLOUR_MODES[0]
        self.options = dict(zip(OPTIONS, OPTION_DEFAULTS))
        self.zoom = 1.0
        self.playing = False
        self.legend = None
        self.axes2d = []
        self.colwidgets = None
        self.fig = plt.figure(figsize=(15, 9))
        try:
            self.fig.canvas.manager.set_window_title('HENDOS viewer')
        except Exception:
            pass
        self.timer = self.fig.canvas.new_timer(interval=120)
        self.timer.add_callback(self.step)
        self.load_run(runid)
        self.build()
        self.redraw()

    # ---- the run ------------------------------------------------------
    def load_run(self, runid):
        self.run = Run(self.datadir, runid)
        self.isnap = 0
        self.series = {}
        self.quantity = list(ha.QUANTITIES)[0]
        self.afile = next(iter(self.run.analysis_files), None)
        self.acolumns = {}
        first = self.run.snapshot(0)
        self.rh0 = ha.half_mass_radius(first)

    @property
    def theme(self):
        return DARK if self.options['dark background'] else LIGHT

    # ---- building the window -----------------------------------------
    def panel(self, rect):
        ax = self.fig.add_axes(rect)
        ax.set_xticks([])
        ax.set_yticks([])
        for s in ax.spines.values():
            s.set_visible(False)
        return ax

    def box(self, rect, head=0.026):
        """A shaded box, and inside it (below room for a title) the axes for buttons."""
        back = self.panel(rect)
        back.set_facecolor(self.theme['panel'])
        inner = self.panel([rect[0], rect[1], rect[2], rect[3] - head])
        inner.set_facecolor(self.theme['panel'])
        return back, inner

    def build(self):
        fig, th = self.fig, self.theme
        # the old buttons must let go of the window before it is cleared
        for w in getattr(self, 'widgets', []):
            w.disconnect_events()
        if self.colwidgets is not None:
            self.colwidgets[1].disconnect_events()
        fig.clear()
        fig.patch.set_facecolor(th['surface'])
        self.legend = None
        self.axes2d = []
        self.colwidgets = None
        lab = dict(color=[th['text']], fontsize=[9.5])
        dot = dict(edgecolor=[th['text2']], s=[42])

        self.ax3d = fig.add_axes(MAIN3D, projection='3d')
        self.ax3d.view_init(elev=22, azim=-58)
        self.cax = fig.add_axes(KEYBAR)

        def titled(rect, title):
            _, ax = self.box(rect)
            fig.text(rect[0] + 0.006, rect[1] + rect[3] - 0.006, title, va='top',
                     color=th['text2'], fontsize=9, fontweight='bold')
            return ax

        # view
        ax = titled([X0, 0.842, W, 0.128], 'VIEW')
        self.w_view = RadioButtons(ax, VIEWS, active=VIEWS.index(self.view),
                                   activecolor=th['series'][0], label_props=lab,
                                   radio_props=dict(dot, s=[42]))
        self.w_view.on_clicked(self.set_view)

        # one of three selectors shares the next panel
        rect = [X0, 0.560, W, 0.270]
        self.t_select = fig.text(rect[0] + 0.006, rect[1] + rect[3] - 0.006, '', va='top',
                                 color=th['text2'], fontsize=9, fontweight='bold')
        self.bg_colour, self.ax_colour = self.box(rect)
        self.w_colour = RadioButtons(self.ax_colour, COLOUR_MODES,
                                     active=COLOUR_MODES.index(self.colour),
                                     activecolor=th['series'][0], label_props=lab,
                                     radio_props=dot)
        self.w_colour.on_clicked(self.set_colour)

        names = list(ha.QUANTITIES)
        self.bg_quant, self.ax_quant = self.box(rect)
        self.w_quant = RadioButtons(self.ax_quant, names, active=names.index(self.quantity),
                                    activecolor=th['series'][0],
                                    label_props=dict(lab, fontsize=[9]), radio_props=dot)
        self.w_quant.on_clicked(self.set_quantity)

        files = list(self.run.analysis_files) or ['(none for this run)']
        high = min(rect[3], 0.040 * (len(files) + 1) + 0.026)
        self.bg_file, self.ax_file = self.box([rect[0], rect[1] + rect[3] - high, rect[2], high])
        self.w_file = RadioButtons(self.ax_file, files,
                                   active=files.index(self.afile) if self.afile in files else 0,
                                   activecolor=th['series'][0], label_props=lab, radio_props=dot)
        self.w_file.on_clicked(self.set_file)

        # options
        ax = titled([X0, 0.262, W, 0.286], 'OPTIONS')
        self.w_opts = CheckButtons(ax, OPTIONS, [self.options[o] for o in OPTIONS],
                                   label_props=lab,
                                   frame_props=dict(edgecolor=[th['text2']], s=[44]),
                                   check_props=dict(facecolor=[th['series'][0]], s=[44]))
        self.w_opts.on_clicked(self.set_option)

        # buttons
        self.buttons = {}
        names = [('Play', self.toggle_play), ('Next run', self.next_run),
                 ('Zoom in', lambda e: self.set_zoom(1 / 1.4)),
                 ('Zoom out', lambda e: self.set_zoom(1.4)),
                 ('Save image', self.save_image), ('Save movie', self.save_movie)]
        for k, (name, func) in enumerate(names):
            bx = fig.add_axes([X0 + (k % 2) * (W / 2 + 0.002), 0.190 - (k // 2) * 0.062,
                               W / 2 - 0.002, 0.052])
            b = Button(bx, name, color=th['panel'], hovercolor=th['grid'])
            b.label.set_color(th['text'])
            b.label.set_fontsize(9.5)
            for s in bx.spines.values():
                s.set_edgecolor(th['grid'])
            b.on_clicked(func)
            self.buttons[name] = b

        # time slider
        sx = fig.add_axes([0.30, 0.035, 0.46, 0.026])
        sx.set_facecolor(th['panel'])
        nmax = max(self.run.nsnap - 1, 1)
        self.w_time = Slider(sx, 'snapshot ', 0, nmax, valinit=min(self.isnap, nmax),
                             valstep=1, color=th['series'][0], track_color=th['panel'],
                             handle_style=dict(facecolor=th['surface'],
                                               edgecolor=th['series'][0], size=12))
        self.w_time.label.set_color(th['text2'])
        self.w_time.valtext.set_visible(False)
        self.w_time.on_changed(self.set_snapshot)
        self.t_time = fig.text(0.775, 0.048, '', va='center', color=th['text'], fontsize=10)

        self.t_title = fig.text(0.21, 0.965, '', va='top', color=th['text'], fontsize=13)
        self.widgets = [self.w_view, self.w_colour, self.w_quant, self.w_file, self.w_opts,
                        self.w_time] + list(self.buttons.values())
        self.t_status = fig.text(X0, 0.012, '', color=th['text2'], fontsize=8.5)

    # ---- things the buttons do ---------------------------------------
    def set_view(self, label):
        self.view = label
        self.redraw()

    def set_colour(self, label):
        self.colour = label
        self.redraw()

    def set_quantity(self, label):
        self.quantity = label
        self.redraw()

    def set_file(self, label):
        if label in self.run.analysis_files:
            self.afile = label
        self.redraw()

    def set_option(self, label):
        self.options[label] = not self.options[label]
        if label == 'dark background':
            elev, azim = self.ax3d.elev, self.ax3d.azim
            self.build()                      # every colour changes
            self.ax3d.view_init(elev=elev, azim=azim)
        self.redraw()

    def set_snapshot(self, value):
        self.isnap = int(round(value))
        self.redraw()

    def set_zoom(self, factor):
        self.zoom *= factor
        self.redraw()

    def toggle_play(self, event=None):
        self.playing = not self.playing
        self.buttons['Play'].label.set_text('Pause' if self.playing else 'Play')
        if self.playing:
            self.timer.start()
        else:
            self.timer.stop()
        self.fig.canvas.draw_idle()

    def step(self):
        if self.options['rotate when playing'] and self.view == VIEWS[0]:
            self.ax3d.view_init(elev=self.ax3d.elev, azim=self.ax3d.azim + 2)
        self.w_time.set_val((self.isnap + 1) % self.run.nsnap)

    def next_run(self, event=None):
        k = (self.runids.index(self.run.runid) + 1) % len(self.runids)
        if self.playing:
            self.toggle_play()
        self.load_run(self.runids[k])
        elev, azim = self.ax3d.elev, self.ax3d.azim
        self.build()
        self.ax3d.view_init(elev=elev, azim=azim)
        self.redraw()

    def status(self, text):
        self.t_status.set_text(text)
        try:
            self.fig.canvas.draw_idle()
            self.fig.canvas.flush_events()
        except Exception:
            pass

    # ---- drawing ------------------------------------------------------
    def halfwidth(self):
        """Half the side of the box: fixed for the whole run, or None to fit each snapshot."""
        if self.options['fixed scale']:
            return 4.0 * self.rh0 * self.zoom
        return None

    def particle_options(self):
        o = self.options
        return dict(colour=self.colour, theme=self.theme, size_by_mass=o['size by mass'],
                    shape_by_type=o['shape by type'], ring_massive=o['ring 10 most massive'],
                    centre=o['centre on c.o.m.'])

    def draw_3d(self, ax, fig, cax, snap):
        hw = self.halfwidth()
        key = draw_particles(ax, snap, halfwidth=hw, **self.particle_options())
        if hw is None and self.zoom != 1.0:
            for get, put in ((ax.get_xlim, ax.set_xlim), (ax.get_ylim, ax.set_ylim),
                             (ax.get_zlim, ax.set_zlim)):
                lo, hi = get()
                mid, half = 0.5 * (lo + hi), 0.5 * (hi - lo) * self.zoom
                put(mid - half, mid + half)
        return draw_key(fig, cax, key, self.theme)

    def curve_groups(self):
        """What the two plot views show: (title, groups for plot_curves)."""
        if self.view == VIEWS[1]:
            name = self.quantity
            if name not in self.series:
                self.series[name] = ha.timeseries(
                    self.run, name,
                    progress=lambda i, n: (self.status('reading snapshots: %d of %d' % (i, n))
                                           if i % 5 == 0 or i == n else None))
                self.status('')
            t, curves = self.series[name]
            _, label, unit = ha.QUANTITIES[name]
            ylabel = label + (' (%s)' % unit if unit else '')
            if list(curves) == ['']:
                curves = {name: curves['']}
            return '%s   %s' % (self.run.runid, label), [(ylabel, t, curves)]
        if self.afile is None:
            return '%s   no run-time analysis files' % self.run.runid, []
        af = self.run.analysis(self.afile)
        chosen = self.acolumns.setdefault(self.afile, self.default_columns(af))
        t = af.data[:, 0]
        groups = {}
        for k, (col, unit) in enumerate(zip(af.columns, af.units)):
            if k > 0 and col in chosen:
                groups.setdefault(unit, {})[col] = af.data[:, k]
        return ('%s   run-time analysis: %s' % (self.run.runid, self.afile),
                [((unit if unit and unit != '-' else ''), t, c) for unit, c in groups.items()])

    @staticmethod
    def default_columns(af):
        """Start with the columns that share the unit of the first one."""
        cols = af.columns[1:]
        units = af.units[1:]
        return [c for c, u in zip(cols, units) if u == units[0]] if cols else []

    def column_boxes(self):
        """The tick boxes for the columns of the run-time file (right-hand side)."""
        if self.colwidgets is not None:
            self.colwidgets[1].disconnect_events()
            for artist in self.colwidgets[0]:
                artist.remove()
            self.colwidgets = None
        if self.view != VIEWS[2] or self.afile is None:
            return
        th = self.theme
        af = self.run.analysis(self.afile)
        cols = af.columns[1:]
        chosen = self.acolumns.setdefault(self.afile, self.default_columns(af))
        high = 0.040 * (len(cols) + 1) + 0.026
        rect = [0.865, 0.87 - high, 0.125, high]
        back, ax = self.box(rect)
        title = self.fig.text(rect[0] + 0.006, rect[1] + rect[3] - 0.006, 'COLUMNS', va='top',
                              color=th['text2'], fontsize=9, fontweight='bold')
        w = CheckButtons(ax, ['%s  (%s)' % (c, u) if u else c for c, u in zip(cols, af.units[1:])],
                         [c in chosen for c in cols],
                         label_props=dict(color=[th['text']], fontsize=[9]),
                         frame_props=dict(edgecolor=[th['text2']], s=[40]),
                         check_props=dict(facecolor=[th['series'][0]], s=[40]))

        def clicked(label):
            col = label.split('  (')[0]
            if col in chosen:
                chosen.remove(col)
            else:
                chosen.append(col)
            self.redraw()
        w.on_clicked(clicked)
        self.colwidgets = ((back, ax, title), w)

    def redraw(self):
        th = self.theme
        run = self.run
        three_d = self.view == VIEWS[0]
        # which selector is showing
        for axes, on in (((self.bg_colour, self.ax_colour), three_d),
                         ((self.bg_quant, self.ax_quant), self.view == VIEWS[1]),
                         ((self.bg_file, self.ax_file), self.view == VIEWS[2])):
            for a in axes:
                a.set_visible(on)
        self.t_select.set_text({VIEWS[0]: 'COLOUR BY', VIEWS[1]: 'QUANTITY',
                                VIEWS[2]: 'FILE'}[self.view])
        if self.legend is not None:
            self.legend.remove()
            self.legend = None
        for ax in self.axes2d:
            ax.remove()
        self.axes2d = []
        self.column_boxes()
        time = run.times[self.isnap]
        self.t_time.set_text('%.4f Myr   (%d of %d)' % (time, self.isnap + 1, run.nsnap))

        if three_d:
            snap = run.snapshot(self.isnap)
            elev, azim = self.ax3d.elev, self.ax3d.azim
            self.ax3d.cla()
            self.ax3d.set_visible(True)
            self.legend = self.draw_3d(self.ax3d, self.fig, self.cax, snap)
            self.ax3d.view_init(elev=elev, azim=azim)
            self.t_title.set_text('%s      t = %.4f Myr      %d particles      %.1f Msun'
                                  % (run.runid, snap.time, snap.n, snap.m.sum()))
        else:
            self.ax3d.set_visible(False)
            self.cax.set_visible(False)
            title, groups = self.curve_groups()
            self.t_title.set_text(title)
            if groups:
                self.axes2d = plot_curves(self.fig, MAIN2D, groups, th,
                                          logy=self.options['log y axis (plots)'], tmark=time)
            else:
                ax = self.panel(MAIN2D)
                ax.set_facecolor(th['surface'])
                ax.text(0.5, 0.5, 'Nothing to show: tick a column on the right'
                        if self.afile else 'This run has no run-time analysis files\n'
                        '(<runid>_ana_<name>.dat)', ha='center', va='center',
                        color=th['text2'], transform=ax.transAxes)
                self.axes2d = [ax]
        self.fig.canvas.draw_idle()

    # ---- saving -------------------------------------------------------
    def offscreen(self, snap):
        """A clean figure of one snapshot, with no buttons, for saving."""
        th = self.theme
        fig = Figure(figsize=(9, 8), facecolor=th['surface'])
        FigureCanvasAgg(fig)
        ax = fig.add_axes([0.0, 0.03, 0.82, 0.9], projection='3d')
        ax.view_init(elev=self.ax3d.elev, azim=self.ax3d.azim)
        cax = fig.add_axes([0.87, 0.34, 0.016, 0.40])
        legend = self.draw_3d(ax, fig, cax, snap)
        if legend is not None:
            legend.set_bbox_to_anchor((0.82, 0.62 if legend.get_title().get_text() != 'type'
                                       or cax.get_visible() is False else 0.26))
        title = fig.text(0.04, 0.955, '', color=th['text'], fontsize=13, va='top')
        return fig, ax, cax, title

    def outname(self, tail):
        return os.path.join(self.datadir, '%s_%s' % (self.run.runid, tail))

    def save_image(self, event=None):
        if self.view == VIEWS[0]:
            snap = self.run.snapshot(self.isnap)
            fig, ax, cax, title = self.offscreen(snap)
            title.set_text('%s      t = %.4f Myr' % (self.run.runid, snap.time))
            name = self.outname('view_%09.4fMyr.png' % snap.time)
        else:
            th = self.theme
            fig = Figure(figsize=(9, 6), facecolor=th['surface'])
            FigureCanvasAgg(fig)
            title, groups = self.curve_groups()
            plot_curves(fig, [0.11, 0.11, 0.78, 0.74], groups, th,
                        logy=self.options['log y axis (plots)'], title=title)
            tag = self.quantity if self.view == VIEWS[1] else 'ana_' + str(self.afile)
            name = self.outname('plot_%s.png' % tag.replace(' ', '_').replace('/', ''))
        fig.savefig(name, dpi=150, facecolor=fig.get_facecolor())
        self.status('saved ' + name)
        return name

    def save_movie(self, event=None, fps=10):
        """Animation of the 3D view through every snapshot, as shown now."""
        if self.playing:
            self.toggle_play()
        run = self.run
        fig, ax, cax, title = self.offscreen(run.snapshot(0))
        azim0, rotate = self.ax3d.azim, self.options['rotate when playing']
        state = {'legend': None}

        def frame(i):
            snap = run.snapshot(i)
            ax.cla()
            for old in list(fig.legends):
                old.remove()
            self.draw_3d(ax, fig, cax, snap)
            for leg in fig.legends:
                leg.set_bbox_to_anchor((0.82, 0.62 if not cax.get_visible() else 0.26))
            ax.view_init(elev=self.ax3d.elev, azim=azim0 + (2 * i if rotate else 0))
            title.set_text('%s      t = %.4f Myr' % (run.runid, snap.time))
            if i % 5 == 0:
                self.status('making movie: frame %d of %d' % (i + 1, run.nsnap))
            return []

        anim = animation.FuncAnimation(fig, frame, frames=run.nsnap, blit=False)
        if animation.writers.is_available('ffmpeg'):
            name = self.outname('movie.mp4')
            anim.save(name, writer=animation.FFMpegWriter(fps=fps, bitrate=2400), dpi=110,
                      savefig_kwargs=dict(facecolor=fig.get_facecolor()))
        else:
            name = self.outname('movie.gif')
            anim.save(name, writer=animation.PillowWriter(fps=fps), dpi=90,
                      savefig_kwargs=dict(facecolor=fig.get_facecolor()))
        self.status('saved ' + name)
        return name


def main():
    args = sys.argv[1:]
    datadir = args[0] if len(args) > 0 else 'data'
    runid = args[1] if len(args) > 1 else None
    Viewer(datadir, runid)
    plt.show()


if __name__ == '__main__':
    main()
